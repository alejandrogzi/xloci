use flate2::{Compression, write::GzEncoder};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};
use tempfile::TempDir;
use twobit::convert::{fasta::FastaReader, to_2bit};
use xloci::{Args, Feature, xloci};

const GENOME_FASTA: &str = ">chr1\nAACCGGTTTACGATCG\n";
const BED_CONTENT: &str = "chr1\t0\t4\tbed_plus\t0\t+\t0\t4\t0,0,0\t1\t4\t0\nchr1\t8\t12\tbed_minus\t0\t-\t8\t12\t0,0,0\t1\t4\t0\n";
const GTF_CONTENT: &str = "chr1\tsrc\ttranscript\t1\t4\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx_plus\";\nchr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx_plus\";\nchr1\tsrc\ttranscript\t9\t12\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_minus\";\nchr1\tsrc\texon\t9\t12\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx_minus\";\n";
const GFF_CONTENT: &str = "##gff-version 3\nchr1\tsrc\tmRNA\t1\t4\t.\t+\t.\tID=tx_plus;Name=tx_plus\nchr1\tsrc\texon\t1\t4\t.\t+\t.\tParent=tx_plus\nchr1\tsrc\tmRNA\t9\t12\t.\t-\t.\tID=tx_minus;Name=tx_minus\nchr1\tsrc\texon\t9\t12\t.\t-\t.\tParent=tx_minus\n";

#[derive(Clone, Copy)]
enum SequenceFormat {
    Fa,
    TwoBit,
}

#[derive(Clone, Copy)]
enum RegionFormat {
    Bed,
    Gtf,
    Gff,
}

struct Case {
    sequence_format: SequenceFormat,
    region_format: RegionFormat,
    region_gz: bool,
}

fn base_args(sequence: PathBuf, regions: PathBuf, outdir: PathBuf) -> Args {
    Args {
        sequence,
        regions,
        outdir,
        chunks: 1,
        upstream_flank: 0,
        downstream_flank: 0,
        feature: Feature::Exon,
        ignore_errors: false,
        level: log::Level::Info,
        prefix: "output".to_string(),
        translate: false,
        split_extraction: false,
        as_tsv: false,
        add_tab: false,
        generic_id: false,
        as_chunk: false,
        include_bed: false,
        compress: false,
        threads: 1,
    }
}

fn run_case(case: Case) {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();

    let fasta_path = root.join("genome.fa");
    write_bytes(&fasta_path, GENOME_FASTA.as_bytes());

    let twobit_path = root.join("genome.2bit");
    write_twobit(&fasta_path, &twobit_path);

    let sequence_path = match case.sequence_format {
        SequenceFormat::Fa => fasta_path,
        SequenceFormat::TwoBit => twobit_path,
    };

    let regions_path = write_regions(root, case.region_format, case.region_gz);
    let outdir = root.join("out");

    xloci(base_args(sequence_path, regions_path, outdir.clone()));

    let records = read_fasta(outdir.join("output.fa"));
    let (plus_name, minus_name) = expected_names(case.region_format);

    assert_eq!(
        records.get(plus_name).map(std::string::String::as_str),
        Some("AACC"),
        "plus-strand sequence mismatch for {}",
        case_name(&case)
    );
    assert_eq!(
        records.get(minus_name).map(std::string::String::as_str),
        Some("CGTA"),
        "minus-strand sequence mismatch for {}",
        case_name(&case)
    );
    assert_eq!(
        records.len(),
        2,
        "unexpected record count for {}",
        case_name(&case)
    );
}

fn run_stdin_case(case: Case) {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();

    let fasta_path = root.join("genome.fa");
    write_bytes(&fasta_path, GENOME_FASTA.as_bytes());

    let twobit_path = root.join("genome.2bit");
    write_twobit(&fasta_path, &twobit_path);

    let stdin_bytes = match case.sequence_format {
        SequenceFormat::Fa => GENOME_FASTA.as_bytes().to_vec(),
        SequenceFormat::TwoBit => std::fs::read(&twobit_path)
            .unwrap_or_else(|e| panic!("failed to read {}: {}", twobit_path.display(), e)),
    };

    let regions_path = write_regions(root, case.region_format, case.region_gz);
    let outdir = root.join("out");
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_xloci"));

    let mut child = Command::new(binary)
        .arg("-r")
        .arg(&regions_path)
        .arg("-o")
        .arg(&outdir)
        .arg("-c")
        .arg("1")
        .arg("-L")
        .arg("error")
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn xloci binary");

    let mut stdin = child.stdin.take().expect("failed to open child stdin");
    stdin
        .write_all(&stdin_bytes)
        .unwrap_or_else(|e| panic!("failed to write stdin for {}: {}", case_name(&case), e));
    drop(stdin);

    let output = child
        .wait_with_output()
        .expect("failed to wait for xloci binary");

    assert!(
        output.status.success(),
        "stdin case {} failed: {}",
        case_name(&case),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = read_fasta(outdir.join("output.fa"));
    let (plus_name, minus_name) = expected_names(case.region_format);

    assert_eq!(
        records.get(plus_name).map(std::string::String::as_str),
        Some("AACC"),
        "plus-strand sequence mismatch for stdin {}",
        case_name(&case)
    );
    assert_eq!(
        records.get(minus_name).map(std::string::String::as_str),
        Some("CGTA"),
        "minus-strand sequence mismatch for stdin {}",
        case_name(&case)
    );
    assert_eq!(
        records.len(),
        2,
        "unexpected record count for stdin {}",
        case_name(&case)
    );
}

fn expected_names(format: RegionFormat) -> (&'static str, &'static str) {
    match format {
        RegionFormat::Bed => ("bed_plus", "bed_minus"),
        RegionFormat::Gtf | RegionFormat::Gff => ("tx_plus", "tx_minus"),
    }
}

fn case_name(case: &Case) -> String {
    let seq = match case.sequence_format {
        SequenceFormat::Fa => "fa",
        SequenceFormat::TwoBit => "2bit",
    };
    let region = match case.region_format {
        RegionFormat::Bed => "bed",
        RegionFormat::Gtf => "gtf",
        RegionFormat::Gff => "gff",
    };
    let suffix = if case.region_gz { ".gz" } else { "" };

    format!("{seq}/{region}{suffix}")
}

fn write_regions(root: &Path, format: RegionFormat, gz: bool) -> PathBuf {
    let (stem, content) = match format {
        RegionFormat::Bed => ("regions.bed", BED_CONTENT),
        RegionFormat::Gtf => ("regions.gtf", GTF_CONTENT),
        RegionFormat::Gff => ("regions.gff", GFF_CONTENT),
    };

    let plain_path = root.join(stem);
    write_bytes(&plain_path, content.as_bytes());

    if !gz {
        return plain_path;
    }

    let gz_path = root.join(format!("{stem}.gz"));
    write_gzip(&gz_path, content.as_bytes());
    gz_path
}

fn write_bytes(path: &Path, content: &[u8]) {
    std::fs::write(path, content)
        .unwrap_or_else(|e| panic!("failed to write {}: {}", path.display(), e));
}

fn write_gzip(path: &Path, content: &[u8]) {
    let file = File::create(path)
        .unwrap_or_else(|e| panic!("failed to create gzip file {}: {}", path.display(), e));
    let mut writer = GzEncoder::new(BufWriter::new(file), Compression::default());

    writer
        .write_all(content)
        .unwrap_or_else(|e| panic!("failed to write gzip content {}: {}", path.display(), e));
    writer
        .finish()
        .unwrap_or_else(|e| panic!("failed to finish gzip file {}: {}", path.display(), e));
}

fn write_twobit(fasta_path: &Path, twobit_path: &Path) {
    let reader = FastaReader::open(fasta_path)
        .unwrap_or_else(|e| panic!("failed to open FASTA {}: {}", fasta_path.display(), e));
    let file = File::create(twobit_path)
        .unwrap_or_else(|e| panic!("failed to create {}: {}", twobit_path.display(), e));
    let mut writer = BufWriter::new(file);

    to_2bit(&mut writer, &reader)
        .unwrap_or_else(|e| panic!("failed to write 2bit {}: {}", twobit_path.display(), e));
    writer
        .flush()
        .unwrap_or_else(|e| panic!("failed to flush {}: {}", twobit_path.display(), e));
}

fn read_fasta(path: PathBuf) -> HashMap<String, String> {
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {}", path.display(), e));

    let mut records = HashMap::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(prev_header) = header.replace(rest.to_string()) {
                records.insert(prev_header, std::mem::take(&mut seq));
            }
        } else {
            seq.push_str(line);
        }
    }

    if let Some(last_header) = header {
        records.insert(last_header, seq);
    }

    records
}

fn read_tsv(path: PathBuf) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {}", path.display(), e));

    text.lines()
        .map(|line| line.split('\t').map(|field| field.to_string()).collect())
        .collect()
}

fn rows_to_map(rows: Vec<Vec<String>>) -> HashMap<String, Vec<String>> {
    rows.into_iter()
        .map(|row| (row[0].clone(), row))
        .collect::<HashMap<_, _>>()
}

#[test]
fn test_2bit_bed() {
    run_case(Case {
        sequence_format: SequenceFormat::TwoBit,
        region_format: RegionFormat::Bed,
        region_gz: false,
    });
}

#[test]
fn test_2bit_gtf() {
    run_case(Case {
        sequence_format: SequenceFormat::TwoBit,
        region_format: RegionFormat::Gtf,
        region_gz: false,
    });
}

#[test]
fn test_2bit_gff_gz() {
    run_case(Case {
        sequence_format: SequenceFormat::TwoBit,
        region_format: RegionFormat::Gff,
        region_gz: true,
    });
}

#[test]
fn test_fa_bed() {
    run_case(Case {
        sequence_format: SequenceFormat::Fa,
        region_format: RegionFormat::Bed,
        region_gz: false,
    });
}

#[test]
fn test_fa_gff() {
    run_case(Case {
        sequence_format: SequenceFormat::Fa,
        region_format: RegionFormat::Gff,
        region_gz: false,
    });
}

#[test]
fn test_fa_gtf_gz() {
    run_case(Case {
        sequence_format: SequenceFormat::Fa,
        region_format: RegionFormat::Gtf,
        region_gz: true,
    });
}

#[test]
fn test_2bit_bed_gz() {
    run_case(Case {
        sequence_format: SequenceFormat::TwoBit,
        region_format: RegionFormat::Bed,
        region_gz: true,
    });
}

#[test]
fn test_fa_bed_gz() {
    run_case(Case {
        sequence_format: SequenceFormat::Fa,
        region_format: RegionFormat::Bed,
        region_gz: true,
    });
}

#[test]
fn test_stdin_fa_bed() {
    run_stdin_case(Case {
        sequence_format: SequenceFormat::Fa,
        region_format: RegionFormat::Bed,
        region_gz: false,
    });
}

#[test]
fn test_stdin_2bit_gff_gz() {
    run_stdin_case(Case {
        sequence_format: SequenceFormat::TwoBit,
        region_format: RegionFormat::Gff,
        region_gz: true,
    });
}

#[test]
fn test_split_extraction_minus_strand_uses_transcript_order() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let outdir = root.join("out");

    write_bytes(&genome_path, b">chr1\nAACCGGTTTACGATCG\n");
    write_bytes(
        &regions_path,
        b"chr1\t0\t12\ttx_minus\t0\t-\t0\t12\t0,0,0\t2\t4,4\t0,8\n",
    );

    let mut args = base_args(genome_path, regions_path, outdir.clone());
    args.split_extraction = true;
    xloci(args);

    let records = read_fasta(outdir.join("output.fa"));
    assert_eq!(
        records
            .get("tx_minus_EXON1")
            .map(std::string::String::as_str),
        Some("CGTA")
    );
    assert_eq!(
        records
            .get("tx_minus_EXON2")
            .map(std::string::String::as_str),
        Some("GGTT")
    );
    assert_eq!(records.len(), 2);
}

#[test]
fn test_split_extraction_transcript_uses_mixed_piece_names() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let outdir = root.join("out");

    write_bytes(&genome_path, b">chr1\nAACCGGTTTACGATCG\n");
    write_bytes(
        &regions_path,
        b"chr1\t0\t12\ttx_plus\t0\t+\t0\t12\t0,0,0\t2\t4,4\t0,8\n",
    );

    let mut args = base_args(genome_path, regions_path, outdir.clone());
    args.feature = Feature::Transcript;
    args.split_extraction = true;
    xloci(args);

    let records = read_fasta(outdir.join("output.fa"));
    assert_eq!(
        records
            .get("tx_plus_EXON1")
            .map(std::string::String::as_str),
        Some("AACC")
    );
    assert_eq!(
        records
            .get("tx_plus_INTRON1")
            .map(std::string::String::as_str),
        Some("GGTT")
    );
    assert_eq!(
        records
            .get("tx_plus_EXON2")
            .map(std::string::String::as_str),
        Some("TACG")
    );
    assert_eq!(records.len(), 3);
}

#[test]
fn test_split_extraction_cds_and_utr_piece_names() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let cds_outdir = root.join("cds_out");
    let utr_outdir = root.join("utr_out");

    write_bytes(&genome_path, b">chr1\nACGTACGTACGTACGTACGT\n");
    write_bytes(
        &regions_path,
        b"chr1\t0\t12\ttx_plus\t0\t+\t3\t9\t0,0,0\t2\t4,4\t0,8\n",
    );

    let mut cds_args = base_args(
        genome_path.clone(),
        regions_path.clone(),
        cds_outdir.clone(),
    );
    cds_args.feature = Feature::CDS;
    cds_args.split_extraction = true;
    xloci(cds_args);

    let cds_records = read_fasta(cds_outdir.join("output.fa"));
    assert_eq!(
        cds_records
            .get("tx_plus_CDS1")
            .map(std::string::String::as_str),
        Some("T")
    );
    assert_eq!(
        cds_records
            .get("tx_plus_CDS2")
            .map(std::string::String::as_str),
        Some("A")
    );

    let mut utr_args = base_args(genome_path, regions_path, utr_outdir.clone());
    utr_args.feature = Feature::UTR;
    utr_args.split_extraction = true;
    xloci(utr_args);

    let utr_records = read_fasta(utr_outdir.join("output.fa"));
    assert_eq!(
        utr_records
            .get("tx_plus_UTR1")
            .map(std::string::String::as_str),
        Some("ACG")
    );
    assert_eq!(
        utr_records
            .get("tx_plus_UTR2")
            .map(std::string::String::as_str),
        Some("CGT")
    );
}

#[test]
fn test_as_tsv_add_tab_with_both_flanks() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let outdir = root.join("out");

    write_bytes(&genome_path, b">chr1\nAACCGGTTTACGATCG\n");
    write_bytes(
        &regions_path,
        b"chr1\t2\t10\ttx_plus\t0\t+\t2\t10\t0,0,0\t2\t2,2\t0,6\nchr1\t2\t10\ttx_minus\t0\t-\t2\t10\t0,0,0\t2\t2,2\t0,6\n",
    );

    let mut args = base_args(genome_path, regions_path, outdir.clone());
    args.as_tsv = true;
    args.add_tab = true;
    args.upstream_flank = 1;
    args.downstream_flank = 2;
    xloci(args);

    let rows = rows_to_map(read_tsv(outdir.join("output.tsv")));
    assert_eq!(
        rows.get("tx_plus"),
        Some(&vec![
            "tx_plus".to_string(),
            "A".to_string(),
            "CCTA".to_string(),
            "CG".to_string()
        ])
    );
    assert_eq!(
        rows.get("tx_minus"),
        Some(&vec![
            "tx_minus".to_string(),
            "CG".to_string(),
            "TAGG".to_string(),
            "T".to_string()
        ])
    );
}

#[test]
fn test_as_tsv_add_tab_with_single_flank_uses_three_columns() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let outdir = root.join("out");

    write_bytes(&genome_path, b">chr1\nAACCGGTTTACGATCG\n");
    write_bytes(
        &regions_path,
        b"chr1\t2\t10\ttx_plus\t0\t+\t2\t10\t0,0,0\t2\t2,2\t0,6\n",
    );

    let mut args = base_args(genome_path, regions_path, outdir.clone());
    args.as_tsv = true;
    args.add_tab = true;
    args.downstream_flank = 2;
    xloci(args);

    let rows = read_tsv(outdir.join("output.tsv"));
    assert_eq!(
        rows,
        vec![vec![
            "tx_plus".to_string(),
            "CG".to_string(),
            "CCTA".to_string()
        ]]
    );
}

#[test]
fn test_split_extraction_as_tsv_add_tab_uses_piece_flanks() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let outdir = root.join("out");

    write_bytes(&genome_path, b">chr1\nAACCGGTTTACGATCG\n");
    write_bytes(
        &regions_path,
        b"chr1\t2\t10\ttx_minus\t0\t-\t2\t10\t0,0,0\t2\t2,2\t0,6\n",
    );

    let mut args = base_args(genome_path, regions_path, outdir.clone());
    args.as_tsv = true;
    args.add_tab = true;
    args.split_extraction = true;
    args.upstream_flank = 1;
    args.downstream_flank = 2;
    xloci(args);

    let rows = rows_to_map(read_tsv(outdir.join("output.tsv")));
    assert_eq!(
        rows.get("tx_minus_EXON1"),
        Some(&vec![
            "tx_minus_EXON1".to_string(),
            "CG".to_string(),
            "TA".to_string(),
            "A".to_string()
        ])
    );
    assert_eq!(
        rows.get("tx_minus_EXON2"),
        Some(&vec![
            "tx_minus_EXON2".to_string(),
            "CC".to_string(),
            "GG".to_string(),
            "T".to_string()
        ])
    );
}

#[test]
fn test_split_extraction_is_incompatible_with_translate() {
    let result = std::panic::catch_unwind(|| {
        xloci(Args {
            translate: true,
            split_extraction: true,
            ..base_args(
                PathBuf::from("unused.fa"),
                PathBuf::from("unused.bed"),
                PathBuf::from("unused"),
            )
        });
    });

    assert!(result.is_err());
}

#[test]
fn test_add_tab_requires_tsv_and_flanks() {
    let missing_tsv = std::panic::catch_unwind(|| {
        xloci(Args {
            add_tab: true,
            upstream_flank: 1,
            ..base_args(
                PathBuf::from("unused.fa"),
                PathBuf::from("unused.bed"),
                PathBuf::from("unused"),
            )
        });
    });
    assert!(missing_tsv.is_err());

    let missing_flanks = std::panic::catch_unwind(|| {
        xloci(Args {
            add_tab: true,
            as_tsv: true,
            ..base_args(
                PathBuf::from("unused.fa"),
                PathBuf::from("unused.bed"),
                PathBuf::from("unused"),
            )
        });
    });
    assert!(missing_flanks.is_err());
}

#[test]
fn test_generic_id_joined_uses_feature_bounds_without_flanks() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let outdir = root.join("out");

    write_bytes(&genome_path, b">chr1\nAACCGGTTTACGATCG\n");
    write_bytes(
        &regions_path,
        b"chr1\t2\t10\ttx_minus\t0\t-\t2\t10\t0,0,0\t2\t2,2\t0,6\n",
    );

    let mut args = base_args(genome_path, regions_path, outdir.clone());
    args.generic_id = true;
    args.upstream_flank = 1;
    args.downstream_flank = 2;
    xloci(args);

    let records = read_fasta(outdir.join("output.fa"));
    assert_eq!(
        records.get("chr1:2-10(-)").map(std::string::String::as_str),
        Some("CGTAGGT")
    );
}

#[test]
fn test_generic_id_split_uses_piece_bounds_without_flanks() {
    let temp = TempDir::new().expect("failed to create temporary directory");
    let root = temp.path();
    let genome_path = root.join("genome.fa");
    let regions_path = root.join("regions.bed");
    let outdir = root.join("out");

    write_bytes(&genome_path, b">chr1\nAACCGGTTTACGATCG\n");
    write_bytes(
        &regions_path,
        b"chr1\t2\t10\ttx_minus\t0\t-\t2\t10\t0,0,0\t2\t2,2\t0,6\n",
    );

    let mut args = base_args(genome_path, regions_path, outdir.clone());
    args.generic_id = true;
    args.split_extraction = true;
    args.upstream_flank = 1;
    args.downstream_flank = 2;
    xloci(args);

    let records = read_fasta(outdir.join("output.fa"));
    assert_eq!(
        records.get("chr1:8-10(-)").map(std::string::String::as_str),
        Some("CGTAA")
    );
    assert_eq!(
        records.get("chr1:2-4(-)").map(std::string::String::as_str),
        Some("CCGGT")
    );
}
