use flate2::{Compression, write::GzEncoder};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
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

    xloci(Args {
        sequence: sequence_path,
        regions: regions_path,
        outdir: outdir.clone(),
        chunks: 1,
        upstream_flank: 0,
        downstream_flank: 0,
        feature: Feature::Exon,
        ignore_errors: false,
        level: log::Level::Info,
        prefix: "output.fa".to_string(),
        translate: false,
        as_chunk: false,
        include_bed: false,
        compress: false,
    });

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
