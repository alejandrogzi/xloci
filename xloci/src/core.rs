use crate::{
    cli::{Args, Feature},
    consts::CODON_TABLE,
};

use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use genepred::{Bed12, GenePred, Gff, Gtf, Reader, ReaderResult, Strand, Writer, bed::BedFormat};
use log::{info, warn};
use rayon::prelude::*;
use twobit::TwoBitFile;

use std::{
    collections::HashMap,
    fmt::Debug,
    fs::{File, create_dir_all},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

/// Main processing function that orchestrates genomic sequence extraction.
pub fn xloci(args: Args) {
    let Args {
        sequence,
        regions,
        outdir,
        chunks,
        upstream_flank,
        downstream_flank,
        feature,
        ignore_errors,
        prefix,
        translate,
        as_chunk,
        include_bed,
        compress,
        ..
    } = args;

    if include_bed && !as_chunk {
        panic!("ERROR: --include-bed requires --as-chunk");
    }

    let genome = get_sequences(sequence);
    create_dir_all(&outdir).unwrap_or_else(|e| panic!("{}", e));

    match detect_region_format(&regions) {
        Some(RegionFormat::Bed) => process_reader::<Bed12>(
            &regions,
            chunks,
            &outdir,
            &genome,
            upstream_flank,
            downstream_flank,
            feature,
            ignore_errors,
            &prefix,
            translate,
            as_chunk,
            include_bed,
            compress,
        ),
        Some(RegionFormat::Gtf) => process_reader::<Gtf>(
            &regions,
            chunks,
            &outdir,
            &genome,
            upstream_flank,
            downstream_flank,
            feature,
            ignore_errors,
            &prefix,
            translate,
            as_chunk,
            include_bed,
            compress,
        ),
        Some(RegionFormat::Gff) => process_reader::<Gff>(
            &regions,
            chunks,
            &outdir,
            &genome,
            upstream_flank,
            downstream_flank,
            feature,
            ignore_errors,
            &prefix,
            translate,
            as_chunk,
            include_bed,
            compress,
        ),
        None => panic!("ERROR: Unsupported file format"),
    }
}

/// Processes genomic regions in parallel chunks and writes output to FASTA file.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
fn process_reader<R>(
    regions: &Path,
    chunks: usize,
    outdir: &Path,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    upstream_flank: usize,
    downstream_flank: usize,
    feature_type: Feature,
    ignore_errors: bool,
    prefix: &str,
    translate: bool,
    as_chunk: bool,
    include_bed: bool,
    compress: bool,
) where
    R: BedFormat + Into<GenePred> + Send,
{
    info!("Processing regions from file {}", regions.display());

    let tmp_dir = outdir.join("tmp");
    create_dir_all(&tmp_dir).unwrap_or_else(|e| panic!("{}", e));

    let collector: Option<Arc<Mutex<Vec<(usize, PathBuf)>>>> = if as_chunk {
        None
    } else {
        Some(Arc::new(Mutex::new(Vec::new())))
    };

    read_regions::<R>(regions)
        .unwrap_or_else(|e| panic!("{}", e))
        .par_chunks(chunks)
        .unwrap_or_else(|e| panic!("{}", e))
        .for_each(|(idx, chunk)| {
            write_chunk(
                idx,
                chunk,
                genome,
                upstream_flank,
                downstream_flank,
                &feature_type,
                ignore_errors,
                translate,
                collector.clone(),
                &tmp_dir,
                include_bed,
                compress && as_chunk,
            );
        });

    if as_chunk {
        info!("Wrote chunk outputs to {}", tmp_dir.display());
        return;
    }

    let mut chunk_paths = collector
        .unwrap_or_else(|| panic!("ERROR: missing collector"))
        .lock()
        .unwrap_or_else(|e| panic!("ERROR: Cannot acquire lock on collector: {}", e))
        .clone();

    chunk_paths.sort_by_key(|(idx, _)| *idx);

    let output_path = with_gzip_extension(outdir.join(prefix), compress);
    let output_file = File::create(&output_path)
        .unwrap_or_else(|e| panic!("ERROR: cannot create {}: {}", output_path.display(), e));

    if compress {
        let mut writer = GzEncoder::new(BufWriter::new(output_file), Compression::default());
        merge_and_cleanup_chunks(&chunk_paths, &mut writer);
        writer
            .finish()
            .unwrap_or_else(|e| panic!("ERROR: cannot finish gzip stream: {}", e));
    } else {
        let mut writer = BufWriter::new(output_file);
        merge_and_cleanup_chunks(&chunk_paths, &mut writer);
        writer
            .flush()
            .unwrap_or_else(|e| panic!("ERROR: cannot flush writer: {}", e));
    }
}

fn read_regions<R>(regions: &Path) -> genepred::ReaderResult<Reader<R>>
where
    R: BedFormat + Into<GenePred> + Send,
{
    if is_compressed_path(regions) {
        Reader::<R>::from_path(regions)
    } else {
        Reader::<R>::from_mmap(regions)
    }
}

fn is_compressed_path(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("gz" | "zst" | "zstd" | "bz2" | "bzip2")
    )
}

fn merge_and_cleanup_chunks<W: Write>(chunk_paths: &[(usize, PathBuf)], writer: &mut W) {
    for (_, path) in chunk_paths {
        let mut chunk_file = File::open(path)
            .unwrap_or_else(|e| panic!("ERROR: cannot open chunk {}: {}", path.display(), e));
        std::io::copy(&mut chunk_file, writer).unwrap_or_else(|e| panic!("{}", e));
        std::fs::remove_file(path)
            .unwrap_or_else(|e| panic!("ERROR: cannot remove chunk {}: {}", path.display(), e));
    }
}

fn with_gzip_extension(mut path: PathBuf, compress: bool) -> PathBuf {
    if compress && path.extension().and_then(|ext| ext.to_str()) != Some("gz") {
        path.as_mut_os_string().push(".gz");
    }

    path
}

/// Processes a chunk of genomic records and writes extracted sequences to a temporary file.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
fn write_chunk(
    idx: usize,
    chunk: Vec<ReaderResult<GenePred>>,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    upstream_flank: usize,
    downstream_flank: usize,
    feature_type: &Feature,
    ignore_errors: bool,
    to_protein: bool,
    collector: Option<Arc<Mutex<Vec<(usize, PathBuf)>>>>,
    outdir: &Path,
    include_bed: bool,
    compress: bool,
) {
    info!("Processing chunk {}", idx);

    let chunk_base = outdir.join(format!("tmp_{}", idx));
    let fasta_path = if compress {
        chunk_base.with_extension("fa.gz")
    } else {
        chunk_base.with_extension("fa")
    };

    let fasta_file = File::create(&fasta_path)
        .unwrap_or_else(|e| panic!("ERROR: cannot create {}: {}", fasta_path.display(), e));
    let mut writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(
            BufWriter::new(fasta_file),
            Compression::default(),
        ))
    } else {
        Box::new(BufWriter::new(fasta_file))
    };

    let mut bed_writer = if include_bed {
        let bed_path = chunk_base.with_extension("bed");
        Some(BufWriter::new(File::create(&bed_path).unwrap_or_else(
            |e| panic!("ERROR: cannot create {}: {}", bed_path.display(), e),
        )))
    } else {
        None
    };

    chunk
        .into_iter()
        .filter_map(|result| match result {
            Ok(record) => Some(record),
            Err(e) => {
                if ignore_errors {
                    eprintln!("WARN: Failed to process record: {}", e);
                    None
                } else {
                    panic!("ERROR: Failed to process record: {}", e);
                }
            }
        })
        .for_each(|record| {
            let seq = genome.get(&record.chrom).unwrap_or_else(|| {
                panic!(
                    "ERROR: Chromosome {} not found!",
                    String::from_utf8_lossy(&record.chrom)
                )
            });

            let target = extract_seq(
                &record,
                seq,
                upstream_flank,
                downstream_flank,
                feature_type,
                ignore_errors,
            );

            if let Some(mut target) = target {
                match &record.strand {
                    Some(Strand::Forward) => {}
                    Some(Strand::Reverse) => {
                        target.reverse();

                        for base in target.iter_mut() {
                            *base = match *base {
                                b'A' => b'T',
                                b'C' => b'G',
                                b'G' => b'C',
                                b'T' => b'A',
                                b'N' => b'N',
                                b'a' => b't',
                                b'c' => b'g',
                                b'g' => b'c',
                                b't' => b'a',
                                b'n' => b'n',
                                _ => panic!("ERROR: Invalid base"),
                            }
                        }
                    }
                    Some(Strand::Unknown) | None => {}
                }

                if to_protein {
                    target = translate(&target);
                }

                if target.is_empty() {
                    warn!("WARN: empty sequence for record {}", record);
                    return;
                }

                writer.write_all(b">").unwrap_or_else(|e| panic!("{}", e));
                writer
                    .write_all(record.name().unwrap())
                    .unwrap_or_else(|e| panic!("{}", e));
                writer.write_all(b"\n").unwrap_or_else(|e| panic!("{}", e));
                writer
                    .write_all(&target)
                    .unwrap_or_else(|e| panic!("{}", e));
                writer.write_all(b"\n").unwrap_or_else(|e| panic!("{}", e));

                if let Some(bed_writer) = &mut bed_writer {
                    Writer::<Bed12>::from_record(&record, bed_writer)
                        .unwrap_or_else(|e| panic!("{}", e));
                }
            }
        });

    writer
        .flush()
        .unwrap_or_else(|e| panic!("ERROR: cannot flush {}: {}", fasta_path.display(), e));

    if let Some(bed_writer) = &mut bed_writer {
        bed_writer
            .flush()
            .unwrap_or_else(|e| panic!("ERROR: cannot flush BED writer: {}", e));
    }

    if let Some(collector) = collector {
        collector
            .lock()
            .unwrap_or_else(|e| panic!("ERROR: Cannot acquire lock on collector: {}", e))
            .push((idx, fasta_path));
    }
}

/// Error type for range calculation failures during sequence extraction.
#[derive(Debug)]
#[allow(dead_code)]
enum RangeError {
    Underflow { feature_coord: usize, flank: usize },
}

impl std::fmt::Display for RangeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RangeError::Underflow {
                feature_coord,
                flank,
            } => write!(
                f,
                "ERROR: Feature coordinate {} is underflowing by {} bases",
                feature_coord, flank
            ),
        }
    }
}

/// Calculates the slice range for an exon with appropriate flanking regions.
fn slice_range_for_exon(
    exon_idx: usize,
    exon_count: usize,
    feature_start: usize,
    feature_end: usize,
    upstream_flank: usize,
    downstream_flank: usize,
) -> Result<std::ops::Range<usize>, RangeError> {
    let is_first = exon_idx == 0;
    let is_last = exon_idx + 1 == exon_count;

    let start = if is_first {
        feature_start
            .checked_sub(upstream_flank)
            .ok_or(RangeError::Underflow {
                feature_coord: feature_start,
                flank: upstream_flank,
            })?
    } else {
        feature_start
    };

    let end = if is_last {
        feature_end
            .checked_add(downstream_flank)
            .ok_or(RangeError::Underflow {
                feature_coord: feature_end,
                flank: downstream_flank,
            })?
    } else {
        feature_end
    };

    Ok(start..end)
}

/// Extends target sequence with a slice or handles out-of-bounds errors gracefully.
fn extend_or_handle_oob(
    record: &GenePred,
    target: &mut Vec<u8>,
    seq: &[u8],
    range: std::ops::Range<usize>,
    exon_idx: usize,
    ignore_errors: bool,
) -> bool {
    if let Some(slice) = seq.get(range.clone()) {
        target.extend_from_slice(slice);
        return true;
    }

    if ignore_errors {
        eprintln!(
            "WARN: out-of-bounds slice for {} exon {}: {:?} (seq_len={})",
            record,
            exon_idx,
            range,
            seq.len()
        );
        false
    } else {
        panic!(
            "ERROR: out-of-bounds slice for {} exon {}: {:?} (seq_len={})",
            record,
            exon_idx,
            range,
            seq.len()
        );
    }
}

/// Extracts genomic sequence for a feature with flanking regions applied.
fn extract_seq(
    record: &GenePred,
    seq: &[u8],
    upstream_flank: usize,
    downstream_flank: usize,
    feature_type: &Feature,
    ignore_errors: bool,
) -> Option<Vec<u8>> {
    let feature = match feature_type {
        Feature::Transcript => {
            let mut f = record.exons();
            f.extend(record.introns());
            f.sort_by(|a, b| a.0.cmp(&b.0));
            f
        }
        Feature::Exon => record.exons(),
        Feature::Intron => record.introns(),
        Feature::CDS => record.coding_exons(),
        Feature::UTR => record.utr_exons(),
    };

    let feature_count = feature.len();
    let mut target = Vec::new();

    for (idx, (start, end)) in feature.iter().enumerate() {
        let feature_start = *start as usize;
        let feature_end = *end as usize;

        let range = match slice_range_for_exon(
            idx,
            feature_count,
            feature_start,
            feature_end,
            upstream_flank,
            downstream_flank,
        ) {
            Ok(r) => r,
            Err(err) => {
                if ignore_errors {
                    eprintln!("WARN: {:?} for record {}", err, record);
                    return None;
                } else {
                    panic!("ERROR: {:?} for record {}", err, record);
                }
            }
        };

        if !extend_or_handle_oob(record, &mut target, seq, range, idx, ignore_errors) {
            return None;
        }
    }

    Some(target)
}

/// Supported genomic annotation file formats.
#[derive(Clone, Copy)]
enum RegionFormat {
    Bed,
    Gtf,
    Gff,
}

/// Detects the genomic annotation format from file extension.
fn detect_region_format(path: &Path) -> Option<RegionFormat> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("bed") => Some(RegionFormat::Bed),
        Some("gtf") => Some(RegionFormat::Gtf),
        Some("gff") => Some(RegionFormat::Gff),
        Some("gz") => {
            let stem = path.file_stem()?.to_str()?;
            if stem.ends_with(".bed") {
                Some(RegionFormat::Bed)
            } else if stem.ends_with(".gtf") {
                Some(RegionFormat::Gtf)
            } else if stem.ends_with(".gff") {
                Some(RegionFormat::Gff)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Loads genome sequences from a file (2bit or FASTA format).
pub fn get_sequences(sequence: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    info!("Reading sequences from file {}", sequence.display());
    match sequence.extension() {
        Some(ext) => match ext.to_str() {
            Some("2bit") => from_2bit(sequence),
            Some("fa") | Some("gz") => from_fa(sequence),
            _ => panic!("ERROR: Unsupported file format"),
        },
        None => panic!("ERROR: No file extension"),
    }
}

/// Loads genome sequences from a 2bit compressed format file.
fn from_2bit(twobit: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut genome = TwoBitFile::open_and_read(&twobit).expect("ERROR: Cannot open 2bit file");

    let mut sequences = HashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .unwrap_or_else(|e| panic!("ERROR: {}", e))
            .as_bytes()
            .to_vec();

        sequences.insert(chr.as_bytes().to_vec(), seq);
    });

    info!(
        "Read {} sequences from file {:}",
        sequences.len(),
        twobit.display()
    );

    sequences
}

/// Loads genome sequences from a FASTA format file (optionally gzipped).
pub fn from_fa<F: AsRef<Path> + Debug>(f: F) -> HashMap<Vec<u8>, Vec<u8>> {
    let path = f.as_ref();
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("ERROR: cannot open FASTA {}: {}", path.display(), e));

    let mut reader: Box<dyn BufRead> = match path.extension().and_then(|ext| ext.to_str()) {
        Some("gz") => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    };

    let mut acc = HashMap::new();
    let mut line = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut seq = Vec::new();

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .unwrap_or_else(|e| panic!("ERROR: cannot read FASTA {}: {}", path.display(), e));

        if bytes_read == 0 {
            break;
        }

        if line.ends_with(b"\n") {
            line.pop();
        }

        if line.ends_with(b"\r") {
            line.pop();
        }

        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if let Some(prev_header) = header.replace(line[1..].to_vec()) {
                acc.insert(prev_header, std::mem::take(&mut seq));
            }
        } else {
            seq.extend_from_slice(&line);
        }
    }

    if let Some(last_header) = header {
        acc.insert(last_header, seq);
    }

    info!("Read {} sequences from file {:#?}", acc.len(), f);

    acc
}

/// Translates a sequence into amino acids.
fn translate(sequence: &[u8]) -> Vec<u8> {
    let mut aa = Vec::new();

    for codon in sequence.chunks(3) {
        if codon.len() != 3 {
            break;
        }

        if codon.iter().any(|&b| !is_unambiguous_dna_base(b)) {
            aa.push(b'X');
            continue;
        }

        let amino_acid = translate_codon(codon);
        if amino_acid == b'X' {
            panic!(
                "ERROR: codon -> {:?} is not a valid codon from sequence -> {:?}!",
                std::str::from_utf8(codon).unwrap(),
                std::str::from_utf8(sequence).unwrap()
            );
        }

        aa.push(amino_acid);
    }

    aa
}

/// Checks if a base is unambiguous DNA.
fn is_unambiguous_dna_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Translates a codon into an amino acid.
fn translate_codon(codon: &[u8]) -> u8 {
    for (table_codon, amino_acid) in &CODON_TABLE {
        if codon == *table_codon {
            return *amino_acid;
        }
    }

    b'X' // INFO: unknown codon
}
