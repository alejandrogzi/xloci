//! get sequences from 2bit/fa using bed/gtf/gff
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This tool provides an easy way to get any sequence (exon, intron, cds, utr, etc.)
//! completely agnostic of the underlying format, either for reference sequences (2bit, fa, fa.gz)
//! or regions (bed, gtf, gff, gz, bz2, zstd)

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
    borrow::Cow,
    collections::HashMap,
    fmt::Debug,
    fs::{File, create_dir_all},
    io::{BufRead, BufReader, BufWriter, Cursor, Read, Seek, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const TWOBIT_MAGIC: [u8; 4] = [0x43, 0x27, 0x41, 0x1a];
const TWOBIT_REV_MAGIC: [u8; 4] = [0x1a, 0x41, 0x27, 0x43];

/// Output format for extracted sequences.
///
/// # Variants
///
/// - `Fasta`: FASTA format
/// - `Tsv`: Tab-separated values
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::OutputFormat;
///
/// let ext = OutputFormat::Fasta.extension();
/// assert_eq!(ext, "fa");
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OutputFormat {
    Fasta,
    Tsv,
}

impl OutputFormat {
    /// Returns the file extension for the output format.
    ///
    /// # Arguments
    ///
    /// - `self`: The output format
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use xloci::core::OutputFormat;
    ///
    /// assert_eq!(OutputFormat::Fasta.extension(), "fa");
    /// assert_eq!(OutputFormat::Tsv.extension(), "tsv");
    /// ```
    fn extension(self) -> &'static str {
        match self {
            Self::Fasta => "fa",
            Self::Tsv => "tsv",
        }
    }
}

/// Configuration for output formatting and extraction behavior.
///
/// # Fields
///
/// - `format`: Output format (FASTA or TSV)
/// - `split_extraction`: Whether to emit one record per feature piece
/// - `add_tab`: Whether to separate flank columns in TSV output
/// - `generic_id`: Whether to use genomic coordinates as identifiers
/// - `translate`: Whether to translate DNA sequences to protein
/// - `unmask`: Whether to convert soft-masked bases to uppercase in output
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{OutputConfig, OutputFormat};
///
/// let config = OutputConfig {
///     format: OutputFormat::Fasta,
///     split_extraction: false,
///     add_tab: false,
///     generic_id: false,
///     translate: false,
///     unmask: false,
/// };
/// ```
#[derive(Clone, Copy, Debug)]
struct OutputConfig {
    format: OutputFormat,
    split_extraction: bool,
    add_tab: bool,
    generic_id: bool,
    translate: bool,
    unmask: bool,
}

/// Type of genomic feature piece being extracted.
///
/// # Variants
///
/// - `Exon`: Exonic region
/// - `Intron`: Intronic region
/// - `Cds`: Coding sequence
/// - `Utr`: Untranslated region
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::PieceKind;
///
/// assert_eq!(PieceKind::Exon.suffix(), b"EXON");
/// assert_eq!(PieceKind::Intron.suffix(), b"INTRON");
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PieceKind {
    Exon,
    Intron,
    Cds,
    Utr,
}

impl PieceKind {
    /// Returns the suffix bytes for this feature kind.
    ///
    /// # Arguments
    ///
    /// - `self`: The piece kind
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use xloci::core::PieceKind;
    ///
    /// assert_eq!(PieceKind::Cds.suffix(), b"CDS");
    /// assert_eq!(PieceKind::Utr.suffix(), b"UTR");
    /// ```
    fn suffix(self) -> &'static [u8] {
        match self {
            Self::Exon => b"EXON",
            Self::Intron => b"INTRON",
            Self::Cds => b"CDS",
            Self::Utr => b"UTR",
        }
    }
}

/// Represents a single extracted genomic feature piece with sequence and flanking regions.
///
/// # Fields
///
/// - `kind`: Type of feature (exon, intron, CDS, UTR)
/// - `start`: Start coordinate (0-indexed)
/// - `end`: End coordinate (exclusive)
/// - `prefix_flank`: Upstream flanking sequence
/// - `sequence`: Main feature sequence
/// - `suffix_flank`: Downstream flanking sequence
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{FeaturePiece, PieceKind};
///
/// let piece = FeaturePiece {
///     kind: PieceKind::Exon,
///     start: 100,
///     end: 200,
///     prefix_flank: b"AT".to_vec(),
///     sequence: b"ACGT".to_vec(),
///     suffix_flank: b"TA".to_vec(),
/// };
/// assert_eq!(piece.kind, PieceKind::Exon);
/// ```
#[derive(Debug)]
struct FeaturePiece {
    kind: PieceKind,
    start: usize,
    end: usize,
    prefix_flank: Vec<u8>,
    sequence: Vec<u8>,
    suffix_flank: Vec<u8>,
}

/// Represents extracted genomic feature data with pieces and flanking regions.
///
/// # Fields
///
/// - `prefix_flank`: Upstream flanking sequence (before first piece)
/// - `pieces`: Vector of feature pieces (exons, introns, etc.)
/// - `suffix_flank`: Downstream flanking sequence (after last piece)
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{ExtractedFeature, FeaturePiece, PieceKind};
///
/// let feature = ExtractedFeature::default();
/// assert!(feature.is_empty());
/// ```
#[derive(Debug, Default)]
struct ExtractedFeature {
    prefix_flank: Vec<u8>,
    pieces: Vec<FeaturePiece>,
    suffix_flank: Vec<u8>,
}

/// A split entry representing a single feature piece for output.
///
/// # Fields
///
/// - `kind`: Type of feature (exon, intron, CDS, UTR)
/// - `ordinal`: Ordinal number within its feature type
/// - `start`: Start coordinate
/// - `end`: End coordinate
/// - `prefix_flank`: Upstream flanking sequence
/// - `sequence`: Feature sequence
/// - `suffix_flank`: Downstream flanking sequence
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::SplitEntry;
/// ```
struct SplitEntry<'a> {
    kind: PieceKind,
    ordinal: usize,
    start: usize,
    end: usize,
    prefix_flank: &'a [u8],
    sequence: &'a [u8],
    suffix_flank: &'a [u8],
}

impl ExtractedFeature {
    /// Checks if the extracted feature contains any pieces.
    ///
    /// # Arguments
    ///
    /// - `self`: Reference to the extracted feature
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use xloci::core::ExtractedFeature;
    ///
    /// let feature = ExtractedFeature::default();
    /// assert!(feature.is_empty());
    /// ```
    fn is_empty(&self) -> bool {
        self.pieces.is_empty()
    }

    /// Joins only the core feature sequences without flanking regions.
    ///
    /// # Arguments
    ///
    /// - `self`: Reference to the extracted feature
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use xloci::core::ExtractedFeature;
    ///
    /// let feature = ExtractedFeature::default();
    /// let core = feature.joined_core();
    /// ```
    fn joined_core(&self) -> Vec<u8> {
        let total_len = self.pieces.iter().map(|piece| piece.sequence.len()).sum();
        let mut sequence = Vec::with_capacity(total_len);

        self.pieces
            .iter()
            .for_each(|piece| sequence.extend_from_slice(&piece.sequence));

        sequence
    }

    /// Joins all sequences including flanking regions.
    ///
    /// # Arguments
    ///
    /// - `self`: Reference to the extracted feature
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use xloci::core::ExtractedFeature;
    ///
    /// let feature = ExtractedFeature::default();
    /// let seq = feature.joined_sequence();
    /// ```
    fn joined_sequence(&self) -> Vec<u8> {
        let total_len = self.prefix_flank.len()
            + self.suffix_flank.len()
            + self
                .pieces
                .iter()
                .map(|piece| piece.sequence.len())
                .sum::<usize>();
        let mut sequence = Vec::with_capacity(total_len);

        sequence.extend_from_slice(&self.prefix_flank);
        self.pieces
            .iter()
            .for_each(|piece| sequence.extend_from_slice(&piece.sequence));
        sequence.extend_from_slice(&self.suffix_flank);

        sequence
    }

    /// Converts pieces into split entries with ordinal indices per feature type.
    ///
    /// # Arguments
    ///
    /// - `self`: Reference to the extracted feature
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use xloci::core::ExtractedFeature;
    ///
    /// let feature = ExtractedFeature::default();
    /// let entries = feature.split_entries();
    /// ```
    fn split_entries(&self) -> Vec<SplitEntry<'_>> {
        let mut exon_idx = 0;
        let mut intron_idx = 0;
        let mut cds_idx = 0;
        let mut utr_idx = 0;

        self.pieces
            .iter()
            .map(|piece| {
                let ordinal = match piece.kind {
                    PieceKind::Exon => {
                        exon_idx += 1;
                        exon_idx
                    }
                    PieceKind::Intron => {
                        intron_idx += 1;
                        intron_idx
                    }
                    PieceKind::Cds => {
                        cds_idx += 1;
                        cds_idx
                    }
                    PieceKind::Utr => {
                        utr_idx += 1;
                        utr_idx
                    }
                };

                SplitEntry {
                    kind: piece.kind,
                    ordinal,
                    start: piece.start,
                    end: piece.end,
                    prefix_flank: &piece.prefix_flank,
                    sequence: &piece.sequence,
                    suffix_flank: &piece.suffix_flank,
                }
            })
            .collect()
    }

    /// Returns the minimum start and maximum end coordinates across all pieces.
    ///
    /// # Arguments
    ///
    /// - `self`: Reference to the extracted feature
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use xloci::core::ExtractedFeature;
    ///
    /// let feature = ExtractedFeature::default();
    /// if let Some((start, end)) = feature.bounds() {
    ///     println!("Bounds: {} - {}", start, end);
    /// }
    /// ```
    fn bounds(&self) -> Option<(usize, usize)> {
        let start = self.pieces.iter().map(|piece| piece.start).min()?;
        let end = self.pieces.iter().map(|piece| piece.end).max()?;
        Some((start, end))
    }
}

/// Main processing function that orchestrates genomic sequence extraction.
///
/// # Arguments
///
/// - `args`: Command-line arguments containing all processing configuration
///
/// # Example
///
/// ```rust,ignore
/// use xloci::{Args, xloci};
/// use clap::Parser;
///
/// let args = Args::parse_from([
///     "xloci",
///     "-s", "genome.2bit",
///     "-r", "regions.gtf",
///     "-o", "output/",
/// ]);
/// xloci(args);
/// ```
pub fn xloci(args: Args) {
    validate_args(&args);

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
        unmask,
        split_extraction,
        as_tsv,
        add_tab,
        generic_id,
        as_chunk,
        include_bed,
        compress,
        ..
    } = args;

    let output = OutputConfig {
        format: if as_tsv {
            OutputFormat::Tsv
        } else {
            OutputFormat::Fasta
        },
        split_extraction,
        add_tab,
        generic_id,
        translate,
        unmask,
    };

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
            output,
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
            output,
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
            output,
            as_chunk,
            include_bed,
            compress,
        ),
        None => panic!("ERROR: Unsupported file format"),
    }
}

/// Validates command-line arguments for mutual compatibility.
///
/// # Arguments
///
/// - `args`: Command-line arguments to validate
///
/// # Example
///
/// ```rust,ignore
/// use xloci::Args;
/// use clap::Parser;
///
/// let args = Args::parse_from([
///     "xloci", "-r", "test.bed", "-o", "out/"
/// ]);
/// // validate_args(&args); // Would panic if args are invalid
/// ```
fn validate_args(args: &Args) {
    if args.include_bed && !args.as_chunk {
        panic!("ERROR: --include-bed requires --as-chunk");
    }

    if args.add_tab && !args.as_tsv {
        panic!("ERROR: --add-tab requires --as-tsv");
    }

    if args.add_tab && args.upstream_flank == 0 && args.downstream_flank == 0 {
        panic!("ERROR: --add-tab requires --upstream-flank and/or --downstream-flank");
    }

    if args.split_extraction && args.translate {
        panic!("ERROR: --split-extraction is not compatible with --translate");
    }
}

/// Processes genomic regions in parallel chunks and writes output to FASTA file.
///
/// # Arguments
///
/// - `regions`: Path to the annotation file (BED, GTF, or GFF)
/// - `chunks`: Number of records per parallel processing chunk
/// - `outdir`: Output directory path
/// - `genome`: HashMap of chromosome names to sequences
/// - `upstream_flank`: Bases to extend upstream of first exon
/// - `downstream_flank`: Bases to extend downstream of last exon
/// - `feature_type`: Type of genomic feature to extract
/// - `ignore_errors`: Whether to continue on errors
/// - `prefix`: Stem for output FASTA file names
/// - `translate`: Whether to translate sequences to protein
/// - `as_chunk`: Keep chunks separate instead of merging
/// - `include_bed`: Also write BED outputs for each chunk
/// - `compress`: Gzip-compress output files
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::process_reader;
/// use xloci::Feature;
/// use std::collections::HashMap;
///
/// let genome: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
/// process_reader::<genepred::Gtf>(
///     std::path::Path::new("regions.gtf"),
///     1000,
///     std::path::Path::new("output/"),
///     &genome,
///     0,
///     0,
///     Feature::Exon,
///     false,
///     "output",
///     false,
///     false,
///     false,
///     false,
/// );
/// ```
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
    output: OutputConfig,
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
                output,
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

    let output_path = output_path(outdir, prefix, output.format, compress);
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

/// Reads genomic regions from an annotation file.
///
/// # Arguments
///
/// - `regions`: Path to the annotation file
///
/// # Example
///
/// ```rust,ignore
/// use genepred::Gtf;
/// let reader = read_regions::<Gtf>(std::path::Path::new("regions.gtf"));
/// ```
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

/// Checks if a file path indicates a compressed file based on extension.
///
/// # Arguments
///
/// - `path`: Path to check for compression extension
///
/// # Example
///
/// ```rust,ignore
/// use std::path::Path;
/// assert!(is_compressed_path(Path::new("file.gz")));
/// assert!(is_compressed_path(Path::new("file.zst")));
/// assert!(!is_compressed_path(Path::new("file.bed")));
/// ```
fn is_compressed_path(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("gz" | "zst" | "zstd" | "bz2" | "bzip2")
    )
}

/// Merges chunk files into a single output and removes the temporary files.
///
/// # Arguments
///
/// - `chunk_paths`: Slice of (index, path) tuples for chunk files
/// - `writer`: Writer to merge chunks into
///
/// # Example
///
/// ```rust,ignore
/// use std::io::BufWriter;
/// use std::fs::File;
///
/// let chunks = vec![(0, PathBuf::from("tmp_0.fa")), (1, PathBuf::from("tmp_1.fa"))];
/// let output = File::create("merged.fa").unwrap();
/// let mut writer = BufWriter::new(output);
/// merge_and_cleanup_chunks(&chunks, &mut writer);
/// ```
fn merge_and_cleanup_chunks<W: Write>(chunk_paths: &[(usize, PathBuf)], writer: &mut W) {
    for (_, path) in chunk_paths {
        let mut chunk_file = File::open(path)
            .unwrap_or_else(|e| panic!("ERROR: cannot open chunk {}: {}", path.display(), e));
        std::io::copy(&mut chunk_file, writer).unwrap_or_else(|e| panic!("{}", e));
        std::fs::remove_file(path)
            .unwrap_or_else(|e| panic!("ERROR: cannot remove chunk {}: {}", path.display(), e));
    }
}

/// Adds .gz extension to a path if compression is enabled and not already present.
///
/// # Arguments
///
/// - `path`: Original file path
/// - `compress`: Whether to add gzip extension
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let path = PathBuf::from("output.fa");
/// assert_eq!(with_gzip_extension(path.clone(), true), PathBuf::from("output.fa.gz"));
/// assert_eq!(with_gzip_extension(path, false), PathBuf::from("output.fa"));
/// ```
fn with_gzip_extension(mut path: PathBuf, compress: bool) -> PathBuf {
    if compress && path.extension().and_then(|ext| ext.to_str()) != Some("gz") {
        path.as_mut_os_string().push(".gz");
    }

    path
}

/// Constructs the output file path from directory, prefix, format, and compression.
///
/// # Arguments
///
/// - `outdir`: Output directory
/// - `prefix`: Filename prefix
/// - `format`: Output format
/// - `compress`: Whether to gzip compress
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{output_path, OutputFormat};
/// use std::path::Path;
///
/// let path = output_path(Path::new("out/"), "sequences", OutputFormat::Fasta, false);
/// assert!(path.to_str().unwrap().ends_with(".fa"));
/// ```
fn output_path(outdir: &Path, prefix: &str, format: OutputFormat, compress: bool) -> PathBuf {
    let stem = prefix
        .strip_suffix(".fa.gz")
        .or_else(|| prefix.strip_suffix(".fa"))
        .or_else(|| prefix.strip_suffix(".tsv.gz"))
        .or_else(|| prefix.strip_suffix(".tsv"))
        .unwrap_or(prefix);
    let path = outdir.join(format!("{stem}.{}", format.extension()));
    with_gzip_extension(path, compress)
}

/// Processes a chunk of genomic records and writes extracted sequences to a temporary file.
///
/// # Arguments
///
/// - `idx`: Chunk index for naming output files
/// - `chunk`: Vector of GenePred records to process
/// - `genome`: HashMap of chromosome names to sequences
/// - `upstream_flank`: Bases to extend upstream of first exon
/// - `downstream_flank`: Bases to extend downstream of last exon
/// - `feature_type`: Type of genomic feature to extract
/// - `ignore_errors`: Whether to continue on errors
/// - `to_protein`: Whether to translate sequences to protein
/// - `collector`: Optional shared collector for chunk paths
/// - `outdir`: Output directory for chunk files
/// - `include_bed`: Whether to write BED output
/// - `compress`: Whether to gzip-compress output
///
/// # Example
///
/// ```rust,ignore
/// use std::sync::{Arc, Mutex};
/// use xloci::Feature;
///
/// let collector = Some(Arc::new(Mutex::new(Vec::new())));
/// write_chunk(
///     0,
///     vec![],
///     &genome,
///     0,
///     0,
///     &Feature::Exon,
///     false,
///     false,
///     collector,
///     Path::new("tmp/"),
///     false,
///     false,
/// );
/// ```
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
    output: OutputConfig,
    collector: Option<Arc<Mutex<Vec<(usize, PathBuf)>>>>,
    outdir: &Path,
    include_bed: bool,
    compress: bool,
) {
    info!("Processing chunk {}", idx);

    let chunk_base = outdir.join(format!("tmp_{}", idx));
    let output_path = if compress {
        chunk_base.with_extension(format!("{}.gz", output.format.extension()))
    } else {
        chunk_base.with_extension(output.format.extension())
    };

    let output_file = File::create(&output_path)
        .unwrap_or_else(|e| panic!("ERROR: cannot create {}: {}", output_path.display(), e));
    let mut writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(
            BufWriter::new(output_file),
            Compression::default(),
        ))
    } else {
        Box::new(BufWriter::new(output_file))
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
                let keys = genome
                    .keys()
                    .map(|k| std::str::from_utf8(k).unwrap())
                    .collect::<Vec<_>>();

                log::error!(
                    "ERROR: Chromosome {} from record {} not found in genome with keys {:?}!",
                    String::from_utf8_lossy(&record.chrom),
                    record,
                    keys
                );
                std::process::exit(1);
            });

            let extracted = extract_feature(
                &record,
                seq,
                upstream_flank,
                downstream_flank,
                feature_type,
                output.split_extraction,
                ignore_errors,
            );

            if let Some(extracted) = extracted {
                write_record_output(&mut writer, &record, &extracted, output);

                if let Some(bed_writer) = &mut bed_writer {
                    Writer::<Bed12>::from_record(&record, bed_writer)
                        .unwrap_or_else(|e| panic!("{}", e));
                }
            }
        });

    writer
        .flush()
        .unwrap_or_else(|e| panic!("ERROR: cannot flush {}: {}", output_path.display(), e));

    if let Some(bed_writer) = &mut bed_writer {
        bed_writer
            .flush()
            .unwrap_or_else(|e| panic!("ERROR: cannot flush BED writer: {}", e));
    }

    if let Some(collector) = collector {
        collector
            .lock()
            .unwrap_or_else(|e| panic!("ERROR: Cannot acquire lock on collector: {}", e))
            .push((idx, output_path));
    }
}

/// Error type for range calculation failures during sequence extraction.
///
/// # Variants
///
/// - `Underflow`: Coordinate minus flank would result in negative value
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::RangeError;
///
/// let err = RangeError::Underflow {
///     feature_coord: 5,
///     flank: 10,
/// };
/// println!("{}", err); // "Feature coordinate 5 is underflowing by 10 bases"
/// ```
#[derive(Debug)]
#[allow(dead_code)]
enum RangeError {
    Underflow { feature_coord: usize, flank: usize },
    Overflow { feature_coord: usize, flank: usize },
}

/// Formats the RangeError as a human-readable error message.
///
/// # Arguments
///
/// - `f`: The formatter to write to
///
/// # Example
///
/// ```rust,ignore
/// let err = RangeError::Underflow { feature_coord: 5, flank: 10 };
/// assert!(err.to_string().contains("underflowing"));
/// ```
impl std::fmt::Display for RangeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RangeError::Underflow {
                feature_coord,
                flank,
            } => write!(
                f,
                "Feature coordinate {} is underflowing by {} bases",
                feature_coord, flank
            ),
            RangeError::Overflow {
                feature_coord,
                flank,
            } => write!(
                f,
                "Feature coordinate {} is overflowing by {} bases",
                feature_coord, flank
            ),
        }
    }
}

/// Extracts genomic intervals for a given feature type from a record.
///
/// # Arguments
///
/// - `record`: GenePred record containing annotation
/// - `feature_type`: Type of feature to extract
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::feature_intervals;
/// use xloci::Feature;
/// use genepred::GenePred;
/// // let intervals = feature_intervals(&record, &Feature::Exon);
/// ```
fn feature_intervals(record: &GenePred, feature_type: &Feature) -> Vec<(PieceKind, usize, usize)> {
    match feature_type {
        Feature::Transcript => {
            let mut intervals = Vec::new();

            record.exons().into_iter().for_each(|(start, end)| {
                intervals.push((PieceKind::Exon, start as usize, end as usize));
            });
            record.introns().into_iter().for_each(|(start, end)| {
                intervals.push((PieceKind::Intron, start as usize, end as usize));
            });

            intervals.sort_by_key(|(_, start, end)| (*start, *end));
            intervals
        }
        Feature::Exon => record
            .exons()
            .into_iter()
            .map(|(start, end)| (PieceKind::Exon, start as usize, end as usize))
            .collect(),
        Feature::Intron => record
            .introns()
            .into_iter()
            .map(|(start, end)| (PieceKind::Intron, start as usize, end as usize))
            .collect(),
        Feature::CDS => record
            .coding_exons()
            .into_iter()
            .map(|(start, end)| (PieceKind::Cds, start as usize, end as usize))
            .collect(),
        Feature::UTR => record
            .utr_exons()
            .into_iter()
            .map(|(start, end)| (PieceKind::Utr, start as usize, end as usize))
            .collect(),
    }
}

/// Extracts a sequence range from a chromosome sequence.
///
/// # Arguments
///
/// - `record`: GenePred record for error messaging
/// - `seq`: Chromosome sequence
/// - `range`: Range to extract
/// - `label`: Label for error messages
/// - `ignore_errors`: Whether to return None on error or panic
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::extract_range;
/// // let seq = extract_range(&record, b"ACGTACGT", 0..4, "test", false);
/// // assert_eq!(seq, Some(b"ACGT".to_vec()));
/// ```
fn extract_range(
    record: &GenePred,
    seq: &[u8],
    range: std::ops::Range<usize>,
    label: &str,
    ignore_errors: bool,
) -> Option<Vec<u8>> {
    if let Some(slice) = seq.get(range.clone()) {
        return Some(slice.to_vec());
    }

    if ignore_errors {
        eprintln!(
            "WARN: out-of-bounds slice for {} {}: {:?} (seq_len={})",
            record,
            label,
            range,
            seq.len()
        );
        None
    } else {
        panic!(
            "ERROR: out-of-bounds slice for {} {}: {:?} (seq_len={})",
            record,
            label,
            range,
            seq.len()
        );
    }
}

/// Extracts genomic features from a chromosome sequence.
///
/// # Arguments
///
/// - `record`: GenePred record containing annotation
/// - `seq`: Chromosome sequence
/// - `upstream_flank`: Bases to extend upstream
/// - `downstream_flank`: Bases to extend downstream
/// - `feature_type`: Type of feature to extract
/// - `split_extraction`: Whether to split extraction per piece
/// - `ignore_errors`: Whether to continue on errors
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::extract_feature;
/// use xloci::Feature;
/// // let result = extract_feature(&record, seq, 100, 50, &Feature::Exon, false, false);
/// ```
fn extract_feature(
    record: &GenePred,
    seq: &[u8],
    upstream_flank: usize,
    downstream_flank: usize,
    feature_type: &Feature,
    split_extraction: bool,
    ignore_errors: bool,
) -> Option<ExtractedFeature> {
    let intervals = feature_intervals(record, feature_type);

    if intervals.is_empty() {
        return Some(ExtractedFeature::default());
    }

    let (_, first_start, _) = intervals
        .first()
        .copied()
        .unwrap_or_else(|| panic!("ERROR: missing first interval"));
    let (_, _, last_end) = intervals
        .last()
        .copied()
        .unwrap_or_else(|| panic!("ERROR: missing last interval"));

    let prefix_flank = if !split_extraction && upstream_flank > 0 {
        let flank_start = match first_start.checked_sub(upstream_flank) {
            Some(start) => start,
            None => {
                let err = RangeError::Underflow {
                    feature_coord: first_start,
                    flank: upstream_flank,
                };
                if ignore_errors {
                    eprintln!("WARN: {} for record {}", err, record);
                    return None;
                } else {
                    panic!("ERROR: {} for record {}", err, record);
                }
            }
        };

        extract_range(
            record,
            seq,
            flank_start..first_start,
            "prefix flank",
            ignore_errors,
        )?
    } else {
        Vec::new()
    };

    let mut pieces = Vec::with_capacity(intervals.len());
    for (idx, (kind, start, end)) in intervals.into_iter().enumerate() {
        let label = match kind {
            PieceKind::Exon => format!("exon {}", idx + 1),
            PieceKind::Intron => format!("intron {}", idx + 1),
            PieceKind::Cds => format!("cds {}", idx + 1),
            PieceKind::Utr => format!("utr {}", idx + 1),
        };
        let prefix_flank = if split_extraction && upstream_flank > 0 {
            let Some(flank_start) = start.checked_sub(upstream_flank) else {
                let err = RangeError::Underflow {
                    feature_coord: start,
                    flank: upstream_flank,
                };
                if ignore_errors {
                    eprintln!(
                        "WARN: skipping {} for record {}: {}; flank could not be reached",
                        label, record, err
                    );
                    continue;
                } else {
                    panic!("ERROR: {} for record {}", err, record);
                }
            };

            match extract_range(
                record,
                seq,
                flank_start..start,
                &format!("{} prefix flank", label),
                ignore_errors,
            ) {
                Some(flank) => flank,
                None => continue,
            }
        } else {
            Vec::new()
        };
        let Some(sequence) = extract_range(record, seq, start..end, &label, ignore_errors) else {
            continue;
        };
        let suffix_flank = if split_extraction && downstream_flank > 0 {
            let Some(flank_end) = end.checked_add(downstream_flank) else {
                let err = RangeError::Overflow {
                    feature_coord: end,
                    flank: downstream_flank,
                };
                if ignore_errors {
                    eprintln!(
                        "WARN: skipping {} for record {}: {}; flank could not be reached",
                        label, record, err
                    );
                    continue;
                } else {
                    panic!("ERROR: {} for record {}", err, record);
                }
            };

            match extract_range(
                record,
                seq,
                end..flank_end,
                &format!("{} suffix flank", label),
                ignore_errors,
            ) {
                Some(flank) => flank,
                None => continue,
            }
        } else {
            Vec::new()
        };
        pieces.push(FeaturePiece {
            kind,
            start,
            end,
            prefix_flank,
            sequence,
            suffix_flank,
        });
    }

    let suffix_flank = if !split_extraction && downstream_flank > 0 {
        let flank_end = match last_end.checked_add(downstream_flank) {
            Some(end) => end,
            None => {
                let err = RangeError::Overflow {
                    feature_coord: last_end,
                    flank: downstream_flank,
                };
                if ignore_errors {
                    eprintln!("WARN: {} for record {}", err, record);
                    return None;
                } else {
                    panic!("ERROR: {} for record {}", err, record);
                }
            }
        };

        extract_range(
            record,
            seq,
            last_end..flank_end,
            "suffix flank",
            ignore_errors,
        )?
    } else {
        Vec::new()
    };

    Some(orient_extracted_feature(
        record,
        ExtractedFeature {
            prefix_flank,
            pieces,
            suffix_flank,
        },
    ))
}

/// Orients extracted features based on strand direction.
///
/// Reverse strand features are reverse complemented and reordered.
///
/// # Arguments
///
/// - `record`: GenePred record containing strand information
/// - `extracted`: Extracted feature data
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::orient_extracted_feature;
/// // let oriented = orient_extracted_feature(&record, extracted);
/// ```
fn orient_extracted_feature(record: &GenePred, extracted: ExtractedFeature) -> ExtractedFeature {
    match record.strand() {
        Some(Strand::Reverse) => ExtractedFeature {
            prefix_flank: reverse_complement(extracted.suffix_flank),
            pieces: extracted
                .pieces
                .into_iter()
                .rev()
                .map(|mut piece| {
                    let prefix_flank = reverse_complement(piece.suffix_flank);
                    let suffix_flank = reverse_complement(piece.prefix_flank);
                    reverse_complement_in_place(&mut piece.sequence);
                    piece.prefix_flank = prefix_flank;
                    piece.suffix_flank = suffix_flank;
                    piece
                })
                .collect(),
            suffix_flank: reverse_complement(extracted.prefix_flank),
        },
        Some(Strand::Forward) | Some(Strand::Unknown) | None => extracted,
    }
}

/// Returns the reverse complement of a DNA sequence.
///
/// # Arguments
///
/// - `sequence`: DNA sequence to reverse complement
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::reverse_complement;
///
/// assert_eq!(reverse_complement(b"ATGC".to_vec()), b"GCAT");
/// ```
fn reverse_complement(mut sequence: Vec<u8>) -> Vec<u8> {
    reverse_complement_in_place(&mut sequence);
    sequence
}

/// Reverses a sequence and replaces each base with its complement in place.
///
/// # Arguments
///
/// - `sequence`: DNA sequence to reverse complement
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::reverse_complement_in_place;
///
/// let mut seq = b"ATGC".to_vec();
/// reverse_complement_in_place(&mut seq);
/// assert_eq!(seq, b"GCAT");
/// ```
fn reverse_complement_in_place(sequence: &mut [u8]) {
    sequence.reverse();
    sequence.iter_mut().for_each(|base| {
        *base = complement_base(*base);
    });
}

/// Returns the complement of a DNA base.
///
/// # Arguments
///
/// - `base`: DNA base byte
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::complement_base;
///
/// assert_eq!(complement_base(b'A'), b'T');
/// assert_eq!(complement_base(b'C'), b'G');
/// ```
fn complement_base(base: u8) -> u8 {
    match base {
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

/// Creates a split feature identifier from name, kind, and ordinal.
///
/// # Arguments
///
/// - `name`: Base name
/// - `kind`: Feature kind
/// - `ordinal`: Ordinal number
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{split_identifier, PieceKind};
///
/// let id = split_identifier(b"gene1", PieceKind::Exon, 1);
/// assert_eq!(id, b"gene1_EXON1");
/// ```
fn split_identifier(name: &[u8], kind: PieceKind, ordinal: usize) -> Vec<u8> {
    let ordinal = ordinal.to_string();
    let suffix = kind.suffix();
    let mut identifier = Vec::with_capacity(name.len() + suffix.len() + ordinal.len() + 1);

    identifier.extend_from_slice(name);
    identifier.push(b'_');
    identifier.extend_from_slice(suffix);
    identifier.extend_from_slice(ordinal.as_bytes());

    identifier
}

/// Returns the strand label character for a given strand.
///
/// # Arguments
///
/// - `strand`: Optional strand information
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::strand_label;
/// use genepred::Strand;
///
/// assert_eq!(strand_label(Some(Strand::Forward)), b'+');
/// assert_eq!(strand_label(Some(Strand::Reverse)), b'-');
/// assert_eq!(strand_label(None), b'.');
/// ```
fn strand_label(strand: Option<Strand>) -> u8 {
    match strand {
        Some(Strand::Forward) => b'+',
        Some(Strand::Reverse) => b'-',
        Some(Strand::Unknown) | None => b'.',
    }
}

/// Creates a generic identifier from genomic coordinates.
///
/// # Arguments
///
/// - `chrom`: Chromosome name
/// - `start`: Start coordinate
/// - `end`: End coordinate
/// - `strand`: Strand information
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::generic_identifier;
///
/// let id = generic_identifier(b"chr1", 1000, 2000, None);
/// assert_eq!(id, b"chr1:1000-2000(.)");
/// ```
fn generic_identifier(chrom: &[u8], start: usize, end: usize, strand: Option<Strand>) -> Vec<u8> {
    let start = start.to_string();
    let end = end.to_string();
    let mut identifier = Vec::with_capacity(chrom.len() + start.len() + end.len() + 5);

    identifier.extend_from_slice(chrom);
    identifier.push(b':');
    identifier.extend_from_slice(start.as_bytes());
    identifier.push(b'-');
    identifier.extend_from_slice(end.as_bytes());
    identifier.push(b'(');
    identifier.push(strand_label(strand));
    identifier.push(b')');

    identifier
}

/// Writes a FASTA record to the output.
///
/// # Arguments
///
/// - `writer`: Output writer
/// - `identifier`: Sequence identifier
/// - `sequence`: DNA sequence
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::write_fasta_record;
/// // write_fasta_record(&mut writer, b"seq1", b"ACGT");
/// ```
fn write_fasta_record(writer: &mut dyn Write, identifier: &[u8], sequence: &[u8]) {
    writer.write_all(b">").unwrap_or_else(|e| panic!("{}", e));
    writer
        .write_all(identifier)
        .unwrap_or_else(|e| panic!("{}", e));
    writer.write_all(b"\n").unwrap_or_else(|e| panic!("{}", e));
    writer
        .write_all(sequence)
        .unwrap_or_else(|e| panic!("{}", e));
    writer.write_all(b"\n").unwrap_or_else(|e| panic!("{}", e));
}

/// Writes a TSV record to the output.
///
/// # Arguments
///
/// - `writer`: Output writer
/// - `columns`: Column values
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::write_tsv_record;
/// // write_tsv_record(&mut writer, &[b"id", b"seq"]);
/// ```
fn write_tsv_record(writer: &mut dyn Write, columns: &[&[u8]]) {
    columns.iter().enumerate().for_each(|(idx, column)| {
        if idx > 0 {
            writer.write_all(b"\t").unwrap_or_else(|e| panic!("{}", e));
        }

        writer.write_all(column).unwrap_or_else(|e| panic!("{}", e));
    });
    writer.write_all(b"\n").unwrap_or_else(|e| panic!("{}", e));
}

/// Returns sequence bytes ready for output.
///
/// # Arguments
///
/// - `sequence`: Sequence bytes to write
/// - `unmask`: Whether to uppercase soft-masked bases
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::output_sequence;
/// // assert_eq!(output_sequence(b"aCgT", true).as_ref(), b"ACGT");
/// ```
fn output_sequence(sequence: &[u8], unmask: bool) -> Cow<'_, [u8]> {
    if unmask {
        Cow::Owned(sequence.to_ascii_uppercase())
    } else {
        Cow::Borrowed(sequence)
    }
}

/// Converts sequence bytes to uppercase when unmasking is enabled.
///
/// # Arguments
///
/// - `sequence`: Sequence bytes to update
/// - `unmask`: Whether to uppercase soft-masked bases
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::unmask_sequence;
/// // let mut seq = b"aCgT".to_vec();
/// // unmask_sequence(&mut seq, true);
/// ```
fn unmask_sequence(sequence: &mut [u8], unmask: bool) {
    if unmask {
        sequence.make_ascii_uppercase();
    }
}

/// Writes joined (non-split) feature output in FASTA or TSV format.
///
/// # Arguments
///
/// - `writer`: Output writer
/// - `record`: GenePred record
/// - `identifier`: Sequence identifier
/// - `extracted`: Extracted feature data
/// - `output`: Output configuration
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{write_joined_output, OutputConfig, OutputFormat};
/// // write_joined_output(&mut writer, &record, b"id", &extracted, config);
/// ```
fn write_joined_output(
    writer: &mut dyn Write,
    record: &GenePred,
    identifier: &[u8],
    extracted: &ExtractedFeature,
    output: OutputConfig,
) {
    if output.add_tab {
        let mut core_sequence = extracted.joined_core();
        unmask_sequence(&mut core_sequence, output.unmask);
        if output.translate {
            core_sequence = translate(&core_sequence);
        }

        if core_sequence.is_empty() {
            warn!("WARN: empty sequence for record {}", record);
            return;
        }

        let prefix_flank = output_sequence(&extracted.prefix_flank, output.unmask);
        let suffix_flank = output_sequence(&extracted.suffix_flank, output.unmask);

        if !extracted.prefix_flank.is_empty() && !extracted.suffix_flank.is_empty() {
            write_tsv_record(
                writer,
                &[
                    identifier,
                    prefix_flank.as_ref(),
                    &core_sequence,
                    suffix_flank.as_ref(),
                ],
            );
        } else {
            let flank = if !extracted.prefix_flank.is_empty() {
                prefix_flank.as_ref()
            } else {
                suffix_flank.as_ref()
            };

            write_tsv_record(writer, &[identifier, flank, &core_sequence]);
        }

        return;
    }

    let mut sequence = extracted.joined_sequence();
    unmask_sequence(&mut sequence, output.unmask);
    if output.translate {
        sequence = translate(&sequence);
    }

    if sequence.is_empty() {
        warn!("WARN: empty sequence for record {}", record);
        return;
    }

    match output.format {
        OutputFormat::Fasta => write_fasta_record(writer, identifier, &sequence),
        OutputFormat::Tsv => write_tsv_record(writer, &[identifier, &sequence]),
    }
}

/// Writes split feature output (one record per piece) in FASTA or TSV format.
///
/// # Arguments
///
/// - `writer`: Output writer
/// - `record`: GenePred record
/// - `base_identifier`: Base identifier for splitting
/// - `extracted`: Extracted feature data
/// - `output`: Output configuration
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{write_split_output, OutputConfig, OutputFormat};
/// // write_split_output(&mut writer, &record, Some(b"gene1"), &extracted, config);
/// ```
fn write_split_output(
    writer: &mut dyn Write,
    record: &GenePred,
    base_identifier: Option<&[u8]>,
    extracted: &ExtractedFeature,
    output: OutputConfig,
) {
    let entries = extracted.split_entries();

    if entries.is_empty() {
        warn!("WARN: empty sequence for record {}", record);
        return;
    }

    let include_prefix_flank = extracted
        .pieces
        .iter()
        .any(|piece| !piece.prefix_flank.is_empty());
    let include_suffix_flank = extracted
        .pieces
        .iter()
        .any(|piece| !piece.suffix_flank.is_empty());

    for entry in entries {
        let split_id = if output.generic_id {
            generic_identifier(record.chrom(), entry.start, entry.end, record.strand())
        } else {
            split_identifier(
                base_identifier.unwrap_or_else(|| {
                    panic!("ERROR: missing split identifier for record {}", record)
                }),
                entry.kind,
                entry.ordinal,
            )
        };

        if output.add_tab {
            let prefix_flank = output_sequence(entry.prefix_flank, output.unmask);
            let sequence = output_sequence(entry.sequence, output.unmask);
            let suffix_flank = output_sequence(entry.suffix_flank, output.unmask);

            if include_prefix_flank && include_suffix_flank {
                write_tsv_record(
                    writer,
                    &[
                        &split_id,
                        prefix_flank.as_ref(),
                        sequence.as_ref(),
                        suffix_flank.as_ref(),
                    ],
                );
            } else {
                let flank = if include_prefix_flank {
                    prefix_flank.as_ref()
                } else {
                    suffix_flank.as_ref()
                };

                write_tsv_record(writer, &[&split_id, flank, sequence.as_ref()]);
            }

            continue;
        }

        let mut sequence = Vec::with_capacity(
            entry.prefix_flank.len() + entry.sequence.len() + entry.suffix_flank.len(),
        );
        sequence.extend_from_slice(entry.prefix_flank);
        sequence.extend_from_slice(entry.sequence);
        sequence.extend_from_slice(entry.suffix_flank);
        unmask_sequence(&mut sequence, output.unmask);

        if sequence.is_empty() {
            warn!("WARN: empty sequence for record {}", record);
            continue;
        }

        match output.format {
            OutputFormat::Fasta => write_fasta_record(writer, &split_id, &sequence),
            OutputFormat::Tsv => write_tsv_record(writer, &[&split_id, &sequence]),
        }
    }
}

/// Writes extracted feature output, delegating to split or joined output.
///
/// # Arguments
///
/// - `writer`: Output writer
/// - `record`: GenePred record
/// - `extracted`: Extracted feature data
/// - `output`: Output configuration
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::{write_record_output, OutputConfig, OutputFormat};
/// // write_record_output(&mut writer, &record, &extracted, config);
/// ```
fn write_record_output(
    writer: &mut dyn Write,
    record: &GenePred,
    extracted: &ExtractedFeature,
    output: OutputConfig,
) {
    if extracted.is_empty() {
        warn!("WARN: empty sequence for record {}", record);
        return;
    }

    if output.split_extraction {
        let base_identifier = if output.generic_id {
            None
        } else {
            Some(
                record
                    .name()
                    .unwrap_or_else(|| panic!("ERROR: missing name for record {}", record)),
            )
        };
        write_split_output(writer, record, base_identifier, extracted, output);
    } else {
        let identifier = if output.generic_id {
            let (start, end) = extracted
                .bounds()
                .unwrap_or_else(|| panic!("ERROR: missing bounds for record {}", record));
            generic_identifier(record.chrom(), start, end, record.strand())
        } else {
            record
                .name()
                .unwrap_or_else(|| panic!("ERROR: missing name for record {}", record))
                .to_vec()
        };
        write_joined_output(writer, record, &identifier, extracted, output);
    }
}

/// Supported genomic annotation file formats.
///
/// # Variants
///
/// - `Bed`: BED format (12-column)
/// - `Gtf`: GTF format
/// - `Gff`: GFF format
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::RegionFormat;
///
/// let format = detect_region_format(Path::new("annotations.gtf"));
/// assert_eq!(format, Some(RegionFormat::Gtf));
/// ```
#[derive(Clone, Copy)]
enum RegionFormat {
    Bed,
    Gtf,
    Gff,
}

/// Detects the genomic annotation format from file extension.
///
/// # Arguments
///
/// - `path`: Path to the annotation file
///
/// # Example
///
/// ```rust,ignore
/// use std::path::Path;
///
/// assert_eq!(detect_region_format(Path::new("file.bed")), Some(RegionFormat::Bed));
/// assert_eq!(detect_region_format(Path::new("file.gtf")), Some(RegionFormat::Gtf));
/// assert_eq!(detect_region_format(Path::new("file.gff")), Some(RegionFormat::Gff));
/// assert_eq!(detect_region_format(Path::new("file.gtf.gz")), Some(RegionFormat::Gtf));
/// assert_eq!(detect_region_format(Path::new("file.txt")), None);
/// ```
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
///
/// # Arguments
///
/// - `sequence`: Path to the genome file (.fa, .fa.gz, or .2bit)
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let genome = get_sequences(PathBuf::from("genome.2bit"));
/// let genome = get_sequences(PathBuf::from("genome.fa"));
/// let genome = get_sequences(PathBuf::from("genome.fa.gz"));
/// ```
pub fn get_sequences(sequence: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    if sequence == *"-" {
        return from_stdin();
    }

    info!("Reading sequences from file {}", sequence.display());
    match sequence.extension() {
        Some(ext) => match ext.to_str() {
            Some("2bit") => from_2bit(sequence),
            Some("fa") | Some("fasta") | Some("fna") | Some("gz") => from_fa(sequence),
            _ => panic!("ERROR: Unsupported file format"),
        },
        None => panic!("ERROR: No file extension"),
    }
}

/// Loads genome sequences from stdin.
///
/// # Example
///
/// ```rust,ignore
/// use std::io::Write;
/// use std::process::{Command, Stdio};
///
/// let mut child = Command::new("cat")
///     .arg("genome.fa")
///     .stdout(Stdio::piped())
///     .spawn()
///     .unwrap_or_else(|e| panic!("ERROR: cannot spawn cat: {}", e));
///
/// let genome = from_stdin();
/// ```
fn from_stdin() -> HashMap<Vec<u8>, Vec<u8>> {
    info!("Reading sequences from stdin");

    let mut input = Vec::new();
    std::io::stdin()
        .read_to_end(&mut input)
        .unwrap_or_else(|e| panic!("ERROR: cannot read stdin: {}", e));

    if input.is_empty() {
        panic!("ERROR: Missing --sequence and stdin is empty");
    }

    if input.starts_with(&GZIP_MAGIC) {
        return parse_fasta_reader(
            BufReader::new(MultiGzDecoder::new(Cursor::new(input))),
            "stdin",
        );
    }

    if input.starts_with(&TWOBIT_MAGIC) || input.starts_with(&TWOBIT_REV_MAGIC) {
        return from_2bit_buf(input, "stdin");
    }

    if input
        .iter()
        .copied()
        .find(|b| !b.is_ascii_whitespace())
        .is_some_and(|b| b == b'>')
    {
        return parse_fasta_reader(BufReader::new(Cursor::new(input)), "stdin");
    }

    panic!("ERROR: Unsupported stdin sequence format");
}

/// Loads genome sequences from a 2bit compressed format file.
///
/// # Arguments
///
/// - `twobit`: Path to the 2bit file
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let sequences = from_2bit(PathBuf::from("genome.2bit"));
/// let chr1 = sequences.get(b"chr1");
/// ```
fn from_2bit(twobit: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    let genome = TwoBitFile::open_and_read(&twobit).expect("ERROR: Cannot open 2bit file");
    let source = format!("file {}", twobit.display());
    collect_2bit_sequences(genome, &source)
}

/// Loads genome sequences from a 2bit buffer.
///
/// # Arguments
///
/// - `buf`: 2bit file contents as bytes
/// - `source`: Source description for error messages
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::from_2bit_buf;
/// // let genome = from_2bit_buf(data, "stdin");
/// ```
fn from_2bit_buf(buf: Vec<u8>, source: &str) -> HashMap<Vec<u8>, Vec<u8>> {
    let genome = TwoBitFile::from_buf(buf)
        .unwrap_or_else(|e| panic!("ERROR: Cannot read 2bit from {}: {}", source, e));
    collect_2bit_sequences(genome, source)
}

/// Collects sequences from a 2bit reader into a HashMap.
///
/// # Arguments
///
/// - `genome`: TwoBitFile reader
/// - `source`: Source description for logging
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::collect_2bit_sequences;
/// ```
fn collect_2bit_sequences<R: Read + Seek>(
    mut genome: TwoBitFile<R>,
    source: &str,
) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut sequences = HashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .unwrap_or_else(|e| panic!("ERROR: {}", e))
            .as_bytes()
            .to_vec();

        sequences.insert(chr.as_bytes().to_vec(), seq);
    });

    info!("Read {} sequences from {}", sequences.len(), source);

    sequences
}

/// Loads genome sequences from a FASTA format file (optionally gzipped).
///
/// # Arguments
///
/// - `f`: Path to the FASTA file (.fa or .fa.gz)
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let sequences = from_fa(PathBuf::from("genome.fa"));
/// let sequences = from_fa(PathBuf::from("genome.fa.gz"));
/// let chr1 = sequences.get(b"chr1");
/// ```
pub fn from_fa<P: AsRef<Path>>(f: P) -> HashMap<Vec<u8>, Vec<u8>> {
    let path = f.as_ref();
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("ERROR: cannot open FASTA {}: {}", path.display(), e));

    let reader: Box<dyn BufRead> = match path.extension().and_then(|ext| ext.to_str()) {
        Some("gz") => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    };

    let source = format!("file {}", path.display());
    parse_fasta_reader(reader, &source)
}

/// Parses FASTA format reader into a HashMap of sequences.
///
/// # Arguments
///
/// - `reader`: Buffered reader containing FASTA data
/// - `source`: Source description for logging
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::parse_fasta_reader;
/// use std::io::BufReader;
/// // let genome = parse_fasta_reader(BufReader::new(file), "file.fasta");
/// ```
fn parse_fasta_reader<R: BufRead>(mut reader: R, source: &str) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut acc = HashMap::new();
    let mut line = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut seq = Vec::new();

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .unwrap_or_else(|e| panic!("ERROR: cannot read FASTA {}: {}", source, e));

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

    info!("Read {} sequences from {}", acc.len(), source);

    acc
}

/// Translates a DNA sequence into amino acids.
///
/// # Arguments
///
/// - `sequence`: DNA sequence as bytes (A, C, G, T)
///
/// # Example
///
/// ```rust,ignore
/// let dna = b"ATGGCT";
/// let protein = translate(dna);
/// assert_eq!(protein, b"MA");
/// ```
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

/// Checks if a base is an unambiguous DNA nucleotide (A, C, G, or T).
///
/// # Arguments
///
/// - `b`: Byte to check
///
/// # Example
///
/// ```rust,ignore
/// assert!(is_unambiguous_dna_base(b'A'));
/// assert!(is_unambiguous_dna_base(b'T'));
/// assert!(!is_unambiguous_dna_base(b'N'));
/// assert!(!is_unambiguous_dna_base(b'a'));
/// ```
fn is_unambiguous_dna_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Translates a single codon into an amino acid.
///
/// # Arguments
///
/// - `codon`: Three-byte slice representing a codon
///
/// # Example
///
/// ```rust,ignore
/// assert_eq!(translate_codon(b"ATG"), b'M'); // Methionine (start)
/// assert_eq!(translate_codon(b"TAA"), b'*'); // Stop codon
/// assert_eq!(translate_codon(b"TTT"), b'F'); // Phenylalanine
/// assert_eq!(translate_codon(b"???"), b'X'); // Unknown
/// ```
fn translate_codon(codon: &[u8]) -> u8 {
    for (table_codon, amino_acid) in &CODON_TABLE {
        if codon == *table_codon {
            return *amino_acid;
        }
    }

    b'X' // INFO: unknown codon
}
