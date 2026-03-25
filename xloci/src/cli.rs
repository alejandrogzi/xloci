//! get sequences from 2bit/fa using bed/gtf/gff
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This tool provides an easy way to get any sequence (exon, intron, cds, utr, etc.)
//! completely agnostic of the underlying format, either for reference sequences (2bit, fa, fa.gz)
//! or regions (bed, gtf, gff, gz, bz2, zstd)

use clap::{ArgAction, Parser, ValueEnum};
use log::Level;

use std::{fmt, path::PathBuf, str::FromStr};

#[derive(Parser, Debug)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    long_about = None
)]
pub struct Args {
    /// Path to genome sequence file (.fa, .fa.gz, or .2bit); reads from stdin when omitted
    #[arg(short = 's', long, default_value = "-", hide_default_value = true)]
    pub sequence: PathBuf,

    /// Path to genomic regions file (BED, GTF, or GFF format)
    #[arg(short = 'r', long)]
    pub regions: PathBuf,

    /// Output directory for extracted sequences
    #[arg(short = 'o', long)]
    pub outdir: PathBuf,

    /// Number of records per parallel processing chunk
    #[arg(short = 'c', long, default_value = "1000")]
    pub chunks: usize,

    /// Bases to extend upstream of features
    #[arg(short = 'u', long, default_value = "0")]
    pub upstream_flank: usize,

    /// Bases to extend downstream of features
    #[arg(short = 'd', long, default_value = "0")]
    pub downstream_flank: usize,

    /// Type of genomic feature to extract
    #[arg(short = 'f', long, value_enum, default_value = "exon")]
    pub feature: Feature,

    /// Continue processing on errors instead of panicking
    #[arg(short = 'I', long, default_value = "false", action = ArgAction::SetTrue)]
    pub ignore_errors: bool,

    /// Logging verbosity level
    #[arg(short = 'L', long, default_value = "info")]
    pub level: Level,

    /// Stem for output FASTA files (writes <prefix>.fa or <prefix>.fa.gz)
    #[arg(short = 'p', long, default_value = "output")]
    pub prefix: String,

    /// Translate sequences to protein
    #[arg(short = 'X', long, default_value = "false", action = ArgAction::SetTrue)]
    pub translate: bool,

    /// Emit one output record per extracted feature piece
    #[arg(short = 'S', long = "split-extraction", default_value = "false", action = ArgAction::SetTrue)]
    pub split_extraction: bool,

    /// Write tab-separated output instead of FASTA
    #[arg(long = "as-tsv", default_value = "false", action = ArgAction::SetTrue)]
    pub as_tsv: bool,

    /// Separate flank columns in TSV output (requires --as-tsv and at least one flank)
    #[arg(long = "add-tab", default_value = "false", action = ArgAction::SetTrue)]
    pub add_tab: bool,

    /// Use genomic coordinates as identifiers instead of record names
    #[arg(short = 'G', long = "generic-id", default_value = "false", action = ArgAction::SetTrue)]
    pub generic_id: bool,

    /// Keep chunk outputs and skip merging into a single file
    #[arg(short = 'A', long = "as-chunk", default_value = "false", action = ArgAction::SetTrue)]
    pub as_chunk: bool,

    /// Also emit chunked BED outputs (requires --as-chunk)
    #[arg(short = 'B', long = "include-bed", requires = "as_chunk", default_value = "false", action = ArgAction::SetTrue)]
    pub include_bed: bool,

    /// Gzip-compress output files
    #[arg(short = 'Z', long, default_value = "false", action = ArgAction::SetTrue)]
    pub compress: bool,

    /// Number of threads to use
    #[clap(
        short = 't',
        long,
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,
}

/// Formats the Args struct as a comma-separated string of key=value pairs.
///
/// # Arguments
///
/// - `f`: The formatter to write to
///
/// # Example
///
/// ```rust,ignore
/// use xloci::Args;
/// let args = Args::parse();
/// println!("{}", args);
/// ```
impl fmt::Display for Args {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let sequence = if self.sequence == *"-" {
            "<stdin>".to_string()
        } else {
            self.sequence.display().to_string()
        };

        write!(
            f,
            "sequence={}, regions={}, outdir={}, chunks={}, upstream_flank={}, downstream_flank={}, feature={:?}, ignore_errors={}, level={}, prefix={}, translate={}, split_extraction={}, as_tsv={}, add_tab={}, generic_id={}, as_chunk={}, include_bed={}, compress={}",
            sequence,
            self.regions.display(),
            self.outdir.display(),
            self.chunks,
            self.upstream_flank,
            self.downstream_flank,
            self.feature,
            self.ignore_errors,
            self.level,
            self.prefix,
            self.translate,
            self.split_extraction,
            self.as_tsv,
            self.add_tab,
            self.generic_id,
            self.as_chunk,
            self.include_bed,
            self.compress,
        )
    }
}

/// Represents the type of genomic feature to extract from annotations.
///
/// # Variants
///
/// - `Transcript`: Full transcript (exons + introns)
/// - `Exon`: Exonic regions only
/// - `Intron`: Intronic regions only
/// - `CDS`: Coding sequences
/// - `UTR`: Untranslated regions
///
/// # Example
///
/// ```rust,ignore
/// use xloci::Feature;
/// use std::str::FromStr;
///
/// let feature = Feature::from_str("exon").unwrap();
/// assert_eq!(feature, Feature::Exon);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, Default)]
pub enum Feature {
    Transcript,
    #[default]
    Exon,
    Intron,
    CDS,
    UTR,
}

/// Parses a string into a Feature variant.
///
/// # Arguments
///
/// - `s`: The string to parse ("transcript", "exon", "intron", "cds", "utr")
///
/// # Example
///
/// ```rust,ignore
/// use xloci::Feature;
/// use std::str::FromStr;
///
/// let feature = Feature::from_str("cds");
/// assert_eq!(feature, Ok(Feature::CDS));
///
/// let invalid = Feature::from_str("invalid");
/// assert!(invalid.is_err());
/// ```
impl FromStr for Feature {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "transcript" => Ok(Feature::Transcript),
            "exon" => Ok(Feature::Exon),
            "intron" => Ok(Feature::Intron),
            "cds" => Ok(Feature::CDS),
            "utr" => Ok(Feature::UTR),
            _ => Err(format!("ERROR: Invalid feature type: {}", s)),
        }
    }
}
