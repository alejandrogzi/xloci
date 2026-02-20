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
    /// Path to genome sequence file (.fa, .fa.gz, or .2bit)
    #[arg(short = 's', long)]
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

    /// Prefix for output files
    #[arg(short = 'p', long, default_value = "output.fa")]
    pub prefix: String,

    /// Translate sequences to protein
    #[arg(short = 'X', long, default_value = "false", action = ArgAction::SetTrue)]
    pub translate: bool,

    /// Keep chunk outputs and skip merging into a single file
    #[arg(short = 'A', long = "as-chunk", default_value = "false", action = ArgAction::SetTrue)]
    pub as_chunk: bool,

    /// Also emit chunked BED outputs (requires --as-chunk)
    #[arg(short = 'B', long = "include-bed", requires = "as_chunk", default_value = "false", action = ArgAction::SetTrue)]
    pub include_bed: bool,

    /// Gzip-compress output files
    #[arg(short = 'Z', long, default_value = "false", action = ArgAction::SetTrue)]
    pub compress: bool,
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
        write!(
            f,
            "sequence={}, regions={}, outdir={}, chunks={}, upstream_flank={}, downstream_flank={}, feature={:?}, ignore_errors={}, level={}, prefix={}, translate={}, as_chunk={}, include_bed={}, compress={}",
            self.sequence.display(),
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
