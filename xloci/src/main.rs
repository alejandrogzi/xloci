//! get sequences from 2bit/fa using bed/gtf/gff
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This tool provides an easy way to get any sequence (exon, intron, cds, utr, etc.)
//! completely agnostic of the underlying format, either for reference sequences (2bit, fa, fa.gz)
//! or regions (bed, gtf, gff, gz, bz2, zstd)

use clap::Parser;
use log::info;
use simple_logger::init_with_level;
use xloci::{Args, xloci};

/// Entry point
///
/// # Example
///
/// ```rust,ignore
/// // Run from command line:
/// // xloci -s genome.fa -r regions.gtf -o output/
/// fn main() {
///     let args = Args::parse();
///     xloci(args);
/// }
/// ```
fn main() {
    let args = Args::parse();

    init_with_level(args.level).unwrap_or_else(|e| panic!("{}", e));
    info!("Starting xloci with args: {}", args);

    xloci(args);
}
