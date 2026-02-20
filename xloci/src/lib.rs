//! get sequences from 2bit/fa using bed/gtf/gff
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This tool provides an easy way to get any sequence (exon, intron, cds, utr, etc.)
//! completely agnostic of the underlying format, either for reference sequences (2bit, fa, fa.gz)
//! or regions (bed, gtf, gff, gz, bz2, zstd)
//!
//! # Usage
//!
//! ```bash
//! Usage: xloci [OPTIONS] --sequence <SEQUENCE> --regions <REGIONS> --outdir <OUTDIR>
//!
//! Options:
//!   -s, --sequence <SEQUENCE>    Path to genome sequence file (.fa, .fa.gz, or .2bit)
//!   -r, --regions <REGIONS>      Path to genomic regions file (BED, GTF, or GFF format)
//!   -o, --outdir <OUTDIR>        Output directory for extracted sequences
//!   -c, --chunks <CHUNKS>        Number of records per parallel processing chunk [default: 1000]
//!   -u, --upstream-flank <UPSTREAM_FLANK>
//!                               Bases to extend upstream of features [default: 0]
//!   -d, --downstream-flank <DOWNSTREAM_FLANK>
//!                               Bases to extend downstream of features [default: 0]
//!   -f, --feature <FEATURE>      Type of genomic feature to extract [default: exon] [possible values: transcript, exon, intron, cds, utr]
//!   -I, --ignore-errors          Continue processing on errors instead of panicking
//!   -L, --level <LEVEL>          Logging verbosity level [default: info]  [possible values: trace, debug, info, warn, error]
//!   -p, --prefix <PREFIX>        Prefix for output files [default: output.fa]
//!   -X, --translate              Translate sequences to protein
//!   -A, --as-chunk               Keep chunk outputs and skip merging into a single file
//!   -B, --include-bed            Also emit chunked BED outputs (requires --as-chunk)
//!   -Z, --compress               Gzip-compress output files
//!   -h, --help                   Print help
//!   -V, --version                Print version
//! ```

pub mod cli;
pub mod consts;
pub mod core;

pub use cli::{Args, Feature};
pub use core::xloci;
