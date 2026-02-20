<p align="center">
  <p align="center" style="margin-bottom: 0;">
    <img width=300 align="center" src="./assets/xloci.png" >
  </p>

  <span>
    <h1 align="center">
        xloci
    </h1>
  </span>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.0.1-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.0.1-green">
    </a>
    <a href="https://crates.io/crates/xloci" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/xloci">
    </a>
    <a href="https://github.com/alejandrogzi/xloci" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/xloci?color=blue">
    </a>
    <a href="https://crates.io/crates/xloci" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/xloci">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> get sequences from 2bit/fa using bed/gtf/gff</span>
        <br>
        <br>
        <a href="https://docs.rs/xloci/0.0.1/xloci/">docs</a> .
        <a href="https://github.com/alejandrogzi/xloci?tab=readme-ov-file#Usage">usage</a> .
        <a href="https://github.com/alejandrogzi/xloci?tab=readme-ov-file#Installation">install</a> .
        <a href="https://github.com/alejandrogzi/xloci/?tab=readme-ov-file#Conda">conda</a> .
    </samp>
  </p>

</p>

## Overview

This tool provides an easy way to get any sequence (exon, intron, cds, utr, etc.) completely agnostic of the underlying format, either for reference sequences (2bit, fa, fa.gz) or regions (bed, gtf, gff, gz, bz2, zstd).

## Quick Start

### Installation
to install xloci on your system follow this steps:
1. get rust: `curl https://sh.rustup.rs -sSf | sh` on unix, or go [here](https://www.rust-lang.org/tools/install) for other options
2. run `cargo install xloci` (make sure `~/.cargo/bin` is in your `$PATH` before running it)
4. use `xloci` with the required arguments
5. enjoy!

### Build
to build xloci from this repo, do:

1. get rust (as described above)
2. run `git clone https://github.com/alejandrogzi/xloci.git && cd xloci`
3. run `cargo run --release -- -i <GTF/GFF> -o <BED>`

### Container image
to build the development container image:
1. run `git clone https://github.com/alejandrogzi/xloci.git && cd xloci/assets`
2. initialize docker with `start docker` or `systemctl start docker`
3. build the image `docker image build --tag xloci .`
4. run `docker run --rm -v "[dir_where_your_gtf_is]:/dir" xloci -s /dir/<SEQUENCE> -r /dir/<REGIONS>`

### Conda
to use xloci through Conda just:
1. `conda install xloci -c bioconda` or `conda create -n xloci -c bioconda xloci`

## Usage

 ```bash
Usage: xloci [OPTIONS] --sequence <SEQUENCE> --regions <REGIONS> --outdir <OUTDIR>

Options:
  -s, --sequence <SEQUENCE>
          Path to genome sequence file (.fa, .fa.gz, or .2bit)
  -r, --regions <REGIONS>
          Path to genomic regions file (BED, GTF, or GFF format)
  -o, --outdir <OUTDIR>
          Output directory for extracted sequences
  -c, --chunks <CHUNKS>
          Number of records per parallel processing chunk [default: 1000]
  -u, --upstream-flank <UPSTREAM_FLANK>
          Bases to extend upstream of features [default: 0]
  -d, --downstream-flank <DOWNSTREAM_FLANK>
          Bases to extend downstream of features [default: 0]
  -f, --feature <FEATURE>
          Type of genomic feature to extract [default: exon] [possible values: transcript, exon, intron, cds, utr]
  -I, --ignore-errors
          Continue processing on errors instead of panicking
  -L, --level <LEVEL>
          Logging verbosity level [default: info]
  -p, --prefix <PREFIX>
          Prefix for output files [default: output.fa]
  -X, --translate
          Translate sequences to protein
  -A, --as-chunk
          Keep chunk outputs and skip merging into a single file
  -B, --include-bed
          Also emit chunked BED outputs (requires --as-chunk)
  -Z, --compress
          Gzip-compress output files
  -h, --help
          Print help
  -V, --version
          Print version
```
