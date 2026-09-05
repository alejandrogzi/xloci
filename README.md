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
    <a href="https://img.shields.io/badge/version-0.0.6-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.0.6-green">
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
        <a href="https://docs.rs/xloci/0.0.6/xloci/">docs</a> .
        <a href="https://github.com/alejandrogzi/xloci?tab=readme-ov-file#Usage">usage</a> .
        <a href="https://github.com/alejandrogzi/xloci?tab=readme-ov-file#Installation">install</a> .
        <a href="https://github.com/alejandrogzi/xloci/?tab=readme-ov-file#Conda">conda</a>
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

### Nextflow
to use xloci through Nextflow just:
1. `nextflow run alejandrogzi/xloci -r <REGIONS> -s <SEQUENCE>`
or borrow the [xloci.nf](https://github.com/alejandrogzi/xloci/blob/main/assets/nextflow/xloci/main.nf) file from this repo

## Usage

 ```plaintext
Usage: xloci [OPTIONS] --regions <REGIONS> --outdir <OUTDIR>

Options:
  -s, --sequence <SEQUENCE>
          Path to genome sequence file (.fa, .fa.gz, or .2bit); reads from stdin when omitted
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
          Skip failing features with a warning instead of panicking
  -L, --level <LEVEL>
          Logging verbosity level [default: info]
  -p, --prefix <PREFIX>
          Stem for output FASTA files (writes <prefix>.fa or <prefix>.fa.gz) [default: output]
  -X, --translate
          Translate sequences to protein
  -U, --unmask
          Convert soft-masked bases to uppercase in output
  -S, --split-extraction
          Emit one output record per extracted feature piece
      --as-tsv
          Write tab-separated output instead of FASTA
      --add-tab
          Separate flank columns in TSV output (requires --as-tsv and at least one flank)
  -G, --generic-id
          Use genomic coordinates as identifiers instead of record names
  -A, --as-chunk
          Keep chunk outputs and skip merging into a single file
  -B, --include-bed
          Also emit chunked BED outputs (requires --as-chunk)
  -Z, --compress
          Gzip-compress output files
  -t, --threads <THREADS>
          Number of threads [default: 16]
  -h, --help
          Print help
  -V, --version
          Print version
```

## Benchmarks

### Feature comparison

| Tool | BED12 support | GTF/GFF support | FASTA support | .2bit support | spliced/exon extraction | strand-aware | notes/limitations | install source | link |
|---|---|---|---|---|---|---|---|---|---|
| xloci | ✅ | ✅ | ✅ | ✅ | yes (-f exon) | yes (default RC on minus) | genepred overhead | cargo/git/docker/bioconda(*) | this repo |
| bedtools getfasta | ✅ | ❌ | ✅ | ❌ | yes (BED12 blocks via -split) | yes (-s) | .2bit not supported; splicing semantics are BED12-block based | Bioconda; upstream repo | [1](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html) |
| gffread | ✅ | ✅ | ✅ | ❌ | yes (-w spliced exons; -x spliced CDS) | unspecified | Documented speedup with FASTA .fai; .2bit not documented | binaries/source/GitHub (official page) | [21](https://ccb.jhu.edu/software/stringtie/gff.shtml) |
| agat_sp_extract_sequences.pl | ❌ | ✅ | ✅ | ❌ | yes (-t exon --merge, --mrna, etc.) | yes (default RC on minus; controls available) | no BED input; .2bit not documented | unspecified | [13](https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html) |
| UCSC twoBitToFa | ✅ | ❌ | ❌ | ✅ | yes (“exclude introns” from BED blocks) | yes (RC on - strand) | Designed for .2bit; for GTF/GFF you must convert to BED first | Bioconda recipe; UCSC docs | [22](https://github.com/ENCODE-DCC/kentUtils/blob/master/src/utils/twoBitToFa/twoBitToFa.c) |
| GenomeTools gt extractfeat | ❌ | ✅ | ✅ | ❌ | join support yes (-join) | unspecified | Uses GFF3 graphs; strand behavior not stated in the short manual excerpt | unspecified in cited excerpts | [16](https://genometools.org/tools/gt_extractfeat.html) |
| gff3_to_fasta (GFF3 Toolkit) | ❌ | ✅ | ✅ | ❌ | yes (-st trans spliced transcripts) | unspecified | GFF3-only per docs; swiss-army script with multiple sequence types | Python script; packaging unspecified | [19](https://gff3toolkit.readthedocs.io/en/latest/gff3_to_fasta.html) |
| TopHat gtf_to_fasta | ❌ | ✅ | ✅ | ❌ | yes (exon concatenation) | no (inferred from shown source excerpt) | Legacy; strand/orientation behavior is not presented as a documented feature here | unspecified | [20](https://github.com/DaehwanKimLab/tophat/blob/master/src/GTFToFasta.cpp) |

### Runtime comparison

> [!NOTE]
> Benchamrk was done with hyperfine on a AMD Ryzen 7 5700X with 128 GB of RAM and 16 cores.
> AGAT was excluded from FASTA + GTF because of extremely long runtimes (over 10 minutes) and poor performance.
> GFF3 Toolkit was not included because of problems with the installation.
> GenomeTools was excluded because of problems with the installation.

#### 2bit + BED

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `xloci -s tmp/hg38.2bit -o output -r tmp/gencode.v44.annotation.bed` | 5.567 ± 0.028 | 5.538 | 5.593 | 1.00 |
| `twoBitToFa -bed=tmp/gencode.v44.annotation.bed tmp/hg38.2bit output.fa` | 58.787 ± 1.001 | 58.075 | 59.932 | 10.56 ± 0.19 |

#### FASTA + BED

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `xloci -s tmp/hg38.fa -o output -r tmp/gencode.v44.annotation.bed` | 4.164 ± 0.012 | 4.152 | 4.176 | 1.00 |
| `bedtools getfasta -fi tmp/hg38.fa -bed tmp/gencode.v44.annotation.bed -split -name -fo output.fa` | 12.167 ± 0.116 | 12.047 | 12.277 | 2.92 ± 0.03 |
| `gffread -w output.fa -g tmp/hg38.fa --in-bed tmp/gencode.v44.annotation.bed` | 6.550 ± 0.101 | 6.485 | 6.666 | 1.57 ± 0.02 |
| `bed2gtf -b tmp/gencode.v44.annotation.bed -o transcripts.gtf -n && agat_sp_extract_sequences.pl -g transcripts.gtf -f tmp/hg38.fa -t exon --merge --cpu 16` | 918 ± 0.151 | 917.123 | 923.312 | 220.67 ± 0.02 |

#### 2bit + GTF

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `xloci -s tmp/hg38.2bit -o output -r tmp/gencode.v44.annotation.gtf` | 10.835 ± 0.011 | 10.827 | 10.847 | 2.08 ± 0.00 |
| `gxf2bed -i tmp/gencode.v44.annotation.gtf -o transcripts.bed && xloci -s tmp/hg38.2bit -o output -r transcripts.bed` | 10.889 ± 0.037 | 10.863 | 10.931 | 2.09 ± 0.01 |
| `gffread tmp/gencode.v44.annotation.gtf --bed -o transcripts.bed && cat -p --color never transcripts.bed \| choose :11  -o '\t' > tmp.bed && xloci -s tmp/hg38.2bit -o output -r tmp.bed` | 9.080 ± 0.031 | 9.045 | 9.103 | 1.74 ± 0.01 |
| `gffread tmp/gencode.v44.annotation.gtf --bed -o transcripts.bed && cat -p --color never transcripts.bed \| choose :11  -o '\t' > tmp.bed && twoBitToFa -bed=tmp.bed tmp/hg38.2bit output.fa` | 5.215 ± 0.008 | 5.209 | 5.224 | 1.00 |
| `gxf2bed -i tmp/gencode.v44.annotation.gtf -o transcripts.bed && sort -k1,1 -k2,2n -k3,3n transcripts.bed > tmp.bed && twoBitToFa -bed=tmp.bed tmp/hg38.2bit output.fa` | 11.061 ± 0.022 | 11.048 | 11.087 | 2.12 ± 0.01 |

#### FASTA + GTF

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `xloci -s tmp/hg38.fa -o output -r tmp/gencode.v44.annotation.gtf` | 9.417 ± 0.024 | 9.403 | 9.445 | 1.80 ± 0.01 |
| `gffread tmp/gencode.v44.annotation.gtf --bed -o transcripts.bed && cat -p --color never transcripts.bed \| choose :11 -o '\t' > tmp.bed && bedtools getfasta -fi tmp/hg38.fa -bed tmp.bed -split -name -fo output.fa` | 5.230 ± 0.008 | 5.221 | 5.237 | 1.00 |
| `gxf2bed -i tmp/gencode.v44.annotation.gtf -o transcripts.bed && bedtools getfasta -fi tmp/hg38.fa -bed transcripts.bed -split -name -fo output.fa` | 17.477 ± 0.045 | 17.436 | 17.525 | 3.34 ± 0.01 |
| `gffread -w output.fa -g tmp/hg38.fa tmp/gencode.v44.annotation.gtf` | 10.368 ± 0.057 | 10.322 | 10.432 | 1.98 ± 0.01 |
| `gtf_to_fasta tmp/gencode.v44.annotation.gtf tmp/hg38.fa output.fa` | 40.49 ± 0.024 | 40.01 | 40.76 | 7.74 ± 0.01 |
