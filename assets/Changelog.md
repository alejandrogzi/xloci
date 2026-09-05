# Changelog

## v0.0.6 - 2026-09-05

- Changed `-I, --ignore-errors` so split extraction skips failing pieces with a warning instead of dropping the whole record.
- Bumped `genepred` to `0.0.16`.
- Bumped the crate version to `0.0.6`.

## v0.0.5 - 2026-06-09

- Added `-U, --unmask` to convert soft-masked output sequence bases to uppercase.
- Added this changelog under `assets/`.
- Bumped the crate version to `0.0.5`.

## v0.0.4 - 2026-03-25

- Added TSV output with `--as-tsv`.
- Added split extraction with `--split-extraction`.
- Added generic coordinate identifiers with `--generic-id`.
- Added separated flank columns with `--add-tab`.

## v0.0.3 - 2026-03-18

- Added stdin sequence input support.
- Added the Nextflow module under `assets/nextflow/`.
- Added benchmark coverage.
- Changed `--prefix` handling to use the provided value as an output stem.
- Updated container support for Nextflow compatibility.

## v0.0.2 - 2026-03-08

- Fixed FASTA input extension handling to cover supported FASTA suffixes.
- Added the project license.

## v0.0.1 - 2026-02-20

- Released the initial `xloci` implementation.
- Added project documentation and assets.
- Added Dockerfile support.
- Updated the Rust toolchain version used by the project.
