//! get sequences from 2bit/fa using bed/gtf/gff
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This tool provides an easy way to get any sequence (exon, intron, cds, utr, etc.)
//! completely agnostic of the underlying format, either for reference sequences (2bit, fa, fa.gz)
//! or regions (bed, gtf, gff, gz, bz2, zstd)

pub const CODON_TABLE: [(&[u8], u8); 64] = [
    // Phenylalanine (F)
    (b"TTT", b'F'),
    (b"TTC", b'F'),
    // Leucine (L)
    (b"TTA", b'L'),
    (b"TTG", b'L'),
    (b"CTA", b'L'),
    (b"CTC", b'L'),
    (b"CTG", b'L'),
    (b"CTT", b'L'),
    // Isoleucine (I)
    (b"ATT", b'I'),
    (b"ATC", b'I'),
    (b"ATA", b'I'),
    // Methionine (M) - Start codon
    (b"ATG", b'M'),
    // Valine (V)
    (b"GTT", b'V'),
    (b"GTC", b'V'),
    (b"GTA", b'V'),
    (b"GTG", b'V'),
    // Serine (S)
    (b"TCT", b'S'),
    (b"TCC", b'S'),
    (b"TCA", b'S'),
    (b"AGT", b'S'),
    (b"AGC", b'S'),
    (b"TCG", b'S'),
    // Proline (P)
    (b"CCT", b'P'),
    (b"CCC", b'P'),
    (b"CCA", b'P'),
    (b"CCG", b'P'),
    // Threonine (T)
    (b"ACT", b'T'),
    (b"ACC", b'T'),
    (b"ACA", b'T'),
    (b"ACG", b'T'),
    // Alanine (A)
    (b"GCT", b'A'),
    (b"GCC", b'A'),
    (b"GCA", b'A'),
    (b"GCG", b'A'),
    // Tyrosine (Y)
    (b"TAT", b'Y'),
    (b"TAC", b'Y'),
    // Stop codons (*)
    (b"TAA", b'*'),
    (b"TAG", b'*'),
    (b"TGA", b'*'),
    // Histidine (H)
    (b"CAT", b'H'),
    (b"CAC", b'H'),
    // Glutamine (Q)
    (b"CAA", b'Q'),
    (b"CAG", b'Q'),
    // Asparagine (N)
    (b"AAT", b'N'),
    (b"AAC", b'N'),
    // Lysine (K)
    (b"AAA", b'K'),
    (b"AAG", b'K'),
    // Aspartic acid (D)
    (b"GAT", b'D'),
    (b"GAC", b'D'),
    // Glutamic acid (E)
    (b"GAA", b'E'),
    (b"GAG", b'E'),
    // Cysteine (C)
    (b"TGT", b'C'),
    (b"TGC", b'C'),
    // Tryptophan (W)
    (b"TGG", b'W'),
    // Arginine (R)
    (b"CGA", b'R'),
    (b"CGC", b'R'),
    (b"CGG", b'R'),
    (b"CGT", b'R'),
    (b"AGA", b'R'),
    (b"AGG", b'R'),
    // Glycine (G)
    (b"GGA", b'G'),
    (b"GGC", b'G'),
    (b"GGG", b'G'),
    (b"GGT", b'G'),
];
