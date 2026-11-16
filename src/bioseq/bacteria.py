"""Bacterial genome analysis: replication origin, genes and codon usage.

Most bacteria replicate their circular chromosome bidirectionally from one
origin (oriC) to the terminus (ter). The leading strand is G-rich because of
mutational bias during replication, so GC skew (G-C)/(G+C) changes sign at
ori and ter. The minimum of the *cumulative* skew marks oriC, the maximum ter
(Lobry 1996; Grigoriev 1998).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from .seqio import CODON_TABLE, reverse_complement, translate

#: E. coli DnaA box consensus; boxes cluster around oriC.
DNAA_BOX = "TTATCCACA"
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG", "GTG", "TTG"}


def gc_skew(seq: str, window: int = 1000, step: int | None = None) -> list[tuple[int, float]]:
    """(window midpoint, skew) pairs along the sequence."""
    seq = seq.upper()
    step = step or window
    out = []
    for start in range(0, max(1, len(seq) - window + 1), step):
        w = seq[start : start + window]
        g, c = w.count("G"), w.count("C")
        out.append((start + len(w) // 2, (g - c) / (g + c) if g + c else 0.0))
    return out


def cumulative_skew(seq: str) -> list[int]:
    """Running sum of +1 for G and -1 for C, sampled at every base."""
    total, out = 0, [0]
    for base in seq.upper():
        total += 1 if base == "G" else -1 if base == "C" else 0
        out.append(total)
    return out


@dataclass(frozen=True)
class OriPrediction:
    ori: int  # position of the cumulative-skew minimum
    ter: int  # position of the cumulative-skew maximum
    dnaa_boxes_near_ori: int
    dnaa_boxes_expected: float

    @property
    def enrichment(self) -> float:
        return self.dnaa_boxes_near_ori / self.dnaa_boxes_expected if self.dnaa_boxes_expected else 0.0


def find_motif(seq: str, motif: str, max_mismatches: int = 0, both_strands: bool = True) -> list[int]:
    seq, motif = seq.upper(), motif.upper()
    targets = {motif, reverse_complement(motif)} if both_strands else {motif}
    hits = []
    for i in range(len(seq) - len(motif) + 1):
        window = seq[i : i + len(motif)]
        if any(sum(a != b for a, b in zip(window, t)) <= max_mismatches for t in targets):
            hits.append(i)
    return hits


def predict_ori(seq: str, box: str = DNAA_BOX, window: int = 1000, max_mismatches: int = 1) -> OriPrediction:
    skew = cumulative_skew(seq)
    ori = min(range(len(skew)), key=skew.__getitem__)
    ter = max(range(len(skew)), key=skew.__getitem__)
    lo, hi = max(0, ori - window // 2), min(len(seq), ori + window // 2)
    near = len(find_motif(seq[lo:hi], box, max_mismatches))
    genome_wide = len(find_motif(seq, box, max_mismatches))
    expected = genome_wide * (hi - lo) / len(seq) if seq else 0.0
    return OriPrediction(ori, ter, near, round(expected, 3))


@dataclass(frozen=True)
class ORF:
    start: int  # 0-based, on the forward strand
    end: int  # exclusive, on the forward strand
    strand: str
    protein: str

    @property
    def length_nt(self) -> int:
        return self.end - self.start


def find_orfs(seq: str, min_protein_length: int = 100, starts: set[str] = START_CODONS) -> list[ORF]:
    """Open reading frames on all six frames, longest ORF per stop codon.

    This is how simple prokaryotic gene finders start; Prodigal adds RBS
    motifs and a trained coding model on top.
    """
    seq = seq.upper()
    n = len(seq)
    orfs: list[ORF] = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            start_pos: int | None = None
            for i in range(frame, n - 2, 3):
                codon = s[i : i + 3]
                if start_pos is None and codon in starts:
                    start_pos = i
                elif codon in STOP_CODONS:
                    if start_pos is not None and (i - start_pos) // 3 >= min_protein_length:
                        prot = "M" + translate(s[start_pos + 3 : i])
                        if strand == "+":
                            orfs.append(ORF(start_pos, i + 3, "+", prot))
                        else:
                            orfs.append(ORF(n - (i + 3), n - start_pos, "-", prot))
                    start_pos = None
    return sorted(orfs, key=lambda o: o.start)


def codon_usage(coding_seqs: list[str]) -> dict[str, float]:
    """Relative frequency of each codon within its amino acid (RSCU-style fractions)."""
    counts: Counter[str] = Counter()
    for cds in coding_seqs:
        cds = cds.upper()
        for i in range(0, len(cds) - 2, 3):
            codon = cds[i : i + 3]
            if codon in CODON_TABLE:
                counts[codon] += 1
    by_aa: Counter[str] = Counter()
    for codon, c in counts.items():
        by_aa[CODON_TABLE[codon]] += c
    return {
        codon: round(counts[codon] / by_aa[aa], 4) if by_aa[aa] else 0.0
        for codon, aa in sorted(CODON_TABLE.items())
    }


def coding_density(orfs: list[ORF], genome_length: int) -> float:
    covered = [False] * genome_length
    for o in orfs:
        for i in range(o.start, o.end):
            covered[i] = True
    return sum(covered) / genome_length if genome_length else 0.0
