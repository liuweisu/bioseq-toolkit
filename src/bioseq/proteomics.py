"""Protein chemistry and bottom-up proteomics.

In a bottom-up (shotgun) experiment proteins are digested with trypsin and
the peptide masses are measured by mass spectrometry. Peptide mass
fingerprinting (PMF) identifies a protein by matching the observed masses to
an *in silico* digest of every protein in a database.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Monoisotopic residue masses (Da)
MONO = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
# Average residue masses (Da)
AVERAGE = {
    "G": 57.0519, "A": 71.0788, "S": 87.0782, "P": 97.1167, "V": 99.1326,
    "T": 101.1051, "C": 103.1388, "L": 113.1594, "I": 113.1594, "N": 114.1038,
    "D": 115.0886, "Q": 128.1307, "K": 128.1741, "E": 129.1155, "M": 131.1926,
    "H": 137.1411, "F": 147.1766, "R": 156.1875, "Y": 163.1760, "W": 186.2132,
}
WATER_MONO, WATER_AVG = 18.010565, 18.01528
PROTON = 1.007276
CARBAMIDOMETHYL = 57.021464  # fixed modification of Cys after iodoacetamide

# pKa values (EMBOSS iep)
PKA_POS = {"Nterm": 8.6, "K": 10.8, "R": 12.5, "H": 6.5}
PKA_NEG = {"Cterm": 3.6, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}

# Kyte & Doolittle (1982) hydropathy
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

ENZYMES = {
    # cleave C-terminal to K/R unless followed by P
    "trypsin": re.compile(r"(?<=[KR])(?!P)"),
    # cleave C-terminal to F/W/Y/L unless followed by P (high specificity variant is F/W/Y only)
    "chymotrypsin": re.compile(r"(?<=[FWYL])(?!P)"),
    "lys-c": re.compile(r"(?<=K)"),
    "glu-c": re.compile(r"(?<=E)"),
}


def _check(protein: str) -> str:
    protein = protein.upper().replace("*", "")
    bad = set(protein) - set(MONO)
    if bad:
        raise ValueError(f"Unknown amino acid(s): {''.join(sorted(bad))}")
    return protein


def peptide_mass(peptide: str, monoisotopic: bool = True, carbamidomethyl: bool = True) -> float:
    """Neutral peptide mass [M]."""
    peptide = _check(peptide)
    table, water = (MONO, WATER_MONO) if monoisotopic else (AVERAGE, WATER_AVG)
    mass = sum(table[a] for a in peptide) + water
    if carbamidomethyl:
        mass += CARBAMIDOMETHYL * peptide.count("C")
    return mass


def mz(peptide: str, charge: int = 1, **kwargs) -> float:
    return (peptide_mass(peptide, **kwargs) + charge * PROTON) / charge


def molecular_weight(protein: str) -> float:
    return peptide_mass(protein, monoisotopic=False, carbamidomethyl=False)


def net_charge(protein: str, ph: float) -> float:
    protein = _check(protein)
    counts = {aa: protein.count(aa) for aa in "KRHDECY"}
    pos = 1 / (1 + 10 ** (ph - PKA_POS["Nterm"]))
    pos += sum(counts[a] / (1 + 10 ** (ph - PKA_POS[a])) for a in "KRH")
    neg = 1 / (1 + 10 ** (PKA_NEG["Cterm"] - ph))
    neg += sum(counts[a] / (1 + 10 ** (PKA_NEG[a] - ph)) for a in "DECY")
    return pos - neg


def isoelectric_point(protein: str, tol: float = 1e-4) -> float:
    """pH where the net charge is zero, found by bisection."""
    lo, hi = 0.0, 14.0
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if net_charge(protein, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def gravy(protein: str) -> float:
    """Grand average of hydropathy; > 0 suggests a hydrophobic (often membrane) protein."""
    protein = _check(protein)
    return round(sum(KYTE_DOOLITTLE[a] for a in protein) / len(protein), 3)


def hydropathy_profile(protein: str, window: int = 19) -> list[float]:
    """Sliding-window Kyte-Doolittle; windows > 1.6 hint at transmembrane helices."""
    protein = _check(protein)
    return [
        round(sum(KYTE_DOOLITTLE[a] for a in protein[i : i + window]) / window, 3)
        for i in range(len(protein) - window + 1)
    ]


@dataclass(frozen=True)
class Peptide:
    sequence: str
    start: int  # 0-based position in the protein
    missed_cleavages: int
    mass: float  # neutral monoisotopic mass

    @property
    def mh(self) -> float:
        """[M+H]+, the singly protonated ion seen in MALDI-TOF."""
        return self.mass + PROTON


def digest(
    protein: str,
    enzyme: str = "trypsin",
    missed_cleavages: int = 1,
    min_length: int = 6,
    max_length: int = 50,
) -> list[Peptide]:
    protein = _check(protein)
    pieces = [p for p in ENZYMES[enzyme].split(protein) if p]
    starts, pos = [], 0
    for p in pieces:
        starts.append(pos)
        pos += len(p)
    peptides = []
    for i in range(len(pieces)):
        for mc in range(missed_cleavages + 1):
            if i + mc >= len(pieces):
                break
            seq = "".join(pieces[i : i + mc + 1])
            if min_length <= len(seq) <= max_length:
                peptides.append(Peptide(seq, starts[i], mc, round(peptide_mass(seq), 5)))
    return peptides


@dataclass(frozen=True)
class PMFHit:
    protein_id: str
    matched: int
    coverage: float
    matched_peptides: list[str]
    score: float


def pmf_search(
    observed_mh: list[float],
    database: dict[str, str],
    tolerance_ppm: float = 20.0,
    missed_cleavages: int = 1,
) -> list[PMFHit]:
    """Rank database proteins by how many observed [M+H]+ peaks their digest explains.

    The score is the number of matches minus the number expected by chance,
    given how many theoretical peptides the protein has; this stops large
    proteins from winning just because they have more peptides.
    """
    observed = sorted(observed_mh)
    all_counts = {pid: digest(seq, missed_cleavages=missed_cleavages) for pid, seq in database.items()}
    span = (max(observed) - min(observed)) if len(observed) > 1 else 1.0
    hits = []
    for pid, peptides in all_counts.items():
        matched: dict[float, Peptide] = {}
        for pep in peptides:
            for m in observed:
                if abs(pep.mh - m) / m * 1e6 <= tolerance_ppm:
                    matched.setdefault(m, pep)
        covered = set()
        for pep in matched.values():
            covered.update(range(pep.start, pep.start + len(pep.sequence)))
        mean_window = sum(m * tolerance_ppm * 2e-6 for m in observed) / len(observed)
        p_random = min(1.0, len(peptides) * mean_window / span)
        expected = p_random * len(observed)
        hits.append(
            PMFHit(
                protein_id=pid,
                matched=len(matched),
                coverage=round(len(covered) / len(database[pid]), 3),
                matched_peptides=sorted(p.sequence for p in matched.values()),
                score=round(len(matched) - expected, 3),
            )
        )
    return sorted(hits, key=lambda h: (h.score, h.coverage), reverse=True)
