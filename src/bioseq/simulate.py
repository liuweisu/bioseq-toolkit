"""Simulated datasets with a known truth.

Real genomes are large and their annotations are themselves predictions, so
every method in this package is tested on data where the answer is known.
All sequences here are random and synthetic; family names such as
``AluLike`` describe the *structure* being imitated, not real consensus
sequences. For real data see the README.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from .assembly import random_genome
from .seqio import CODON_TABLE, Record, reverse_complement

BASES = "ACGT"
AA_FREQ = {  # approximate UniProt amino acid composition (%)
    "A": 8.3, "R": 5.5, "N": 4.1, "D": 5.5, "C": 1.4, "Q": 3.9, "E": 6.7, "G": 7.1, "H": 2.3,
    "I": 5.9, "L": 9.7, "K": 5.8, "M": 2.4, "F": 3.9, "P": 4.7, "S": 6.6, "T": 5.4, "W": 1.1,
    "Y": 2.9, "V": 6.9,
}


def mutate(seq: str, divergence: float, rng: random.Random) -> str:
    """Apply point substitutions at the given per-site rate."""
    return "".join(
        rng.choice([b for b in BASES if b != c]) if rng.random() < divergence else c for c in seq
    )


@dataclass
class TruthTE:
    family: str
    te_class: str
    start: int
    end: int
    strand: str
    tsd: str
    ltr_divergence: float | None = None


@dataclass
class TEGenome:
    sequence: str
    truth: list[TruthTE]
    library: list[Record] = field(default_factory=list)


def _ltr_family(name: str, superfamily: str, rng: random.Random, ltr_len: int, internal_len: int):
    ltr = "TG" + random_genome(ltr_len - 4, 0.40, rng.random()) + "CA"
    internal = random_genome(internal_len, 0.42, rng.random())
    return Record(f"{name}#LTR/{superfamily}", ltr + internal + ltr), ltr, internal


def _tir_family(name: str, superfamily: str, rng: random.Random, tir_len: int, internal_len: int):
    tir = random_genome(tir_len, 0.45, rng.random())
    seq = tir + random_genome(internal_len, 0.38, rng.random()) + reverse_complement(tir)
    return Record(f"{name}#DNA/{superfamily}", seq)


def plant_genome_with_tes(length: int = 80_000, seed: int = 7) -> TEGenome:
    """A plant-like genome fragment: GC-poor background with LTR retrotransposons and TIR elements.

    LTR copies get a 5' and 3' LTR that diverge in proportion to a random
    insertion age, which is what LTR dating relies on.
    """
    rng = random.Random(seed)
    families = []
    for name, sf, lt, il in (("RLC_Tekay", "Copia", 450, 3500), ("RLG_Retand", "Gypsy", 700, 4500)):
        rec, ltr, internal = _ltr_family(name, sf, rng, lt, il)
        families.append(("LTR", rec, ltr, internal))
    for name, sf, tl, il in (("DTT_Stowaway", "TcMar-Stowaway", 30, 250), ("DTA_hAT1", "hAT", 18, 2200)):
        families.append(("TIR", _tir_family(name, sf, rng, tl, il), None, None))

    background = random_genome(length, 0.36, rng.random())
    plan = []  # (family index, divergence)
    for _ in range(4):
        plan.append((0, rng.uniform(0.0, 0.04)))
        plan.append((1, rng.uniform(0.0, 0.04)))
    for _ in range(5):
        plan.append((2, rng.uniform(0.02, 0.08)))
    for _ in range(3):
        plan.append((3, rng.uniform(0.02, 0.08)))
    rng.shuffle(plan)

    slots = sorted(rng.sample(range(500, length - 500, 100), len(plan)))
    pieces: list[str] = []
    truth: list[TruthTE] = []
    cursor = 0
    out_len = 0
    for (fam_idx, div), pos in zip(plan, slots):
        kind, rec, ltr, internal = families[fam_idx]
        pieces.append(background[cursor:pos])
        out_len += pos - cursor
        family, _, te_class = rec.id.partition("#")
        if kind == "LTR":
            ltr5 = mutate(ltr, div / 2, rng)
            ltr3 = mutate(ltr, div / 2, rng)
            element = ltr5 + mutate(internal, 0.03, rng) + ltr3
            tsd_len = 5
            ltr_div = sum(a != b for a, b in zip(ltr5, ltr3)) / len(ltr)
        else:
            element = mutate(rec.seq, div, rng)
            tsd_len = 2 if "Stowaway" in family else 8
            ltr_div = None
        strand = rng.choice("+-")
        if strand == "-":
            element = reverse_complement(element)
        tsd = background[pos - tsd_len : pos] if kind == "LTR" or tsd_len == 8 else "TA"
        if tsd == "TA":
            pieces.append("TA")
            out_len += 2
        start = out_len
        pieces.append(element)
        out_len += len(element)
        truth.append(TruthTE(family, te_class, start, out_len, strand, tsd, ltr_div))
        pieces.append(tsd)
        out_len += len(tsd)
        cursor = pos
    pieces.append(background[cursor:])
    return TEGenome("".join(pieces), truth, [f[1] for f in families])


@dataclass
class GeneCase:
    reference: str
    patient: str
    exons: list[tuple[int, int]]
    library: list[Record]
    true_position: int
    true_family: str


def human_gene_with_te_insertion(seed: int = 11) -> GeneCase:
    """A reference gene region and a patient allele with a new SINE (Alu-like) insertion in exon 3."""
    rng = random.Random(seed)
    ref = random_genome(6000, 0.41, rng.random())
    exons = [(400, 560), (1500, 1680), (2800, 3010), (4100, 4300), (5200, 5450)]

    # Alu-like: ~280 bp, two GC-rich monomers joined by an A-rich linker, poly(A) tail
    mono = random_genome(120, 0.62, rng.random())
    alu_like = mono + "AAAAATACAAAAA" + mutate(mono, 0.15, rng) + random_genome(31, 0.6, rng.random()) + "A" * 22
    l1_like = random_genome(1800, 0.42, rng.random()) + "A" * 18
    sva_like = ("CCCTCT" * 10) + random_genome(900, 0.6, rng.random())
    library = [
        Record("AluLike#SINE/Alu", alu_like),
        Record("L1Like#LINE/L1", l1_like),
        Record("SVALike#Retroposon/SVA", sva_like),
    ]

    site = 2900  # inside exon 3
    tsd = ref[site - 12 : site]
    insert = mutate(alu_like, 0.03, rng)
    patient = ref[:site] + insert + tsd + ref[site:]
    return GeneCase(ref, patient, exons, library, site - 12, "AluLike")


@dataclass
class BacterialGenome:
    sequence: str
    ori: int
    ter: int
    genes: list[tuple[int, int, str]]  # start, end, strand


def _random_cds(n_codons: int, rng: random.Random, skew: float = 0.0) -> str:
    """Random coding sequence; ``skew`` > 0 favours G over C on the coding strand."""
    sense = [c for c, aa in CODON_TABLE.items() if aa != "*"]
    weights = [math.exp(skew * (c.count("G") - c.count("C"))) for c in sense]
    body = "".join(rng.choices(sense, weights=weights, k=n_codons - 2))
    return "ATG" + body + rng.choice(["TAA", "TAG", "TGA"])


def bacterial_genome(length: int = 120_000, ori: int = 30_000, seed: int = 3) -> BacterialGenome:
    """A circular chromosome with replication-associated GC skew, DnaA boxes at oriC and planted genes."""
    rng = random.Random(seed)
    ter = (ori + length // 2) % length
    seq = []
    for i in range(length):
        # leading strand of the forward-strand replichore (ori -> ter clockwise) is G-rich
        clockwise = (i - ori) % length < length // 2
        g, c = (0.29, 0.21) if clockwise else (0.21, 0.29)
        seq.append(rng.choices(BASES, weights=[0.25, c, g, 0.25])[0])
    genome = "".join(seq)
    # DnaA box cluster around oriC
    for offset in (-300, -180, -60, 80, 210):
        p = ori + offset
        box = "TTATCCACA" if offset < 0 else reverse_complement("TTATCCACA")
        genome = genome[:p] + box + genome[p + 9 :]
    genes = []
    pos = 1000
    while pos < length - 4000:
        n_codons = rng.randint(150, 900)
        strand = rng.choice("+-")
        # genes share the replication skew of the replichore they sit in
        forward_sign = 1 if (pos - ori) % length < length // 2 else -1
        cds = _random_cds(n_codons, rng, skew=0.3 * forward_sign * (1 if strand == "+" else -1))
        if abs(pos - ori) < 1500 or pos + len(cds) > length - 1000:
            pos += 2000
            continue
        if strand == "-":
            cds = reverse_complement(cds)
        genome = genome[:pos] + cds + genome[pos + len(cds) :]
        genes.append((pos, pos + len(cds), strand))
        pos += len(cds) + rng.randint(50, 400)
    return BacterialGenome(genome, ori, ter, genes)


def random_protein(length: int, rng: random.Random) -> str:
    aas, weights = zip(*AA_FREQ.items())
    return "M" + "".join(rng.choices(aas, weights=weights, k=length - 1))


@dataclass
class PMFDataset:
    database: dict[str, str]
    observed_mh: list[float]
    true_protein: str


def pmf_dataset(n_proteins: int = 200, seed: int = 5, noise_peaks: int = 15) -> PMFDataset:
    """A protein database and a noisy MALDI peak list from one of its proteins."""
    from .proteomics import digest  # local import avoids a cycle at import time

    rng = random.Random(seed)
    database = {f"PROT{i:04d}": random_protein(rng.randint(150, 700), rng) for i in range(n_proteins)}
    target = rng.choice(sorted(database))
    peptides = [p for p in digest(database[target], missed_cleavages=0) if 700 < p.mh < 3500]
    seen = rng.sample(peptides, max(3, int(0.6 * len(peptides))))
    observed = [p.mh * (1 + rng.gauss(0, 5e-6)) for p in seen]
    observed += [rng.uniform(700, 3500) for _ in range(noise_peaks)]
    return PMFDataset(database, sorted(round(m, 4) for m in observed), target)
