"""A small de Bruijn graph genome assembler and assembly statistics.

The steps are the same as in real short-read assemblers (Velvet, SPAdes,
MEGAHIT), just without their many refinements:

1. Count k-mers in the reads and in their reverse complements.
2. Drop k-mers seen fewer than ``min_count`` times; most of them contain
   sequencing errors.
3. Build the graph: (k-1)-mers are nodes, k-mers are edges.
4. Walk maximal non-branching paths (unitigs); these are the contigs.
5. Each contig is found on both strands, so keep one copy of each pair.
"""

from __future__ import annotations

import random
from collections import Counter
from dataclasses import asdict, dataclass

from .seqio import Record, gc_content, reverse_complement

BASES = "ACGT"


def count_kmers(seqs: list[str], k: int, both_strands: bool = True) -> Counter[str]:
    counts: Counter[str] = Counter()
    for s in seqs:
        variants = (s, reverse_complement(s)) if both_strands else (s,)
        for v in variants:
            for i in range(len(v) - k + 1):
                kmer = v[i : i + k]
                if "N" not in kmer:
                    counts[kmer] += 1
    return counts


def kmer_spectrum(counts: Counter[str]) -> dict[int, int]:
    """Histogram: multiplicity -> number of distinct k-mers with that multiplicity.

    Error k-mers pile up at multiplicity 1-2; genomic k-mers form a peak near
    the k-mer coverage. The valley in between is a good ``min_count``.
    """
    return dict(sorted(Counter(counts.values()).items()))


class DeBruijnGraph:
    def __init__(self, kmers: set[str]):
        if not kmers:
            raise ValueError("No k-mers left to build a graph; lower min_count or k")
        self.kmers = kmers
        self.k = len(next(iter(kmers)))

    @classmethod
    def from_reads(cls, reads: list[str], k: int = 31, min_count: int = 2) -> DeBruijnGraph:
        counts = count_kmers(reads, k)
        return cls({km for km, c in counts.items() if c >= min_count})

    def out_edges(self, node: str) -> list[str]:
        return [node + b for b in BASES if node + b in self.kmers]

    def in_edges(self, node: str) -> list[str]:
        return [b + node for b in BASES if b + node in self.kmers]

    def _is_simple(self, node: str) -> bool:
        return len(self.in_edges(node)) == 1 and len(self.out_edges(node)) == 1

    def unitigs(self) -> list[tuple[str, bool]]:
        """Return all maximal non-branching paths as (sequence, is_cycle) pairs."""
        used: set[str] = set()
        paths: list[tuple[str, bool]] = []

        def extend(edge: str) -> str:
            seq = edge
            used.add(edge)
            node = edge[1:]
            while self._is_simple(node):
                nxt = self.out_edges(node)[0]
                if nxt in used:
                    break
                used.add(nxt)
                seq += nxt[-1]
                node = nxt[1:]
            return seq

        for edge in sorted(self.kmers):
            if edge in used:
                continue
            if not self._is_simple(edge[:-1]):
                paths.append((extend(edge), False))
        # Whatever is left belongs to cycles where every node is 1-in-1-out,
        # e.g. a circular bacterial chromosome or plasmid assembled end to end.
        # The walk returns to its first node, so drop the repeated k-1 bases.
        for edge in sorted(self.kmers):
            if edge not in used:
                walk = extend(edge)
                paths.append((walk[: len(walk) - (self.k - 1)], True))
        return paths


def assemble(
    reads: list[str], k: int = 31, min_count: int = 2, min_length: int = 0
) -> list[Record]:
    """Assemble reads into contigs, longest first."""
    graph = DeBruijnGraph.from_reads(reads, k=k, min_count=min_count)
    # Every contig is found twice, once per strand. Remember the k-mers of the
    # contigs we keep and skip a path whose first k-mer is already on the other
    # strand of a kept one (works for linear contigs and rotated cycles alike).
    kept_kmers: set[str] = set()
    contigs: list[tuple[str, bool]] = []
    for seq, circular in graph.unitigs():
        if len(seq) < min_length:
            continue
        if reverse_complement(seq[:k]) in kept_kmers or seq[:k] in kept_kmers:
            continue
        ring = seq + seq[: k - 1] if circular else seq
        kept_kmers.update(ring[i : i + k] for i in range(len(ring) - k + 1))
        contigs.append((seq, circular))
    contigs.sort(key=lambda c: len(c[0]), reverse=True)
    return [
        Record(f"contig_{i + 1}", c, f"len={len(c)} k={k}" + (" circular" if circ else ""))
        for i, (c, circ) in enumerate(contigs)
    ]


@dataclass(frozen=True)
class AssemblyStats:
    n_contigs: int
    total_length: int
    largest: int
    n50: int
    l50: int
    gc: float
    ng50: int | None = None

    def as_dict(self) -> dict:
        return asdict(self)


def assembly_stats(seqs: list[str], genome_size: int | None = None) -> AssemblyStats:
    """N50/L50 (and NG50 when the genome size is known), like QUAST reports.

    N50: the length L such that contigs of length >= L contain half the assembly.
    L50: how many contigs that takes.
    NG50: the same, but relative to the expected genome size, so it cannot be
    inflated by throwing away short contigs.
    """
    lengths = sorted((len(s) for s in seqs), reverse=True)
    total = sum(lengths)
    if not lengths:
        return AssemblyStats(0, 0, 0, 0, 0, 0.0, 0 if genome_size else None)

    def nx(target: float) -> tuple[int, int]:
        running = 0
        for i, length in enumerate(lengths, 1):
            running += length
            if running >= target:
                return length, i
        return 0, 0

    n50, l50 = nx(total / 2)
    ng50 = nx(genome_size / 2)[0] if genome_size else None
    return AssemblyStats(
        n_contigs=len(lengths),
        total_length=total,
        largest=lengths[0],
        n50=n50,
        l50=l50,
        gc=round(gc_content("".join(seqs)), 4),
        ng50=ng50,
    )


def random_genome(length: int, gc: float = 0.5, seed: int | None = None) -> str:
    rng = random.Random(seed)
    at, gcw = (1 - gc) / 2, gc / 2
    return "".join(rng.choices(BASES, weights=[at, gcw, gcw, at], k=length))


def simulate_reads(
    genome: str,
    read_length: int = 100,
    coverage: float = 30.0,
    error_rate: float = 0.005,
    seed: int | None = None,
    circular: bool = False,
) -> list[Record]:
    """Simulate Illumina-like single-end reads from both strands with substitution errors."""
    rng = random.Random(seed)
    template = genome + genome[: read_length - 1] if circular else genome
    n_starts = len(genome) if circular else len(genome) - read_length + 1
    if n_starts <= 0:
        raise ValueError("Genome is shorter than the read length")
    n_reads = int(coverage * len(genome) / read_length)
    reads: list[Record] = []
    for i in range(n_reads):
        start = rng.randrange(n_starts)
        frag = list(template[start : start + read_length])
        quals = []
        for j in range(len(frag)):
            if rng.random() < error_rate:
                frag[j] = rng.choice([b for b in BASES if b != frag[j]])
                quals.append("+")  # Q10
            else:
                quals.append("I")  # Q40
        seq = "".join(frag)
        if rng.random() < 0.5:
            seq, quals = reverse_complement(seq), quals[::-1]
        reads.append(Record(f"read_{i + 1}", seq, f"pos={start}", "".join(quals)))
    return reads
