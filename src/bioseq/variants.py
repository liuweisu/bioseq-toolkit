"""Detect large insertions (for example new TE insertions) in a sample sequence.

Medical relevance: new insertions of active human TEs (Alu, L1, SVA) are a
known cause of genetic disease, e.g. an Alu insertion in *NF1* causing
neurofibromatosis type 1 (Wallace et al. 1991, Nature) and L1 insertions in
the factor VIII gene causing haemophilia A (Kazazian et al. 1988, Nature).

Method: unique k-mers of the reference act as anchors. Between two
consecutive anchors the sample should have the same distance as the
reference; a jump means an indel. The breakpoint is then refined by
extending the left anchor base by base, and a TSD is looked for at the
right breakpoint.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class Insertion:
    ref_position: int  # 0-based position in the reference before which the insertion sits
    length: int  # inserted bases, including one TSD copy
    sequence: str  # inserted bases without the TSD
    tsd: str  # target site duplication ("" if none found)

    def in_interval(self, start: int, end: int) -> bool:
        return start <= self.ref_position < end


def _anchors(ref: str, alt: str, k: int) -> list[tuple[int, int]]:
    counts = Counter(ref[i : i + k] for i in range(len(ref) - k + 1))
    ref_pos = {ref[i : i + k]: i for i in range(len(ref) - k + 1) if counts[ref[i : i + k]] == 1}
    alt_counts = Counter(alt[i : i + k] for i in range(len(alt) - k + 1))
    chain: list[tuple[int, int]] = []
    for a in range(len(alt) - k + 1):
        kmer = alt[a : a + k]
        r = ref_pos.get(kmer)
        if r is None or alt_counts[kmer] != 1:
            continue
        if not chain or (r > chain[-1][0] and a > chain[-1][1]):
            chain.append((r, a))
    return chain


def find_insertions(ref: str, alt: str, k: int = 25, min_size: int = 50) -> list[Insertion]:
    ref, alt = ref.upper(), alt.upper()
    chain = _anchors(ref, alt, k)
    found: list[Insertion] = []
    for (r0, a0), (r1, a1) in zip(chain, chain[1:]):
        extra = (a1 - a0) - (r1 - r0)
        if extra < min_size:
            continue
        rb, ab = r0, a0
        while rb < len(ref) and ab < len(alt) and ref[rb] == alt[ab]:
            rb += 1
            ab += 1
        inserted = alt[ab : ab + extra]
        tsd = ""
        for n in range(25, 3, -1):
            if rb - n >= 0 and ref[rb - n : rb] == inserted[-n:]:
                tsd = ref[rb - n : rb]
                break
        body = inserted[: len(inserted) - len(tsd)] if tsd else inserted
        found.append(Insertion(rb - len(tsd), extra, body, tsd))
    return found


def frameshift(insertion_length: int) -> bool:
    return insertion_length % 3 != 0
