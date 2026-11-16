"""Pairwise alignment: Needleman-Wunsch (global) and Smith-Waterman (local).

Linear gap penalties, plain dynamic programming. Fine for genes, TE consensus
sequences and peptides; use minimap2/BLAST for whole genomes.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Alignment:
    score: int
    a: str  # aligned query with '-' for gaps
    b: str  # aligned target with '-' for gaps
    start_a: int  # 0-based start in the original query
    end_a: int  # exclusive end in the original query
    start_b: int
    end_b: int

    @property
    def length(self) -> int:
        return len(self.a)

    @property
    def identity(self) -> float:
        if not self.a:
            return 0.0
        return sum(x == y and x != "-" for x, y in zip(self.a, self.b)) / len(self.a)

    def pretty(self) -> str:
        mid = "".join("|" if x == y and x != "-" else " " for x, y in zip(self.a, self.b))
        return f"{self.a}\n{mid}\n{self.b}"


def _score(x: str, y: str, match: int, mismatch: int) -> int:
    return match if x == y else mismatch


def needleman_wunsch(
    a: str, b: str, match: int = 1, mismatch: int = -1, gap: int = -2
) -> Alignment:
    n, m = len(a), len(b)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        H[i][0] = i * gap
    for j in range(1, m + 1):
        H[0][j] = j * gap
    for i in range(1, n + 1):
        ai, row, prev = a[i - 1], H[i], H[i - 1]
        for j in range(1, m + 1):
            row[j] = max(
                prev[j - 1] + _score(ai, b[j - 1], match, mismatch),
                prev[j] + gap,
                row[j - 1] + gap,
            )
    al_a, al_b = _traceback(a, b, H, n, m, match, mismatch, gap, local=False)
    return Alignment(H[n][m], al_a, al_b, 0, n, 0, m)


def smith_waterman(
    a: str, b: str, match: int = 2, mismatch: int = -3, gap: int = -5
) -> Alignment:
    n, m = len(a), len(b)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        ai, row, prev = a[i - 1], H[i], H[i - 1]
        for j in range(1, m + 1):
            v = max(
                0,
                prev[j - 1] + _score(ai, b[j - 1], match, mismatch),
                prev[j] + gap,
                row[j - 1] + gap,
            )
            row[j] = v
            if v > best:
                best, bi, bj = v, i, j
    al_a, al_b, si, sj = _traceback(a, b, H, bi, bj, match, mismatch, gap, local=True)
    return Alignment(best, al_a, al_b, si, bi, sj, bj)


def _traceback(a, b, H, i, j, match, mismatch, gap, local):
    out_a: list[str] = []
    out_b: list[str] = []
    while i > 0 or j > 0:
        if local and H[i][j] == 0:
            break
        if i > 0 and j > 0 and H[i][j] == H[i - 1][j - 1] + _score(a[i - 1], b[j - 1], match, mismatch):
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and H[i][j] == H[i - 1][j] + gap:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    al_a, al_b = "".join(reversed(out_a)), "".join(reversed(out_b))
    if local:
        return al_a, al_b, i, j
    return al_a, al_b


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Hamming distance needs sequences of equal length")
    return sum(x != y for x, y in zip(a, b))
