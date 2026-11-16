"""Transposable element (TE) discovery and annotation.

Three structural signals are used, the same ones that tools such as
LTR_FINDER, LTRharvest, EDTA and EarlGrey build on:

* **LTR retrotransposons** (most of a typical plant genome): two long,
  near-identical *direct* repeats a few kb apart, usually starting with TG and
  ending with CA, flanked by a short target site duplication (TSD, often 5 bp).
  Because both LTRs are identical when the element inserts, their divergence
  is a molecular clock for the insertion age.
* **DNA transposons (TIR elements)**: short *inverted* repeats at both ends,
  flanked by a TSD (TA for Tc1/mariner, 8 bp for hAT, ...).
* **High-copy k-mers**: any interspersed repeat, regardless of structure,
  shows up as k-mers that occur far more often than the single-copy genome.

Everything here is seeded with exact k-mer matches and then extended along
diagonals, which keeps it fast enough for a few Mb in pure Python.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass

from .align import smith_waterman
from .seqio import Record, reverse_complement

#: Neutral substitution rate used for plant LTR dating (rice; Ma & Bennetzen 2004).
PLANT_MUTATION_RATE = 1.3e-8


@dataclass(frozen=True)
class LTRCandidate:
    start: int  # 0-based start of the 5' LTR
    end: int  # exclusive end of the 3' LTR
    ltr_length: int
    identity: float
    tsd: str  # target site duplication if one was found, else ""
    has_tg_ca: bool

    @property
    def length(self) -> int:
        return self.end - self.start

    def insertion_age(self, mutation_rate: float = PLANT_MUTATION_RATE) -> float:
        """Estimated years since insertion: T = K / (2 * r).

        K is the divergence between the two LTRs (Jukes-Cantor corrected) and r
        the substitution rate per site per year. Use :func:`kimura_2p` on the
        aligned LTRs when you have them.
        """
        k = jukes_cantor(1 - self.identity)
        return k / (2 * mutation_rate)


@dataclass(frozen=True)
class TIRCandidate:
    start: int
    end: int
    tir_length: int
    identity: float
    tsd: str

    @property
    def length(self) -> int:
        return self.end - self.start


def jukes_cantor(p: float) -> float:
    """Correct an observed p-distance for multiple hits (Jukes-Cantor)."""
    if p >= 0.75:
        return math.inf
    return -0.75 * math.log(1 - 4 * p / 3)


def kimura_2p(a: str, b: str) -> float:
    """Kimura 2-parameter distance between two aligned, gap-free sequences."""
    purines = set("AG")
    transitions = transversions = n = 0
    for x, y in zip(a.upper(), b.upper()):
        if x not in "ACGT" or y not in "ACGT":
            continue
        n += 1
        if x != y:
            if (x in purines) == (y in purines):
                transitions += 1
            else:
                transversions += 1
    if n == 0:
        return 0.0
    P, Q = transitions / n, transversions / n
    try:
        return -0.5 * math.log((1 - 2 * P - Q) * math.sqrt(1 - 2 * Q))
    except ValueError:
        return math.inf


def _kmer_index(seq: str, k: int) -> dict[str, list[int]]:
    index: dict[str, list[int]] = defaultdict(list)
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" not in kmer:
            index[kmer].append(i)
    return index


def _merge_hits(points: list[int], max_gap: int) -> list[tuple[int, int]]:
    """Merge sorted seed start positions into (first, last) runs."""
    runs: list[tuple[int, int]] = []
    for p in sorted(points):
        if runs and p - runs[-1][1] <= max_gap:
            runs[-1] = (runs[-1][0], p)
        else:
            runs.append((p, p))
    return runs


def _identity(a: str, b: str) -> float:
    if not a:
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / len(a)


_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _extend(pairs, x_drop: int = 8) -> int:
    """Extend an alignment one base at a time; return how many steps to keep.

    ``pairs`` yields (base_a, base_b) moving away from the seed. Matches score
    +1, mismatches -2; we stop once the score falls ``x_drop`` below the best
    seen and keep the extension up to that best point (BLAST-style X-drop).
    """
    score = best = best_n = 0
    for n, (x, y) in enumerate(pairs, 1):
        score += 1 if x == y else -2
        if score > best:
            best, best_n = score, n
        elif best - score > x_drop:
            break
    return best_n


def _snap_ltr(seq: str, start: int, end: int, slack: int = 3) -> tuple[int, int]:
    """Move fuzzy LTR boundaries (by up to ``slack`` bp) onto TG...CA with a TSD.

    X-drop extension stops at a mismatch or runs on through chance matches in
    the flanks, so ends are often off by a base or two. Real element ends are
    recognisable, so look for them nearby; keep the input if nothing fits.
    """
    best = None
    for ds in range(-slack, slack + 1):
        for de in range(-slack, slack + 1):
            a, b = start + ds, end + de
            if seq[a : a + 2] == "TG" and seq[b - 2 : b] == "CA" and _find_tsd(seq, a, b):
                cost = abs(ds) + abs(de)
                if best is None or cost < best[0]:
                    best = (cost, a, b)
    return (best[1], best[2]) if best else (start, end)


def _snap_tir(seq: str, start: int, end: int, slack: int = 5) -> tuple[int, str]:
    """Trim both ends symmetrically to the position with the longest TSD.

    A 'TA' TSD next to a TIR that starts with TA (Tc1/mariner, Stowaway MITEs)
    looks like part of the inverted repeat, so extension runs through it.
    """
    best = (0, _find_tsd(seq, start, end, range(10, 1, -1)))
    for t in range(-slack, slack + 1):
        tsd = _find_tsd(seq, start + t, end - t, range(10, 1, -1))
        rank = (len(tsd) + (0.5 if tsd == "TA" else 0), -abs(t))
        best_rank = (len(best[1]) + (0.5 if best[1] == "TA" else 0), -abs(best[0]))
        if tsd and rank > best_rank:
            best = (t, tsd)
    return best


def _find_tsd(seq: str, left: int, right: int, sizes=range(10, 3, -1)) -> str:
    """Look for identical short sequences just outside [left, right)."""
    for n in sizes:
        if left - n < 0 or right + n > len(seq):
            continue
        a, b = seq[left - n : left], seq[right : right + n]
        if a == b:
            return a
    return ""


def find_ltr_retrotransposons(
    seq: str,
    k: int = 16,
    min_ltr: int = 100,
    max_ltr: int = 3000,
    min_distance: int = 1000,
    max_distance: int = 20000,
    min_identity: float = 0.80,
    max_kmer_copies: int = 50,
    strict: bool = True,
) -> list[LTRCandidate]:
    """Find pairs of long direct repeats with LTR-retrotransposon structure.

    With ``strict`` (the default) a candidate needs both the TG...CA termini and
    a target site duplication. Without it, pairs of unrelated copies of the
    same family (e.g. two neighbouring MITEs) are reported too; that is the
    main source of false positives in real LTR searches as well.
    """
    seq = seq.upper()
    index = _kmer_index(seq, k)
    diagonals: dict[int, list[int]] = defaultdict(list)
    for positions in index.values():
        if len(positions) < 2 or len(positions) > max_kmer_copies:
            continue
        for x in range(len(positions)):
            for y in range(x + 1, len(positions)):
                d = positions[y] - positions[x]
                if min_distance <= d <= max_distance:
                    diagonals[d].append(positions[x])

    raw: list[LTRCandidate] = []
    for d, starts in diagonals.items():
        for first, last in _merge_hits(starts, max_gap=60):
            end5 = last + k
            # refine both ends of the direct repeat along its diagonal
            left_ext = _extend((seq[i], seq[i + d]) for i in range(first - 1, -1, -1))
            right_ext = _extend(
                (seq[i], seq[i + d]) for i in range(end5, min(first + d, len(seq) - d))
            )
            first, end5 = first - left_ext, end5 + right_ext
            ltr_len = end5 - first
            if not (min_ltr <= ltr_len <= max_ltr) or end5 > first + d:
                continue
            left = seq[first:end5]
            right = seq[first + d : end5 + d]
            ident = _identity(left, right)
            if ident < min_identity:
                continue
            start, end = _snap_ltr(seq, first, end5 + d)
            cand = LTRCandidate(
                start=start,
                end=end,
                ltr_length=end5 - start,
                identity=round(ident, 4),
                tsd=_find_tsd(seq, start, end),
                has_tg_ca=seq[start : start + 2] == "TG" and seq[end - 2 : end] == "CA",
            )
            if strict and not (cand.tsd and cand.has_tg_ca):
                continue
            raw.append(cand)
    return _remove_overlaps(raw, key=lambda c: c.ltr_length * c.identity)


def find_tir_transposons(
    seq: str,
    k: int = 12,
    min_tir: int = 15,
    max_tir: int = 200,
    min_length: int = 200,
    max_length: int = 6000,
    min_identity: float = 0.80,
    max_kmer_copies: int = 50,
    strict: bool = True,
) -> list[TIRCandidate]:
    """Find elements bounded by terminal inverted repeats.

    Short inverted repeats are common by chance, so with ``strict`` a
    candidate must also have a TSD (2-10 bp). ``max_tir`` stops two copies of
    a family in opposite orientation from being called one element whose
    "TIRs" are the copies themselves; most real TIRs are well under 200 bp.
    """
    seq = seq.upper()
    index = _kmer_index(seq, k)
    anti: dict[int, list[int]] = defaultdict(list)
    for kmer, left_positions in index.items():
        if len(left_positions) > max_kmer_copies:
            continue
        right_positions = index.get(reverse_complement(kmer), [])
        if not right_positions or len(right_positions) > max_kmer_copies:
            continue
        for i in left_positions:
            for j in right_positions:
                # left arm seed at i, right arm seed at j (forward coords) with j > i.
                # Along one inverted repeat i + j stays constant.
                span = j + k - i
                if j > i and min_length <= span <= max_length:
                    anti[i + j].append(i)

    raw: list[TIRCandidate] = []
    for s_, starts in anti.items():
        for first, last in _merge_hits(starts, max_gap=8):
            start = first
            end = s_ - first + k  # right arm end (exclusive)
            inner_l = last + k  # left arm end (exclusive)
            inner_r = s_ - last  # right arm start
            # extend outward: base at start-1 pairs with the complement of base at end
            out = _extend(
                (seq[start - 1 - n], _COMP.get(seq[end + n], "N"))
                for n in range(min(start, len(seq) - end))
            )
            # extend inward
            inn = _extend(
                (seq[inner_l + n], _COMP.get(seq[inner_r - 1 - n], "N"))
                for n in range(max(0, (inner_r - inner_l) // 2))
            )
            start, end = start - out, end + out
            tir_len = inner_l + inn - start
            if not (min_length <= end - start <= max_length) or not (min_tir <= tir_len <= max_tir):
                continue
            left = seq[start : start + tir_len]
            right = reverse_complement(seq[end - tir_len : end])
            ident = _identity(left, right)
            if ident < min_identity:
                continue
            trim, tsd = _snap_tir(seq, start, end)
            if tir_len - trim < min_tir:
                continue
            cand = TIRCandidate(start + trim, end - trim, tir_len - trim, round(ident, 4), tsd)
            if strict and not cand.tsd:
                continue
            raw.append(cand)
    # Between two copies of one family the left TIR of one and the right TIR of
    # the other also form an inverted repeat; prefer the most compact element.
    return _remove_overlaps(raw, key=lambda c: (c.tir_length * c.identity + 3 * len(c.tsd), -c.length))


def _remove_overlaps(cands, key):
    """Keep the best-scoring candidate among overlapping ones."""
    kept = []
    for c in sorted(cands, key=key, reverse=True):
        if all(c.end <= o.start or c.start >= o.end for o in kept):
            kept.append(c)
    return sorted(kept, key=lambda c: c.start)


@dataclass(frozen=True)
class RepeatMask:
    masked: str  # soft-masked sequence (repeats in lowercase)
    intervals: list[tuple[int, int]]

    @property
    def repeat_fraction(self) -> float:
        covered = sum(e - s for s, e in self.intervals)
        return covered / len(self.masked) if self.masked else 0.0


def kmer_repeat_mask(
    seq: str, k: int = 15, min_copies: int = 3, min_length: int = 50, merge_gap: int = 50
) -> RepeatMask:
    """Soft-mask regions made of k-mers that occur at least ``min_copies`` times.

    Copies on the opposite strand count too (canonical k-mers), as TEs insert in
    either orientation.
    """
    seq = seq.upper()
    counts: Counter[str] = Counter()
    canon = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        c = min(kmer, reverse_complement(kmer))
        canon.append(c)
        counts[c] += 1
    flag = [False] * len(seq)
    for i, c in enumerate(canon):
        if counts[c] >= min_copies and "N" not in c:
            for j in range(i, i + k):
                flag[j] = True
    intervals: list[tuple[int, int]] = []
    i = 0
    while i < len(seq):
        if flag[i]:
            j = i
            while j < len(seq) and flag[j]:
                j += 1
            if intervals and i - intervals[-1][1] <= merge_gap:
                intervals[-1] = (intervals[-1][0], j)
            else:
                intervals.append((i, j))
            i = j
        else:
            i += 1
    intervals = [(a, b) for a, b in intervals if b - a >= min_length]
    chars = list(seq)
    for s, e in intervals:
        for j in range(s, e):
            chars[j] = chars[j].lower()
    return RepeatMask("".join(chars), intervals)


@dataclass(frozen=True)
class TEHit:
    family: str
    te_class: str
    identity: float
    score: int
    strand: str
    coverage: float  # fraction of the library consensus covered by the alignment


def classify(element: str, library: list[Record], k: int = 11) -> TEHit | None:
    """Assign an element to the most similar family in a TE library.

    A quick shared-k-mer filter picks candidates, then Smith-Waterman scores
    them on both strands. The class is read from the library header, in the
    RepeatMasker style ``family#Class/Superfamily``.
    """
    element = element.upper()
    el_kmers = {element[i : i + k] for i in range(len(element) - k + 1)}
    el_kmers |= {reverse_complement(x) for x in el_kmers}
    scored = []
    for rec in library:
        lib_kmers = {rec.seq[i : i + k] for i in range(len(rec.seq) - k + 1)}
        shared = len(lib_kmers & el_kmers) / max(1, len(lib_kmers))
        scored.append((shared, rec))
    scored.sort(key=lambda t: t[0], reverse=True)
    best: TEHit | None = None
    for shared, rec in scored[:3]:
        if shared == 0:
            continue
        for strand, query in (("+", element), ("-", reverse_complement(element))):
            aln = smith_waterman(query, rec.seq)
            if best is None or aln.score > best.score:
                family, _, te_class = rec.id.partition("#")
                best = TEHit(
                    family=family,
                    te_class=te_class or "Unknown",
                    identity=round(aln.identity, 4),
                    score=aln.score,
                    strand=strand,
                    coverage=round((aln.end_b - aln.start_b) / len(rec.seq), 4),
                )
    return best


def poly_a_tail(seq: str, min_length: int = 10, window: int = 40) -> bool:
    """True if the 3' end (either strand) carries an A-rich tail, typical of SINEs/LINEs."""
    seq = seq.upper()
    for s in (seq, reverse_complement(seq)):
        tail = s[-window:]
        run = best = 0
        for c in tail:
            run = run + 1 if c == "A" else 0
            best = max(best, run)
        if best >= min_length:
            return True
    return False
