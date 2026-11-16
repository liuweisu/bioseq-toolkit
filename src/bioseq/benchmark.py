"""Compare predicted element coordinates with a known truth set."""

from __future__ import annotations

from dataclasses import dataclass


def reciprocal_overlap(a: tuple[int, int], b: tuple[int, int]) -> float:
    inter = min(a[1], b[1]) - max(a[0], b[0])
    if inter <= 0:
        return 0.0
    return min(inter / (a[1] - a[0]), inter / (b[1] - b[0]))


@dataclass(frozen=True)
class Score:
    true_positives: int
    predicted: int
    truth: int
    exact: int  # predictions with both boundaries exactly right

    @property
    def precision(self) -> float:
        return self.true_positives / self.predicted if self.predicted else 0.0

    @property
    def recall(self) -> float:
        return self.true_positives / self.truth if self.truth else 0.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if p + r else 0.0


def score_intervals(
    predicted: list[tuple[int, int]], truth: list[tuple[int, int]], min_overlap: float = 0.9
) -> Score:
    """A prediction is correct if it overlaps an unmatched true element reciprocally by ``min_overlap``."""
    unmatched = list(truth)
    tp = 0
    for p in predicted:
        match = next((t for t in unmatched if reciprocal_overlap(p, t) >= min_overlap), None)
        if match is not None:
            unmatched.remove(match)
            tp += 1
    exact = len(set(predicted) & set(truth))
    return Score(tp, len(predicted), len(truth), exact)
