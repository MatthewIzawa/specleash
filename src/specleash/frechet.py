"""Exact discrete Fréchet distance between two polylines.

The recurrence is that of Eiter and Mannila (1994). With d[i][j] the Euclidean
distance between point i of P and point j of Q,

    c[0][0] = d[0][0]
    c[i][j] = max(min(c[i-1][j], c[i-1][j-1], c[i][j-1]), d[i][j])

and the distance is c[n-1][m-1]. It is computed iteratively, so curve length
is limited by memory, never by a recursion limit.

Each cell depends on its left neighbour in the same row, so a row cannot be
computed from the row above in one vectorised step: doing so silently drops a
predecessor. The loop is plain Python over lists for that reason.
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FrechetResult:
    """The distance, the matching that achieves it, and its longest link."""

    distance: float
    path: np.ndarray          # (k, 2) int index pairs, from (0, 0) to (n-1, m-1)
    link_lengths: np.ndarray  # (k,) distance between the matched points
    bottleneck: int           # index into `path` of the first longest link


def _as_curve(x, name):
    a = np.asarray(x, dtype=float)
    if a.ndim != 2 or a.shape[1] != 2 or len(a) == 0:
        raise ValueError(f"{name} must be an (N, 2) array with N >= 1; got shape {a.shape}")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains non-finite values")
    return a


def _table(p, q, keep):
    d = np.hypot(p[:, None, 0] - q[None, :, 0], p[:, None, 1] - q[None, :, 1])
    d_rows = d.tolist()
    n, m = d.shape
    prev = d_rows[0][:]
    for j in range(1, m):
        if prev[j] < prev[j - 1]:
            prev[j] = prev[j - 1]
    rows = [prev] if keep else None
    for i in range(1, n):
        d_row = d_rows[i]
        row = [0.0] * m
        row[0] = prev[0] if prev[0] > d_row[0] else d_row[0]
        left, diag = row[0], prev[0]
        for j in range(1, m):
            up = prev[j]
            best = up if up < diag else diag
            if left < best:
                best = left
            value = d_row[j]
            left = best if best > value else value
            row[j] = left
            diag = up
        if keep:
            rows.append(row)
        prev = row
    return d, (rows if keep else [prev])


def discrete_frechet(p, q) -> float:
    """Discrete Fréchet distance between two (N, 2) polylines."""
    p, q = _as_curve(p, "p"), _as_curve(q, "q")
    _, rows = _table(p, q, keep=False)
    return float(rows[-1][-1])


def frechet_coupling(p, q) -> FrechetResult:
    """The distance and a matching that achieves it.

    The path is recovered by stepping back from (n-1, m-1) to the predecessor
    with the smallest table value, preferring the diagonal on ties. Along such a
    path the longest link equals the distance.
    """
    p, q = _as_curve(p, "p"), _as_curve(q, "q")
    d, rows = _table(p, q, keep=True)
    c = np.array(rows)
    i, j = len(p) - 1, len(q) - 1
    path = [(i, j)]
    while i > 0 or j > 0:
        if i == 0:
            j -= 1
        elif j == 0:
            i -= 1
        else:
            _, i, j = min(((c[i - 1, j - 1], i - 1, j - 1),
                           (c[i - 1, j], i - 1, j),
                           (c[i, j - 1], i, j - 1)), key=lambda t: t[0])
        path.append((i, j))
    path = np.array(path[::-1], dtype=int)
    links = d[path[:, 0], path[:, 1]]
    return FrechetResult(float(c[-1, -1]), path, links, int(np.argmax(links)))
