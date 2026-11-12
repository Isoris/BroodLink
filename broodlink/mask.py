"""Callable mask: regions to keep (BED or soft-masked FASTA), regions to drop (BED).

The mask is the denominator of every per-base quantity. Rule for point markers: the position is
callable. Rule for structural variants (new in BroodLink): at least a set fraction of the span and
both breakpoints are callable.
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd

from .io import _open, InputError, log


def _merge(s, e):
    o = np.argsort(s, kind="stable"); s, e = s[o], e[o]
    if len(s) == 0:
        return s, e
    cm = np.maximum.accumulate(e)
    new = np.r_[True, s[1:] > cm[:-1]]
    idx = np.flatnonzero(new)
    return s[idx], np.maximum.reduceat(e, idx)


def read_bed(path):
    d = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2], names=["c", "s", "e"], dtype={"c": str})
    return {c: _merge(g.s.to_numpy(np.int64), g.e.to_numpy(np.int64)) for c, g in d.groupby("c", sort=False)}


def read_softmasked_fasta(path):
    """Uppercase A/C/G/T = callable (lowercase repeats and N excluded), as in the original mask."""
    out, name, chunks = {}, None, []
    lut = np.zeros(256, bool)
    for ch in b"ACGT":
        lut[ch] = True

    def flush():
        if name is None:
            return
        seq = np.frombuffer(b"".join(chunks), dtype=np.uint8)
        ok = lut[seq]
        d = np.diff(np.r_[0, ok.astype(np.int8), 0])
        out[name] = (np.flatnonzero(d == 1).astype(np.int64), np.flatnonzero(d == -1).astype(np.int64))
    with (open(path, "rb") if not path.endswith(".gz") else __import__("gzip").open(path, "rb")) as fh:
        for line in fh:
            if line.startswith(b">"):
                flush(); name = line[1:].split()[0].decode(); chunks = []
            else:
                chunks.append(line.strip())
    flush()
    return out


def is_fasta(path):
    with _open(path) as fh:
        return fh.readline().startswith(">")


def intersect(a, b):
    out = {}
    for c in a:
        if c not in b:
            continue
        s1, e1 = a[c]; s2, e2 = b[c]
        i = j = 0; S, E = [], []
        while i < len(s1) and j < len(s2):
            lo, hi = max(s1[i], s2[j]), min(e1[i], e2[j])
            if lo < hi:
                S.append(lo); E.append(hi)
            if e1[i] < e2[j]:
                i += 1
            else:
                j += 1
        out[c] = (np.array(S, np.int64), np.array(E, np.int64))
    return out


def subtract(a, b):
    out = {}
    for c, (s1, e1) in a.items():
        if c not in b:
            out[c] = (s1, e1); continue
        s2, e2 = b[c]; S, E = [], []
        j = 0
        for lo, hi in zip(s1, e1):
            cur = lo
            while j < len(s2) and e2[j] <= cur:
                j += 1
            k = j
            while k < len(s2) and s2[k] < hi:
                if s2[k] > cur:
                    S.append(cur); E.append(s2[k])
                cur = max(cur, e2[k]); k += 1
            if cur < hi:
                S.append(cur); E.append(hi)
        out[c] = (np.array(S, np.int64), np.array(E, np.int64))
    return out


class Mask:
    def __init__(self, keep=None, genome=None, description="unmasked"):
        self.keep = keep                     # dict chrom -> (starts, ends), 0-based half-open; None = unmasked
        self.genome = genome or {}
        self.description = description
        self._cum = {}
        if keep is not None:
            for c, (s, e) in keep.items():
                self._cum[c] = np.r_[0, np.cumsum(e - s)]

    @property
    def active(self):
        return self.keep is not None

    def span(self):
        return int(sum(int((e - s).sum()) for s, e in self.keep.values())) if self.active else None

    def assembly(self):
        return int(sum(self.genome.values())) if self.genome else None

    def label(self):
        """Text placed in the header of every per-base table."""
        if not self.active:
            return "callable span: unmasked" + (f" (assembly {self.assembly():,} bp)" if self.genome else "")
        sp, asm = self.span(), self.assembly()
        return f"callable span: {sp:,} bp" + (f" ({100 * sp / asm:.2f}% of the {asm:,}-bp assembly)" if asm else " (assembly length not given: pass --genome)")

    def contains(self, chrom, pos1):
        """1-based positions -> bool."""
        if not self.active:
            return np.ones(len(pos1), bool)
        out = np.zeros(len(pos1), bool)
        chrom = np.asarray(chrom)
        for c in pd.unique(chrom):
            m = chrom == c
            if c not in self.keep:
                continue
            s, e = self.keep[c]; p = np.asarray(pos1)[m] - 1
            i = np.searchsorted(s, p, side="right") - 1
            out[m] = (i >= 0) & (p < e[np.clip(i, 0, None)])
        return out

    def covered(self, chrom, start1, end1):
        """Callable bases inside 1-based inclusive spans."""
        if not self.active:
            return (np.asarray(end1) - np.asarray(start1) + 1).astype(np.int64)
        out = np.zeros(len(start1), np.int64)
        chrom = np.asarray(chrom)
        for c in pd.unique(chrom):
            m = chrom == c
            if c not in self.keep:
                continue
            s, e = self.keep[c]; cum = self._cum[c]

            def upto(x):                     # callable bases in [0, x)
                i = np.searchsorted(s, x, side="right") - 1
                ii = np.clip(i, 0, None)
                return np.where(i >= 0, cum[ii] + np.minimum(x, e[ii]) - s[ii], 0)
            out[m] = upto(np.asarray(end1)[m]) - upto(np.asarray(start1)[m] - 1)
        return out


def build_mask(keep_paths, exclude_paths, genome):
    if not keep_paths and not exclude_paths:
        return Mask(None, genome)
    keep, desc = None, []
    for p in keep_paths or []:
        if not os.path.exists(p):
            raise InputError(f"mask file not found: {p}")
        k = read_softmasked_fasta(p) if is_fasta(p) else read_bed(p)
        desc.append(("soft-masked FASTA, uppercase A/C/G/T" if is_fasta(p) else "BED") + f" {os.path.basename(p)}")
        keep = k if keep is None else intersect(keep, k)
    if keep is None:
        if not genome:
            raise InputError("--mask-exclude without --mask needs --genome (FAI) to define the callable span")
        keep = {c: (np.array([0], np.int64), np.array([L], np.int64)) for c, L in genome.items()}
        desc.append("whole assembly")
    for p in exclude_paths or []:
        keep = subtract(keep, read_bed(p)); desc.append(f"minus BED {os.path.basename(p)}")
    m = Mask(keep, genome, "; ".join(desc))
    log(f"mask: {m.description}; {m.label()}")
    return m


def apply_mask(ms, mask, cfg):
    """Restrict a MarkerSet to the mask; returns (masked set, report row)."""
    row = {"input": ms.name, "marker_class": ms.klass, "markers_total": ms.n, "markers_inside_mask": ms.n,
           "markers_outside_mask": 0, "rule": "unmasked", "warning": ""}
    if not mask.active or ms.klass == "HAPLOTYPE":
        return ms, row
    if ms.klass in ("SNP", "INDEL", "SNP_INDEL"):
        keep = mask.contains(ms.chrom, ms.pos); row["rule"] = "position callable"
    else:
        fr = cfg["mask"]["sv_min_callable_fraction"]
        span = np.maximum(ms.end - ms.pos + 1, 1)
        cov = mask.covered(ms.chrom, ms.pos, np.maximum(ms.end, ms.pos))
        keep = (cov / span >= fr) & mask.contains(ms.chrom, ms.pos) & mask.contains(ms.chrom, np.maximum(ms.end, ms.pos))
        row["rule"] = f">= {fr} of the span and both breakpoints callable"
    row["markers_inside_mask"] = int(keep.sum()); row["markers_outside_mask"] = int((~keep).sum())
    if ms.klass not in ("SNP", "INDEL", "SNP_INDEL") and ms.n and (~keep).mean() > cfg["mask"]["sv_class_loss_warning"]:
        row["warning"] = (f"the mask removes {100 * (~keep).mean():.0f}% of the {ms.klass} calls: structural variants in repeats "
                          "(e.g. transposable-element insertions) are being dropped; check that this is intended")
        log("WARNING: " + row["warning"])
    return ms.subset_markers(keep, "callable mask"), row
