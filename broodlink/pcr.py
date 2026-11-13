"""Step 10: PCR markers of the large structural haplotypes.

For every structural-haplotype interval the step screens the deletions inside it, keeps those
whose allele frequency separates the two homokaryotype classes and whose size can be read on a
gel, designs one flanking primer pair on the sequence across the deletion junction, checks the
pair, and writes the primer tables. One pair and one reaction per individual then read the
arrangement as a short band (allele with the deletion), a long band (allele without) or both.

Everything is an in-silico design: no primer was tested in a laboratory.

Where the rules come from is written in docs/EXTRACTED_METHOD.md (section "PCR markers").
In short: the published primers of the original study were designed with Primer3 (the design
script was not found; the settings below return the published pairs and are recorded as
"recovered by reproduction"); the screening rules are those of the published screening table;
the built-in design routine is the one of the later deletion toolkit (`pick_flanking_primers`),
kept as the fallback when the optional package primer3-py is not installed.
"""
from __future__ import annotations
import math
import os
import re

import numpy as np
import pandas as pd

from .io import InputError, detail, log

try:                                                            # optional dependency (GPL-2): pip install broodlink[primers]
    import primer3 as _p3
    HAVE_PRIMER3 = True
    PRIMER3_VERSION = getattr(_p3, "__version__", "unknown")
except Exception:                                               # not installed, or a broken build
    _p3 = None
    HAVE_PRIMER3 = False
    PRIMER3_VERSION = ""

# ----------------------------------------------------------------------------- sequence helpers
_COMP = str.maketrans("ACGTacgtNn", "TGCATGCANN")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def gc(s):
    s = s.upper()
    return 100.0 * (s.count("G") + s.count("C")) / max(1, len(s))


def max_run(s):
    s = s.upper(); best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best


# ----------------------------------------------------------------------------- melting temperature
# SantaLucia (1998) unified nearest-neighbour parameters: dH kcal/mol, dS cal/(mol K)
_NN_H = {"AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2, "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,
         "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2, "CG": -10.6, "GC": -9.8, "GG": -8.0, "CC": -8.0}
_NN_S = {"AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4,
         "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2, "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9}
TM_DEFAULT = {"mv_mM": 50.0, "dv_mM": 1.5, "dntp_mM": 0.6, "dna_nM": 50.0}


def tm_formula(seq, mv_mM=50.0, dv_mM=1.5, dntp_mM=0.6, dna_nM=50.0):
    """Nearest-neighbour Tm in pure Python with the model of Primer3's `calc_tm` defaults:
    SantaLucia (1998) parameters and salt correction, divalent cations converted to a monovalent
    equivalent (120 * sqrt([Mg] - [dNTP])), strand concentration dna_nM. Equal to
    `primer3.calc_tm` within rounding error (tested). None for sequences with other letters."""
    s = seq.upper()
    if len(s) < 2 or any(b not in "ACGT" for b in s):
        return None
    sym = s == revcomp(s)
    dh, ds = 0.0, (-1.4 if sym else 0.0)
    for e in (s[0], s[-1]):
        if e in "AT":
            dh += 2.3; ds += 4.1
        else:
            dh += 0.1; ds += -2.8
    for i in range(len(s) - 1):
        dh += _NN_H[s[i:i + 2]]; ds += _NN_S[s[i:i + 2]]
    k = mv_mM + (120.0 * math.sqrt(dv_mM - dntp_mM) if dv_mM > dntp_mM else 0.0)
    ds += 0.368 * (len(s) - 1) * math.log(k / 1000.0)
    return dh * 1000.0 / (ds + 1.987 * math.log(dna_nM / (1e9 if sym else 4e9))) - 273.15


def tm_toolkit_fallback(seq, primer_nM=250.0, na_mM=50.0):
    """The fallback Tm of the original deletion toolkit (`deltools.nn_tm`), kept for the record
    only: about 10 degrees C below `calc_tm` (250 nM, 16.6 log10 salt term, no divalent cations).
    BroodLink does not design with it."""
    s = seq.upper()
    if len(s) < 8 or any(b not in "ACGT" for b in s):
        return None
    dh = ds = 0.0
    for e in (s[0], s[-1]):
        if e in "GC":
            dh += 0.1; ds += -2.8
        else:
            dh += 2.3; ds += 4.1
    for i in range(len(s) - 1):
        dh += _NN_H[s[i:i + 2]]; ds += _NN_S[s[i:i + 2]]
    t = (dh * 1000.0) / (ds + 1.987 * math.log(primer_nM * 1e-9 / 4.0)) - 273.15
    return float(t + 16.6 * math.log10(max(na_mM, 1e-3) / 1000.0))


def tm(seq, cond=None):
    """Primer Tm as in the original tables: `primer3.calc_tm` when primer3-py is installed, else
    the built-in formula with the same model and conditions."""
    c = {**TM_DEFAULT, **(cond or {})}
    if HAVE_PRIMER3:
        try:
            return float(_p3.calc_tm(seq, mv_conc=c["mv_mM"], dv_conc=c["dv_mM"], dntp_conc=c["dntp_mM"], dna_conc=c["dna_nM"]))
        except Exception:
            pass
    return tm_formula(seq, c["mv_mM"], c["dv_mM"], c["dntp_mM"], c["dna_nM"])


def tm_conditions(cond=None):
    c = {**TM_DEFAULT, **(cond or {})}
    how = (f"primer3-py {PRIMER3_VERSION} calc_tm" if HAVE_PRIMER3 else "built-in formula equal to Primer3 calc_tm (primer3-py not installed)")
    return (f"Tm: nearest neighbour, SantaLucia (1998) parameters and salt correction, {c['mv_mM']:g} mM monovalent, {c['dv_mM']:g} mM divalent, "
            f"{c['dntp_mM']:g} mM dNTP, {c['dna_nM']:g} nM primer ({how})")


def thermo(fwd, rev, cond=None):
    """Hairpin and dimer Tm / dG of a pair (new in BroodLink; needs primer3-py). {} when absent."""
    if not HAVE_PRIMER3:
        return {}
    c = {**TM_DEFAULT, **(cond or {})}
    kw = dict(mv_conc=c["mv_mM"], dv_conc=c["dv_mM"], dntp_conc=c["dntp_mM"], dna_conc=c["dna_nM"])
    out = {}
    try:
        for tag, r in (("hairpin_forward", _p3.calc_hairpin(fwd, **kw)), ("hairpin_reverse", _p3.calc_hairpin(rev, **kw)),
                       ("self_dimer_forward", _p3.calc_homodimer(fwd, **kw)), ("self_dimer_reverse", _p3.calc_homodimer(rev, **kw)),
                       ("cross_dimer", _p3.calc_heterodimer(fwd, rev, **kw))):
            found = bool(getattr(r, "structure_found", False))
            out[tag + "_tm"] = round(float(r.tm), 1) if found else np.nan
            out[tag + "_dg_kcal"] = round(float(r.dg) / 1000.0, 2) if found else np.nan
    except Exception:
        return {}
    return out


THERMO_COLS = [f"{t}_{k}" for t in ("hairpin_forward", "hairpin_reverse", "self_dimer_forward", "self_dimer_reverse", "cross_dimer") for k in ("tm", "dg_kcal")]


# ----------------------------------------------------------------------------- design engines
def gel_resolvable_design(band_deleted, band_undeleted):
    """Toolkit rule on a designed pair (deltools.pick_flanking_primers, key 'gel_resolvable')."""
    return bool(abs(band_undeleted - band_deleted) >= 60 and 60 <= band_deleted <= 1200 and band_undeleted <= 6000)


def _result(fwd, rev, f0, r_last, bp, svlen, tm_f, tm_r, engine, extra=None):
    band_del = int(r_last - f0 + 1)
    out = {"status": "ok", "fwd": fwd, "rev": rev, "len_F": len(fwd), "len_R": len(rev), "Tm_F": round(float(tm_f), 1), "Tm_R": round(float(tm_r), 1),
           "GC_F": round(gc(fwd), 1), "GC_R": round(gc(rev), 1), "fpos": int(f0), "rpos": int(r_last), "consbp": int(bp),
           "band_deleted": band_del, "band_undeleted": int(band_del + int(svlen)), "engine": engine}
    out.update(extra or {})
    return out


def design_builtin(template, bp, svlen, c=None, cond=None):
    """The routine of the deletion toolkit (`deltools.pick_flanking_primers`), unchanged: on either
    side of the breakpoint every window of len_min..len_max nt whose 3' end lies gap_min..gap_max
    nt from the breakpoint is scored |Tm - target| + 0.6 when the 3' base is A or T; windows
    outside the Tm or GC range or with a longer homopolymer are dropped; the best window per
    side is kept. Positions are 0-based on the template; rpos is the last base of the reverse site."""
    c = {"tm_lo": 54.5, "tm_hi": 56.6, "tm_target": 55.5, "len_min": 18, "len_max": 26, "gc_lo": 30.0, "gc_hi": 60.0,
         "gap_min": 40, "gap_max": 160, "max_homopolymer": 4, **(c or {})}
    if not isinstance(template, str) or len(template) < 60 or bp is None:
        return {"status": "no_consensus"}
    bp = int(bp); template = template.upper()
    if bp < c["len_min"] + c["gap_min"] or bp > len(template) - (c["len_min"] + c["gap_min"]):
        return {"status": "breakpoint_too_close_to_end"}

    def score(seq):
        t = tm(seq, cond)
        if t is None or not (c["tm_lo"] <= t <= c["tm_hi"]):
            return None
        g = gc(seq)
        if not (c["gc_lo"] <= g <= c["gc_hi"]) or max_run(seq) > c["max_homopolymer"]:
            return None
        return abs(t - c["tm_target"]) + (0.0 if seq[-1] in "GC" else 0.6), t

    best_f = best_r = None
    for gap in range(c["gap_min"], c["gap_max"] + 1):
        end = bp - gap
        for n in range(c["len_min"], c["len_max"] + 1):
            start = end - n
            if start < 0:
                continue
            seq = template[start:end]; sc = score(seq)
            if sc and (best_f is None or sc[0] < best_f[0]):
                best_f = (sc[0], seq, start, sc[1])
    for gap in range(c["gap_min"], c["gap_max"] + 1):
        start = bp + gap
        for n in range(c["len_min"], c["len_max"] + 1):
            end = start + n
            if end > len(template):
                continue
            seq = revcomp(template[start:end]); sc = score(seq)
            if sc and (best_r is None or sc[0] < best_r[0]):
                best_r = (sc[0], seq, end - 1, sc[1])
    if best_f is None or best_r is None:
        return {"status": "no_primer_in_window"}
    return _result(best_f[1], best_r[1], best_f[2], best_r[2], bp, svlen, best_f[3], best_r[3], "built-in routine")


def primer3_args(c, cond=None):
    """Primer3 global arguments from the config block pcr.primer3 and the Tm conditions."""
    t = {**TM_DEFAULT, **(cond or {})}
    g = {"PRIMER_NUM_RETURN": 1, "PRIMER_OPT_TM": float(c["opt_tm"]), "PRIMER_MIN_TM": float(c["min_tm"]), "PRIMER_MAX_TM": float(c["max_tm"]),
         "PRIMER_OPT_SIZE": int(c["opt_size"]), "PRIMER_MIN_SIZE": int(c["min_size"]), "PRIMER_MAX_SIZE": int(c["max_size"]),
         "PRIMER_MIN_GC": float(c["min_gc"]), "PRIMER_MAX_GC": float(c["max_gc"]), "PRIMER_MAX_POLY_X": int(c["max_poly_x"]),
         "PRIMER_PRODUCT_SIZE_RANGE": [[int(c["product_min"]), int(c["product_max"])]],
         "PRIMER_SALT_MONOVALENT": t["mv_mM"], "PRIMER_SALT_DIVALENT": t["dv_mM"], "PRIMER_DNTP_CONC": t["dntp_mM"], "PRIMER_DNA_CONC": t["dna_nM"]}
    g.update(c.get("extra") or {})
    return g


def design_primer3(template, bp, svlen, c, cond=None):
    """Primer3 (primer3-py `design_primers`): the best pair that contains the breakpoint (target =
    the two bases on either side of it). Positions follow Primer3: fpos = first base of the forward
    primer, rpos = last base of the reverse site, both 0-based on the template."""
    if not HAVE_PRIMER3:
        return {"status": "primer3_not_installed"}
    if not isinstance(template, str) or len(template) < 60 or bp is None:
        return {"status": "no_consensus"}
    bp = int(bp); template = template.upper()
    if not (0 < bp < len(template)):
        return {"status": "breakpoint_outside_sequence"}
    tl, tn = int(c.get("target_left", 1)), int(c.get("target_length", 2))
    try:
        r = _p3.design_primers({"SEQUENCE_TEMPLATE": template, "SEQUENCE_TARGET": [max(bp - tl, 0), tn]}, primer3_args(c, cond))
    except Exception as exc:
        return {"status": "primer3_error", "detail": str(exc)[:200]}
    if not r.get("PRIMER_PAIR_NUM_RETURNED"):
        return {"status": "no_primer_pair", "detail": "; ".join(str(r.get(k, "")) for k in ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN"))[:300]}
    f0, _fl = r["PRIMER_LEFT_0"]; r_last, _rl = r["PRIMER_RIGHT_0"]
    return _result(r["PRIMER_LEFT_0_SEQUENCE"], r["PRIMER_RIGHT_0_SEQUENCE"], f0, r_last, bp, svlen, r["PRIMER_LEFT_0_TM"], r["PRIMER_RIGHT_0_TM"],
                   f"Primer3 (primer3-py {PRIMER3_VERSION})", {"primer3_pair_penalty": round(float(r.get("PRIMER_PAIR_0_PENALTY", np.nan)), 3)})


def resolve_engine(cfg, cli=None):
    """(engine actually used, message). 'primer3' falls back to the built-in routine, with a
    message, when primer3-py cannot be imported; nothing fails."""
    want = (cli or cfg["pcr"].get("engine") or "primer3").lower()
    if want not in ("primer3", "builtin"):
        raise InputError(f"pcr.engine must be 'primer3' or 'builtin' (got '{want}')")
    if want == "primer3" and not HAVE_PRIMER3:
        return "builtin", ("primer design: Primer3 was requested (the original procedure) but the optional package primer3-py is not installed "
                           "(pip install broodlink[primers]); the built-in routine was used instead. It does not return the same primers as Primer3")
    return want, ""


def design(template, bp, svlen, engine, cfg):
    c = cfg["pcr"]
    res = design_primer3(template, bp, svlen, c["primer3"], c["tm"]) if engine == "primer3" else design_builtin(template, bp, svlen, c["builtin"], c["tm"])
    if res["status"] == "ok":
        res["gel_resolvable_design"] = gel_resolvable_design(res["band_deleted"], res["band_undeleted"])
    return res


# ----------------------------------------------------------------------------- reference FASTA
class Fasta:
    """Small indexed FASTA reader in pure Python. Uses <fasta>.fai when present, otherwise scans the
    file once and keeps the index in memory (nothing is written next to the input). Plain text
    FASTA only (not gzip)."""

    def __init__(self, path):
        if str(path).endswith((".gz", ".bgz")):
            raise InputError("--reference must be an uncompressed FASTA (random access is needed); decompress it first")
        if not os.path.exists(path):
            raise InputError(f"input file not found: {path}")
        self.path = path; self.index = {}; self.indexed_from = "scan of the file"
        fai = str(path) + ".fai"
        if os.path.exists(fai) and os.path.getmtime(fai) >= os.path.getmtime(path) - 1:
            with open(fai, encoding="utf-8") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 5:
                        self.index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
            self.indexed_from = ".fai"
        if not self.index:
            self._scan()
        self._fh = open(path, "rb")

    def _scan(self):
        name = None; off = 0; n = 0; lb = lw = 0; start = 0; first = True
        with open(self.path, "rb") as fh:
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        self.index[name] = (n, start, lb, lw)
                    name = raw[1:].split()[0].decode(); n = 0; lb = lw = 0; start = off + len(raw); first = True
                elif name is not None:
                    b = len(raw.rstrip(b"\r\n"))
                    if first and b:
                        lb, lw, first = b, len(raw), False
                    n += b
                off += len(raw)
        if name is not None:
            self.index[name] = (n, start, lb, lw)

    def __contains__(self, chrom):
        return chrom in self.index

    def length(self, chrom):
        return self.index[chrom][0]

    def fetch(self, chrom, start, end):
        """Sequence of chrom[start:end] (0-based, end exclusive), upper case; clipped to the sequence."""
        n, off, lb, lw = self.index[chrom]
        start = max(0, int(start)); end = min(n, int(end))
        if end <= start or lb <= 0:
            return ""
        b0 = off + (start // lb) * lw + start % lb
        b1 = off + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        self._fh.seek(b0)
        return self._fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"").decode("ascii", "replace").upper()

    def chromosome(self, chrom):
        return self.fetch(chrom, 0, self.length(chrom))

    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass


def junction_from_reference(fa, chrom, pos, end, flank):
    """Sequence across the junction of a deletion cut from the reference: `flank` bases up to and
    including POS (the last base before the deletion, VCF convention) joined to `flank` bases
    after END (the last deleted base). Returns (sequence, breakpoint index) or (None, None)."""
    if fa is None or chrom not in fa:
        return None, None
    left = fa.fetch(chrom, int(pos) - flank, int(pos)); right = fa.fetch(chrom, int(end), int(end) + flank)
    if len(left) < 30 or len(right) < 30:
        return None, None
    return left + right, len(left)


# ----------------------------------------------------------------------------- in-silico PCR
def _hits(arr, kcache, query, three_prime_last, max_mm):
    """Start positions (0-based) and mismatch counts of sites matching `query` (+ strand letters as
    codes 0-3) with at most max_mm mismatches; the 3' terminal base must match. Pigeonhole seeds:
    the query is cut in max_mm + 1 pieces, one of which must match exactly."""
    n = len(arr); L = len(query); k = L // (max_mm + 1)
    if k < 4 or n < L:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    cand = []
    for i in range(max_mm + 1):
        off = i * k; piece = query[off:off + k] if i < max_mm else query[off:]
        kk = len(piece)
        if kk not in kcache:
            c = np.zeros(n - kk + 1, np.int64); bad = np.zeros(n - kk + 1, bool)
            for j in range(kk):
                seg = arr[j:n - kk + 1 + j]; c = c * 4 + seg; bad |= seg > 3
            c[bad] = -1; kcache[kk] = c
        code = 0
        for v in piece:
            code = code * 4 + int(v)
        p = np.flatnonzero(kcache[kk] == code) - off
        cand.append(p[(p >= 0) & (p <= n - L)])
    cp = np.unique(np.concatenate(cand)) if cand else np.zeros(0, np.int64)
    if not len(cp):
        return cp, cp
    W = np.stack([arr[cp + j] for j in range(L)], axis=1); mm = (W != query[None, :]).sum(1)
    ok = (mm <= max_mm) & ((W[:, -1] == query[-1]) if three_prime_last else (W[:, 0] == query[0]))
    return cp[ok], mm[ok]


def in_silico_pcr(fa, pairs, c):
    """Search every primer pair in the reference (ported from the in-silico check of the original
    study): sites with at most `max_mismatches` mismatches per primer and a matching 3' base;
    products = forward site on the + strand and reverse site on the - strand of one sequence,
    `product_min`..`product_max` bp. A product is "at the locus" when it overlaps the deletion
    +- `locus_window_bp`. pairs: list of dicts with pair, chrom, start, end, fwd, rev.
    Returns one dict per pair."""
    MAP = np.full(256, 4, np.uint8)
    for i, ch in enumerate("ACGT"):
        MAP[ord(ch)] = i
    q = []
    for p in pairs:
        f = p["fwd"].upper(); r = revcomp(p["rev"].upper())
        if any(ch not in "ACGT" for ch in f + r):
            continue
        q.append((p["pair"], MAP[np.frombuffer(f.encode(), np.uint8)].astype(np.int64), MAP[np.frombuffer(r.encode(), np.uint8)].astype(np.int64)))
    prods = {p["pair"]: [] for p in pairs}; nsite = {p["pair"]: [0, 0] for p in pairs}
    mm_max, lo, hi = int(c["max_mismatches"]), int(c["product_min"]), int(c["product_max"])
    for chrom in fa.index:
        seq = fa.chromosome(chrom)
        if not seq:
            continue
        arr = MAP[np.frombuffer(seq.encode("ascii", "replace"), np.uint8)].astype(np.int64); del seq
        kc = {}
        for name, fq, rq in q:
            fp, fm = _hits(arr, kc, fq, True, mm_max); rp, rm = _hits(arr, kc, rq, False, mm_max)
            nsite[name][0] += len(fp); nsite[name][1] += len(rp)
            if not len(fp) or not len(rp):
                continue
            rend = rp + len(rq)                                 # exclusive end of the reverse site
            for a, m1 in zip(fp, fm):
                j0 = np.searchsorted(rp, a + 1); j1 = np.searchsorted(rend, a + hi, side="right")
                for j in range(j0, j1):
                    size = int(rend[j] - a)
                    if lo <= size <= hi:
                        prods[name].append((chrom, int(a) + 1, int(rend[j]), size, int(m1 + rm[j])))
    out = []
    for p in pairs:
        pr = prods[p["pair"]]; w = int(c["locus_window_bp"])
        on = [x for x in pr if x[0] == p["chrom"] and x[1] <= p["end"] + w and x[2] >= p["start"] - w]
        off = [x for x in pr if x not in on]
        best = min(on, key=lambda x: (x[4], abs(x[3] - p.get("expected", x[3])))) if on else None
        out.append({"pair": p["pair"], "sites_forward": nsite[p["pair"]][0], "sites_reverse": nsite[p["pair"]][1], "products_at_locus": len(on),
                    "product_at_locus_bp": best[3] if best else np.nan, "product_at_locus_mismatches": best[4] if best else np.nan,
                    "products_off_target": len(off), "off_target_products": ";".join(f"{x[0]}:{x[1]}-{x[2]}({x[3]}bp,{x[4]}mm)" for x in off[:6])})
    return out


# ----------------------------------------------------------------------------- screening
SCREEN_COLS = ["interval", "chrom", "start", "end", "variant_id", "length_bp", "arrangement_with_deletion", "deletion_frequency_class_0",
               "deletion_frequency_heterokaryotypes", "deletion_frequency_class_2", "frequency_difference", "carrier_frequency_class_0",
               "carrier_frequency_class_2", "n_called_class_0", "n_called_heterokaryotypes", "n_called_class_2", "screen_band_deleted_bp",
               "screen_band_undeleted_bp", "diagnostic", "gel_resolvable", "step06_tag", "sequence_source", "rank_in_interval", "design_status", "selected"]


def parse_interval(name):
    m = re.match(r"(?:scan_)?(.+):(\d+)-(\d+)$", str(name))
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


def screen(intervals, geno, samples, dels, cfg, tags=None):
    """One row per interval and deletion lying inside it (POS and END inside, as in the original
    screen), with the deletion allele frequency in the two homokaryotype classes and in the
    heterokaryotypes, and the two flags of the published screening table:
      diagnostic      frequency difference between the homokaryotype classes >= min_frequency_difference
                      and heterokaryotype frequency inside het_frequency_min..het_frequency_max;
      gel_resolvable  deletion length inside min_deletion_bp..max_deletion_bp (with the generic
                      screening pair of 2 x screen_flank_bp the long band stays <= 2 kb).
    intervals: DataFrame interval, chrom, start, end. geno: {interval: int8 array over samples}."""
    s = cfg["pcr"]["selection"]; mc = int(s["min_called_per_class"]); rows = []
    mm = dels.subset_samples(samples); six = {x: i for i, x in enumerate(samples)}
    col = np.array([six[x] for x in mm.samples], int)
    tagset = set()
    if tags is not None and len(tags):
        tagset = set(zip(tags.interval.astype(str), tags.variant_id.astype(str), tags.start.astype(int)))
    flank = int(s["screen_flank_bp"])
    for r in intervals.itertuples():
        if r.interval not in geno:
            continue
        g = np.asarray(geno[r.interval])[col]
        m = np.flatnonzero((mm.chrom == r.chrom) & (mm.pos >= r.start) & (mm.end <= r.end))
        if not len(m):
            continue
        G = mm.G[m]
        st = {}
        for k in (0, 1, 2):
            sub = G[:, g == k]; n = (sub >= 0).sum(1)
            st[k] = (n, np.where(sub > 0, sub, 0).sum(1) / (2.0 * np.maximum(n, 1)), (sub > 0).sum(1) / np.maximum(n, 1))
        ok = (st[0][0] >= mc) & (st[2][0] >= mc)
        for i in np.flatnonzero(ok):
            f0, f2 = float(st[0][1][i]), float(st[2][1][i]); nh = int(st[1][0][i]); fh = float(st[1][1][i]) if nh >= mc else np.nan
            d = abs(f2 - f0); ln = int(max(mm.end[m[i]] - mm.pos[m[i]], 0))
            diag = bool(d >= s["min_frequency_difference"] and (np.isnan(fh) or s["het_frequency_min"] <= fh <= s["het_frequency_max"]))
            rows.append({"interval": r.interval, "chrom": r.chrom, "start": int(mm.pos[m[i]]), "end": int(mm.end[m[i]]), "variant_id": str(mm.ids[m[i]]), "length_bp": ln,
                         "arrangement_with_deletion": "class 2" if f2 > f0 else "class 0", "deletion_frequency_class_0": round(f0, 3),
                         "deletion_frequency_heterokaryotypes": round(fh, 3) if not np.isnan(fh) else np.nan, "deletion_frequency_class_2": round(f2, 3),
                         "frequency_difference": round(d, 3), "_d": d, "carrier_frequency_class_0": round(float(st[0][2][i]), 3), "carrier_frequency_class_2": round(float(st[2][2][i]), 3),
                         "n_called_class_0": int(st[0][0][i]), "n_called_heterokaryotypes": nh, "n_called_class_2": int(st[2][0][i]),
                         "screen_band_deleted_bp": 2 * flank, "screen_band_undeleted_bp": 2 * flank + ln, "diagnostic": diag,
                         "gel_resolvable": bool(s["min_deletion_bp"] <= ln <= s["max_deletion_bp"]),
                         "step06_tag": (str(r.interval), str(mm.ids[m[i]]), int(mm.pos[m[i]])) in tagset,
                         "sequence_source": "", "rank_in_interval": 0, "design_status": "", "selected": False})
    S = pd.DataFrame(rows, columns=SCREEN_COLS + ["_d"])
    if len(S):
        # order of the published screening table: largest frequency difference first; ties by position
        S = S.sort_values(["interval", "_d", "start"], ascending=[True, False, True], kind="stable").reset_index(drop=True)
        S["rank_in_interval"] = S.groupby("interval", sort=False).cumcount() + 1
    return S


# ----------------------------------------------------------------------------- the step
PRIMER_COLS = ["marker", "interval", "variant_id", "chrom", "start", "end", "length_bp", "arrangement_with_deletion", "deletion_frequency_class_0",
               "deletion_frequency_heterokaryotypes", "deletion_frequency_class_2", "frequency_difference", "carrier_frequency_class_0", "carrier_frequency_class_2",
               "forward_primer", "reverse_primer", "length_forward", "length_reverse", "tm_forward", "tm_reverse", "gc_forward_percent", "gc_reverse_percent",
               "product_with_deletion_bp", "product_without_deletion_bp", "product_difference_bp", "design_source", "engine", "design_sequence",
               "breakpoint_index", "forward_start_index", "reverse_end_index",
               "check_primers_exact_in_design_sequence", "check_primers_flank_breakpoint", "check_products_recomputed", "check_frequency_difference",
               "check_gel_resolvable_design", "tm_difference_pair", "min_gc_percent", "forward_3prime_gc", "reverse_3prime_gc", "longest_homopolymer",
               "in_silico_products_at_locus", "in_silico_product_at_locus_bp", "in_silico_product_minus_design_bp", "in_silico_products_off_target",
               "in_silico_off_target_products", "check_in_silico_product_within_5bp"] + THERMO_COLS
ORDER_COLS = ["name", "sequence", "marker", "interval", "variant_id", "length", "tm"]


def read_picks(path):
    """Optional table of deletions to design on instead of BroodLink's own choice: TSV with the
    columns interval and variant_id (other columns are ignored)."""
    d = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    low = {c.lower(): c for c in d.columns}
    ci = next((low[k] for k in ("interval", "inversion", "inv") if k in low), None)
    cv = next((low[k] for k in ("variant_id", "deletion", "sv_id", "id", "del_id") if k in low), None)
    if ci is None or cv is None:
        raise InputError(f"{os.path.basename(path)}: the table of marker picks needs the columns 'interval' and 'variant_id'")
    out = {}
    for a, b in zip(d[ci], d[cv]):
        out.setdefault(str(a), []).append(str(b))
    return out


def _checks(row, res, template, s):
    f, r = res["fwd"], res["rev"]; t = template.upper()
    fs = t.find(f); rs = t.find(revcomp(r))
    row["check_primers_exact_in_design_sequence"] = bool(fs >= 0 and rs >= 0)
    row["check_primers_flank_breakpoint"] = bool(fs >= 0 and rs >= 0 and fs + len(f) <= res["consbp"] <= rs)
    amp = (rs + len(r) - fs) if (fs >= 0 and rs >= 0) else -1
    row["check_products_recomputed"] = bool(amp == res["band_deleted"] and res["band_undeleted"] == res["band_deleted"] + row["length_bp"])
    row["check_frequency_difference"] = bool(row["frequency_difference"] >= s["min_frequency_difference"])
    row["check_gel_resolvable_design"] = bool(res.get("gel_resolvable_design"))
    row["tm_difference_pair"] = round(abs(res["Tm_F"] - res["Tm_R"]), 1)
    row["min_gc_percent"] = min(res["GC_F"], res["GC_R"])
    row["forward_3prime_gc"] = f[-1] in "GC"; row["reverse_3prime_gc"] = r[-1] in "GC"
    row["longest_homopolymer"] = max(max_run(f), max_run(r))


def run_step(run):
    """Returns (tables {key: (DataFrame, file name, description)}, summary, status tuple)."""
    cfg = run.cfg; c = cfg["pcr"]; s = c["selection"]
    summ = {"ran": False}
    W = run.t.get("haplotypes", {}).get("genotypes")
    dls = [m for m in run.sets if m.klass == "DEL"]
    if W is None or W.shape[1] < 2:
        return {}, summ, ("no", "", "skipped", "no structural-haplotype genotypes (step 06 produced none: give --inversions or --haplotypes, or let the scan find intervals)")
    if not dls:
        return {}, summ, ("no", "", "skipped", "no deletion genotypes (give a deletion VCF or matrix with --gt)")
    dels = max(dls, key=lambda m: m.n)
    cons = dict(dels.cons or {})
    extra_cons = getattr(run, "consensus_extra", None) or {}
    for k, v in extra_cons.items():
        cons.setdefault(k, v)
    fa = getattr(run, "reference", None)
    if not cons and fa is None:
        return {}, summ, ("no", "", "skipped", "no sequence to design primers on: the deletion input carries no split-read consensus (INFO CONSENSUS / CONSBP, as written by DELLY) "
                                                "and no --reference FASTA was given")
    engine, msg = resolve_engine(cfg, getattr(run.a, "primer_engine", None))
    if msg:
        run.note(msg)
    names = [x for x in W.columns if x != "individual"]
    iv = pd.DataFrame([{"interval": n, "chrom": p[0], "start": p[1], "end": p[2]} for n in names for p in [parse_interval(n)] if p is not None and p[2] > 0])
    A = W.individual.astype(str).tolist()
    geno = {n: W[n].to_numpy(np.int8) for n in names}
    if not len(iv):
        return {}, summ, ("no", "", "skipped", "the structural-haplotype genotypes carry no interval coordinates (names of the form chrom:start-end are needed to find the deletions inside)")
    tags = run.t.get("haplotypes", {}).get("tags")
    S = screen(iv, geno, A, dels, cfg, tags[tags.svtype == "DEL"] if tags is not None and len(tags) else None)
    picks = read_picks(run.a.pcr_picks) if getattr(run.a, "pcr_picks", None) else None
    npi = int(c["markers_per_interval"]); flank = int(c["reference_flank_bp"])

    def sequence(chrom, start, end, vid):
        v = cons.get((chrom, int(start), vid))
        if v is not None and isinstance(v[0], str) and len(v[0]) >= 60:
            return v[0], int(v[1]), "split-read consensus (SV VCF)"
        t, bp = junction_from_reference(fa, chrom, start, end, flank)
        if t is not None:
            return t, bp, f"reference flanks ({flank} bp on each side)"
        return None, None, ""

    prim, per_interval = [], []
    by_iv = {k: g for k, g in S.groupby("interval", sort=False)} if len(S) else {}
    cls_n = {n: (int((geno[n] == 0).sum()), int((geno[n] == 1).sum()), int((geno[n] == 2).sum())) for n in iv.interval}
    k_marker = 0
    for r in iv.itertuples():
        g = by_iv.get(r.interval)
        row = {"interval": r.interval, "chrom": r.chrom, "interval_start": r.start, "interval_end": r.end, "n_class_0": cls_n[r.interval][0],
               "n_heterokaryotypes": cls_n[r.interval][1], "n_class_2": cls_n[r.interval][2],
               "n_deletions_screened": 0 if g is None else int(len(g)), "n_diagnostic": 0 if g is None else int(g.diagnostic.sum()),
               "n_diagnostic_gel_resolvable": 0 if g is None else int((g.diagnostic & g.gel_resolvable).sum()),
               "n_deletion_on_class_0": 0 if g is None else int((g.diagnostic & g.gel_resolvable & (g.arrangement_with_deletion == "class 0")).sum()),
               "n_deletion_on_class_2": 0 if g is None else int((g.diagnostic & g.gel_resolvable & (g.arrangement_with_deletion == "class 2")).sum()),
               "n_primer_pairs": 0, "reason_no_marker": ""}
        made = 0; tried = 0; statuses = []
        if g is not None:
            if picks is not None:
                want = picks.get(r.interval, [])
                cand = g[g.variant_id.isin(want)]
                missing = [w for w in want if w not in set(g.variant_id)]
            else:
                cand = g[g.diagnostic & g.gel_resolvable]; missing = []
            for idx, q in cand.iterrows():
                if made >= npi and picks is None:
                    break
                tried += 1
                t, bp, src = sequence(q.chrom, q.start, q.end, q.variant_id)
                S.at[idx, "sequence_source"] = src
                if t is None:
                    S.at[idx, "design_status"] = "no sequence"; statuses.append("no sequence"); continue
                res = design(t, bp, q.length_bp, engine, cfg)
                S.at[idx, "design_status"] = res["status"]; statuses.append(res["status"])
                if res["status"] != "ok":
                    continue
                made += 1; k_marker += 1; S.at[idx, "selected"] = True
                p = {"marker": f"M{k_marker:02d}", **{k: q[k] for k in ("interval", "variant_id", "chrom", "start", "end", "length_bp", "arrangement_with_deletion",
                                                                           "deletion_frequency_class_0", "deletion_frequency_heterokaryotypes", "deletion_frequency_class_2",
                                                                           "frequency_difference", "carrier_frequency_class_0", "carrier_frequency_class_2")},
                     "forward_primer": res["fwd"], "reverse_primer": res["rev"], "length_forward": res["len_F"], "length_reverse": res["len_R"],
                     "tm_forward": res["Tm_F"], "tm_reverse": res["Tm_R"], "gc_forward_percent": res["GC_F"], "gc_reverse_percent": res["GC_R"],
                     "product_with_deletion_bp": res["band_deleted"], "product_without_deletion_bp": res["band_undeleted"],
                     "product_difference_bp": res["band_undeleted"] - res["band_deleted"], "design_source": src, "engine": res["engine"], "design_sequence": t,
                     "breakpoint_index": bp, "forward_start_index": res["fpos"], "reverse_end_index": res["rpos"]}
                _checks(p, res, t, s)
                p.update(thermo(res["fwd"], res["rev"], c["tm"]))
                prim.append(p)
        row["n_primer_pairs"] = made
        if not made:
            if g is None:
                n0, _n1, n2 = cls_n[r.interval]
                row["reason_no_marker"] = (f"a homokaryotype class has fewer than {s['min_called_per_class']} individuals" if min(n0, n2) < s["min_called_per_class"]
                                           else f"no deletion inside the interval is genotyped in at least {s['min_called_per_class']} individuals of each homokaryotype class")
            elif picks is not None:
                row["reason_no_marker"] = ("no deletion listed for this interval in the table of picks" if not picks.get(r.interval) else
                                           ("listed deletion not among the screened deletions: " + ", ".join(missing) if missing and not tried else "design failed: " + ", ".join(sorted(set(statuses)))))
            elif not row["n_diagnostic"]:
                row["reason_no_marker"] = f"no diagnostic deletion (frequency difference >= {s['min_frequency_difference']} with heterokaryotypes at {s['het_frequency_min']}-{s['het_frequency_max']})"
            elif not row["n_diagnostic_gel_resolvable"]:
                row["reason_no_marker"] = f"no diagnostic deletion of {s['min_deletion_bp']}-{s['max_deletion_bp']} bp (gel-resolvable)"
            else:
                row["reason_no_marker"] = f"no primer pair for any of the {tried} diagnostic gel-resolvable deletions: " + ", ".join(f"{k} ({statuses.count(k)})" for k in sorted(set(statuses)))
        per_interval.append(row)
    P = pd.DataFrame(prim, columns=PRIMER_COLS)
    for col in THERMO_COLS:
        if col in P and not HAVE_PRIMER3:
            P[col] = np.nan
    # ---- in-silico PCR against the reference
    ins_ran = False
    if fa is not None and len(P) and c["in_silico"]["run"]:
        log(f"  in-silico PCR of {len(P)} primer pairs against --reference ({len(fa.index)} sequences)")
        res = in_silico_pcr(fa, [{"pair": q.marker, "chrom": q.chrom, "start": int(q.start), "end": int(q.end), "fwd": q.forward_primer, "rev": q.reverse_primer,
                                  "expected": int(q.product_without_deletion_bp)} for q in P.itertuples()], c["in_silico"])
        R = pd.DataFrame(res).set_index("pair")
        P["in_silico_products_at_locus"] = P.marker.map(R.products_at_locus).astype(int)
        P["in_silico_product_at_locus_bp"] = P.marker.map(R.product_at_locus_bp)
        P["in_silico_product_minus_design_bp"] = P.in_silico_product_at_locus_bp - P.product_without_deletion_bp
        P["in_silico_products_off_target"] = P.marker.map(R.products_off_target).astype(int)
        P["in_silico_off_target_products"] = P.marker.map(R.off_target_products)
        P["check_in_silico_product_within_5bp"] = (P.in_silico_products_at_locus > 0) & ((P.in_silico_product_minus_design_bp).abs() <= 5)
        ins_ran = True
    S = S.drop(columns=["_d"]) if len(S) else pd.DataFrame(columns=SCREEN_COLS)
    # ---- one table per interval: the marker counts joined with the primers
    I = pd.DataFrame(per_interval)
    keep = ["marker", "variant_id", "start", "end", "length_bp", "arrangement_with_deletion", "deletion_frequency_class_0", "deletion_frequency_heterokaryotypes",
            "deletion_frequency_class_2", "frequency_difference", "forward_primer", "reverse_primer", "tm_forward", "tm_reverse", "gc_forward_percent",
            "gc_reverse_percent", "product_with_deletion_bp", "product_without_deletion_bp", "design_source", "engine", "in_silico_products_at_locus",
            "in_silico_product_at_locus_bp", "in_silico_products_off_target"]
    C = I.merge(P[["interval"] + keep].rename(columns={"start": "deletion_start", "end": "deletion_end"}), on="interval", how="left")
    C = C[[x for x in C.columns if x != "reason_no_marker"] + ["reason_no_marker"]]
    for col in ("marker", "variant_id", "arrangement_with_deletion", "forward_primer", "reverse_primer", "design_source", "engine"):
        C[col] = C[col].fillna("")
    O = pd.DataFrame([{"name": f"{q.marker}_{w}", "sequence": seq, "marker": q.marker, "interval": q.interval, "variant_id": q.variant_id, "length": len(seq), "tm": t}
                      for q in P.itertuples() for w, seq, t in (("F", q.forward_primer, q.tm_forward), ("R", q.reverse_primer, q.tm_reverse))], columns=ORDER_COLS)
    tmc = tm_conditions(c["tm"])
    src_counts = P.design_source.value_counts().to_dict() if len(P) else {}
    eng_txt = ("Primer3" if engine == "primer3" else "built-in routine")
    new_txt = "Hairpin and dimer Tm / dG are new in BroodLink (primer3-py thermodynamic alignment), empty without primer3-py."
    tables = {
        "candidates": (S, "pcr_candidate_deletions.tsv", "Screened deletions: one row per interval and deletion lying inside it and genotyped in at least "
                       f"{s['min_called_per_class']} individuals of each homokaryotype class. Deletion allele frequency per class, their difference, the screening bands of a generic pair "
                       f"(2 x {s['screen_flank_bp']} bp of flank), diagnostic (difference >= {s['min_frequency_difference']}, heterokaryotypes at {s['het_frequency_min']}-{s['het_frequency_max']}), "
                       f"gel_resolvable (deletion of {s['min_deletion_bp']}-{s['max_deletion_bp']} bp), rank inside the interval (largest difference first), the status of the primer design "
                       "where it was tried, and selected."),
        "primers": (P, "pcr_primers.tsv", f"Designed primer pairs, one row per pair ({eng_txt}). {tmc}. Products: with the deletion = distance between the primers across the junction; "
                    "without = that plus the deletion length. design_source says which sequence each pair was designed on. Checks: primers are exact substrings of the design sequence and flank "
                    "the breakpoint; products recomputed from the primer positions; in-silico PCR against --reference when given (<= "
                    f"{c['in_silico']['max_mismatches']} mismatches per primer, 3' base matching, products {c['in_silico']['product_min']}-{c['in_silico']['product_max']} bp, at the locus = within "
                    f"{c['in_silico']['locus_window_bp']} bp of the deletion). {new_txt} In-silico design only: no laboratory test."),
        "per_interval": (C, "pcr_markers_per_interval.tsv", "One row per interval (and per marker when several are designed): individuals per class, deletions screened, diagnostic, diagnostic and "
                         f"gel-resolvable, on which class, the selected deletion with its frequencies, the primer pair with Tm, GC and both product sizes, and the reason when an interval has no marker. {tmc}. "
                         "In-silico design only: no laboratory test."),
        "order": (O, "pcr_primer_order_list.tsv", f"Primer list ready to order: name (marker_F / marker_R), sequence 5'-3'. {tmc}. Not tested in a laboratory."),
    }
    n_iv = int(len(I)); n_with = int((I.n_primer_pairs > 0).sum())
    summ = {"ran": True, "engine": engine, "engine_requested": (getattr(run.a, "primer_engine", None) or c.get("engine") or "primer3"), "primer3_installed": HAVE_PRIMER3, "primer3_version": PRIMER3_VERSION,
            "tm_conditions": tmc, "n_intervals": n_iv, "n_intervals_with_marker": n_with, "n_intervals_without_marker": n_iv - n_with,
            "intervals_without_marker": {r.interval: r.reason_no_marker for r in I.itertuples() if not r.n_primer_pairs},
            "n_deletions_screened": int(len(S)), "n_diagnostic": int(S.diagnostic.sum()) if len(S) else 0,
            "n_diagnostic_gel_resolvable": int((S.diagnostic & S.gel_resolvable).sum()) if len(S) else 0, "n_primer_pairs": int(len(P)),
            "design_sources": src_counts, "picks_supplied": picks is not None, "markers_per_interval": npi,
            "reference_given": fa is not None, "reference_index": getattr(fa, "indexed_from", "") if fa is not None else "", "in_silico_ran": ins_ran,
            "n_consensus_available": int(len(cons)), "deletion_input": dels.name}
    if len(P):
        summ.update({"tm_min": float(min(P.tm_forward.min(), P.tm_reverse.min())), "tm_max": float(max(P.tm_forward.max(), P.tm_reverse.max())),
                     "length_min": int(min(P.length_forward.min(), P.length_reverse.min())), "length_max": int(max(P.length_forward.max(), P.length_reverse.max())),
                     "product_with_deletion_min": int(P.product_with_deletion_bp.min()), "product_with_deletion_max": int(P.product_with_deletion_bp.max()),
                     "n_checks_exact": int(P.check_primers_exact_in_design_sequence.sum()), "n_checks_flank": int(P.check_primers_flank_breakpoint.sum()),
                     "n_checks_products": int(P.check_products_recomputed.sum())})
        if ins_ran:
            summ.update({"n_in_silico_product_at_locus": int((P.in_silico_products_at_locus > 0).sum()), "n_in_silico_off_target": int((P.in_silico_products_off_target > 0).sum()),
                         "n_in_silico_within_5bp": int(P.check_in_silico_product_within_5bp.sum())})
    src_txt = "; ".join(f"{k}: {v} pairs" for k, v in src_counts.items()) or "no pair designed"
    run.note(f"PCR markers: {len(P)} primer pairs for {n_with} of {n_iv} intervals, designed with {eng_txt} on: {src_txt}. In-silico design only; no laboratory test.")
    conf = "benchmarked" if engine == "primer3" else "provisional"
    why = ((f"Primer3 with the settings recovered from the published primer table of the original study (33 of its 35 pairs returned exactly); " if engine == "primer3" else
            "built-in routine (primer3-py not installed or engine set to builtin): not the published procedure; ") +
           "screening rules of the published screening table; in-silico design only, no laboratory test")
    return tables, summ, ("yes", dels.name + ("; --reference" if fa is not None else ""), conf, why)
