"""Input readers with automatic detection of file type and marker class.

Everything is read with numpy / pandas only (cyvcf2 or pysam are not required).
Genotypes are stored as int8, markers x samples, -1 = missing.
"""
from __future__ import annotations
import glob
import gzip
import io as _io
import os
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

GT_CODE = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/0": 1, "1/1": 2, "1|1": 2}
SV_CLASSES = ("DEL", "DUP", "INV", "INS")


class InputError(SystemExit):
    """An input cannot be read or matched: the only reason QC stops a run."""

    def __init__(self, msg):
        super().__init__(f"ERROR: {msg}")


import time as _time

# 0 = errors and the final line only, 1 = step headers and results (--quiet),
# 2 = also details and progress of long loops (default)
VERBOSE = 2
_T0 = _time.time()


def set_verbosity(level):
    global VERBOSE, _T0
    VERBOSE = int(level); _T0 = _time.time()


def _clock():
    s = int(_time.time() - _T0)
    return f"{s // 60:d}:{s % 60:02d}"


def log(msg, level=1):
    """Print one line with the time since the start of the run. level 1 = always (unless
    silenced), level 2 = detail, hidden with --quiet."""
    if level <= VERBOSE or str(msg).startswith(("WARNING", "ERROR")):
        print(f"[broodlink {_clock()}] {msg}", flush=True)


def detail(msg):
    log("    " + str(msg), level=2)


class Progress:
    """Progress of a long loop, printed about every 10% (detail level), with the time left.
    Nothing is printed until the loop has run for `quiet_s` seconds, so loops that finish
    quickly stay silent."""

    def __init__(self, label, total, min_total=20, quiet_s=2.0):
        self.label, self.total, self.quiet_s = label, int(total), float(quiet_s)
        self.on = VERBOSE >= 2 and self.total >= min_total
        self.t0 = _time.time(); self.next = 0.1

    def update(self, i):
        if not self.on:
            return
        f = (i + 1) / self.total
        if f >= self.next - 1e-9 or i + 1 == self.total:
            el = _time.time() - self.t0
            while self.next <= f + 1e-9:
                self.next += 0.1
            if el < self.quiet_s:                 # short loops print nothing at all
                return
            left = el / f * (1 - f)
            log(f"    {self.label}: {i + 1:,}/{self.total:,} ({100 * f:.0f}%)" + (f", about {left:.0f} s left" if f < 1 else f", {el:.1f} s"), level=2)


def _open(path):
    if str(path).endswith((".gz", ".bgz")):
        return _io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


@dataclass
class MarkerSet:
    name: str                 # label used in every table, e.g. "DEL (gt:calls.vcf)"
    klass: str                # SNP, SNP_INDEL, DEL, DUP, INV, INS, SV, HAPLOTYPE
    flag: str                 # command-line input it came from
    source: str               # file name (base name only, never a full path)
    samples: list
    chrom: np.ndarray
    pos: np.ndarray
    end: np.ndarray
    ids: np.ndarray
    G: np.ndarray             # int8 markers x samples, -1 missing
    D: np.ndarray | None = None   # uint8 expected dosage x 100 (255 missing), only from likelihoods
    filters: list = field(default_factory=list)   # (step, markers_before, markers_after)
    notes: list = field(default_factory=list)
    cons: dict = field(default_factory=dict)      # (chrom, pos, id) -> (split-read consensus, breakpoint index); DELLY INFO CONSENSUS / CONSBP

    @property
    def n(self):
        return self.G.shape[0]

    def subset_markers(self, keep, step=None):
        keep = np.asarray(keep)
        before = self.n
        out = MarkerSet(self.name, self.klass, self.flag, self.source, list(self.samples), self.chrom[keep],
                        self.pos[keep], self.end[keep], self.ids[keep], self.G[keep],
                        None if self.D is None else self.D[keep], list(self.filters), list(self.notes), self.cons)
        if step:
            out.filters.append((step, before, out.n))
        return out

    def subset_samples(self, samples):
        ix = {s: i for i, s in enumerate(self.samples)}
        cols = [ix[s] for s in samples if s in ix]
        kept = [s for s in samples if s in ix]
        return MarkerSet(self.name, self.klass, self.flag, self.source, kept, self.chrom, self.pos, self.end,
                         self.ids, self.G[:, cols], None if self.D is None else self.D[:, cols],
                         list(self.filters), list(self.notes), self.cons)

    def dosage(self):
        """float32 dosage, NaN missing: expected dosage when likelihoods exist, else the genotype."""
        if self.D is not None:
            d = self.D.astype(np.float32) / 100.0
            d[self.D == 255] = np.nan
            return d
        d = self.G.astype(np.float32)
        d[self.G < 0] = np.nan
        return d


def allele_freq(G):
    called = G >= 0
    n = called.sum(1)
    return np.where(n > 0, np.where(called, G, 0).sum(1) / (2.0 * np.maximum(n, 1)), np.nan)


# --------------------------------------------------------------------------- detection
def sniff(path):
    """Return one of: beagle, vcf, gt_matrix, unknown."""
    with _open(path) as fh:
        first = fh.readline()
        if first.startswith("##fileformat=VCF") or first.startswith("#CHROM"):
            return "vcf"
        f = first.rstrip("\n").split("\t")
        if len(f) >= 6 and f[0].lower() == "marker" and f[1].lower().startswith("allele"):
            return "beagle"
        if len(f) >= 4 and f[0].upper() in ("CHROM", "#CHROM") and f[1].upper() == "POS":
            return "gt_matrix"
    return "unknown"


def expand(path_or_list):
    """A file, a folder or a glob -> sorted list of files."""
    out = []
    for p in ([path_or_list] if isinstance(path_or_list, str) else list(path_or_list or [])):
        if os.path.isdir(p):
            out += sorted(f for f in glob.glob(os.path.join(p, "*")) if os.path.isfile(f))
        elif any(c in p for c in "*?["):
            out += sorted(glob.glob(p))
        else:
            out.append(p)
    for p in out:
        if not os.path.exists(p):
            raise InputError(f"input file not found: {p}")
    return out


# --------------------------------------------------------------------------- beagle
def read_beagle(paths, sample_names=None, step=1, hard=0.90, flag="bisnp", max_markers=None):
    """ANGSD beagle genotype likelihoods (marker, allele1, allele2, 3 columns per individual)."""
    Gs, Ds, chroms, poss = [], [], [], []
    names = None
    for path in paths:
        if sniff(path) != "beagle":
            continue
        with _open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
        n = (len(header) - 3) // 3
        hn = [header[3 + 3 * i] for i in range(n)]
        if names is None:
            names = hn
        elif hn != names:
            raise InputError(f"beagle files disagree on the individuals: {os.path.basename(path)}")
        skip = (lambda i: i != 0 and i % step != 1) if step > 1 else None
        for B in pd.read_csv(path, sep="\t", chunksize=20000, skiprows=skip, engine="c"):
            mk = B.iloc[:, 0].astype(str)
            sp = mk.str.rsplit("_", n=1)
            gl = B.iloc[:, 3:].to_numpy(dtype=np.float32).reshape(len(B), n, 3)
            s = gl.sum(2, keepdims=True)
            flat = (gl.max(2) - gl.min(2)) < 1e-6           # uninformative (no reads)
            gl = gl / np.where(s > 0, s, 1)
            G = np.where((gl.max(2) > hard) & ~flat, gl.argmax(2), -1).astype(np.int8)
            dos = gl[:, :, 1] + 2 * gl[:, :, 2]
            D = np.where(flat, 255, np.clip(np.rint(dos * 100), 0, 200)).astype(np.uint8)
            Gs.append(G); Ds.append(D); chroms.append(sp.str[0].to_numpy()); poss.append(sp.str[1].astype(np.int64).to_numpy())
    if names is None:
        raise InputError("no beagle genotype-likelihood file could be read from --bisnp")
    G = np.vstack(Gs); D = np.vstack(Ds); chrom = np.concatenate(chroms); pos = np.concatenate(poss)
    if sample_names is not None and len(sample_names) == len(names):
        names = list(sample_names)
    ms = MarkerSet("SNP (bisnp)", "SNP", flag, os.path.basename(paths[0]) + (f" (+{len(paths) - 1} files)" if len(paths) > 1 else ""),
                   list(names), chrom, pos, pos.copy(), np.char.add(np.char.add(chrom.astype(str), "_"), pos.astype(str)), G, D)
    ms.notes.append(f"genotype likelihoods; hard call when the best likelihood > {hard}; every {step} SNP(s) read")
    if max_markers and ms.n > max_markers:
        keep = np.zeros(ms.n, bool); keep[np.linspace(0, ms.n - 1, max_markers).astype(int)] = True
        ms = ms.subset_markers(keep, "marker cap")
    return [ms]


# --------------------------------------------------------------------------- VCF
_INFO_END = re.compile(r"(?:^|;)END=(\d+)")
_INFO_SVTYPE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
_INFO_SVLEN = re.compile(r"(?:^|;)SVLEN=(-?\d+)")
_INFO_CONS = re.compile(r"(?:^|;)CONSENSUS=([ACGTNacgtn]+)")
_INFO_CONSBP = re.compile(r"(?:^|;)CONSBP=(\d+)")


def _gt_code(s):
    c = GT_CODE.get(s)
    if c is not None:
        return c
    if "." in s:
        return -1
    a = re.split(r"[/|]", s)
    try:
        return min(2, sum(int(x) > 0 for x in a)) if len(a) == 2 else -1
    except ValueError:
        return -1


def _lik_from_field(val, kind):
    try:
        x = np.array([float(v) for v in val.split(",")[:3]])
    except ValueError:
        return None
    if len(x) != 3:
        return None
    x = x - x.max() if kind == "GL" else -(x - x.min()) / 10.0
    p = np.power(10.0, x)
    return p / p.sum()


def read_vcf(path, flag, hard=0.90, want_lik=False):
    """Any VCF / gVCF with GT (or GL / PL when want_lik). Returns one MarkerSet per marker class."""
    base = os.path.basename(path)
    rows = {}            # klass -> lists
    samples = None
    n_multi = n_ref = 0
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            if samples is None or len(f) < 10:
                continue
            alt = f[4]
            if alt in (".", "<NON_REF>", "<*>"):
                n_ref += 1
                continue
            alts = [a for a in alt.split(",") if a not in ("<NON_REF>", "<*>")]
            if len(alts) != 1:
                n_multi += 1
                continue
            info = f[7]
            m = _INFO_SVTYPE.search(info)
            pos = int(f[1])
            if m:
                k = m.group(1).upper()
                k = "DUP" if k.startswith("DUP") else k
                k = k if k in SV_CLASSES else "SV"
                e = _INFO_END.search(info)
                if e:
                    end = int(e.group(1))
                else:
                    sl = _INFO_SVLEN.search(info)
                    end = pos + abs(int(sl.group(1))) if sl else pos
            elif alts[0].startswith("<"):
                k = alts[0].strip("<>").split(":")[0].upper()
                k = k if k in SV_CLASSES else "SV"
                e = _INFO_END.search(info)
                end = int(e.group(1)) if e else pos
            else:
                k = "SNP" if (len(f[3]) == 1 and len(alts[0]) == 1) else "INDEL"
                end = pos + len(f[3]) - 1
            fmt = f[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else -1
            li, lkind = (-1, None)
            if want_lik:
                for kind in ("GL", "PL"):
                    if kind in fmt:
                        li, lkind = fmt.index(kind), kind
                        break
            g = np.full(len(samples), -1, np.int8)
            d = np.full(len(samples), 255, np.uint8) if li >= 0 else None
            for i, cell in enumerate(f[9:]):
                parts = cell.split(":")
                if li >= 0 and li < len(parts):
                    p = _lik_from_field(parts[li], lkind)
                    if p is not None and (p.max() - p.min()) > 1e-6:
                        d[i] = int(round(100 * (p[1] + 2 * p[2])))
                        g[i] = int(p.argmax()) if p.max() > hard else -1
                        continue
                if li < 0 and gi >= 0 and gi < len(parts):
                    g[i] = _gt_code(parts[gi])
            r = rows.setdefault(k, {"chrom": [], "pos": [], "end": [], "id": [], "G": [], "D": [], "cons": {}})
            r["chrom"].append(f[0]); r["pos"].append(pos); r["end"].append(end)
            r["id"].append(f[2] if f[2] != "." else f"{f[0]}_{pos}_{k}"); r["G"].append(g)
            if m and "CONSENSUS=" in info:                       # split-read consensus of the junction (DELLY): kept for primer design
                cs, cb = _INFO_CONS.search(info), _INFO_CONSBP.search(info)
                if cs and cb:
                    r["cons"][(f[0], pos, r["id"][-1])] = (cs.group(1).upper(), int(cb.group(1)))
            if d is not None:
                r["D"].append(d)
    if samples is None:
        raise InputError(f"no #CHROM header line in VCF {base}")
    out = []
    if flag == "snp_indel" and ("SNP" in rows or "INDEL" in rows):       # small variants kept as one class
        merged = {"chrom": [], "pos": [], "end": [], "id": [], "G": [], "D": []}
        for k in ("SNP", "INDEL"):
            for key in merged:
                merged[key] += rows.get(k, {}).get(key, [])
        rows = {k: v for k, v in rows.items() if k not in ("SNP", "INDEL")}
        rows["SNP_INDEL"] = merged
    for k, r in rows.items():
        if not r["G"]:
            continue
        o = np.lexsort((np.array(r["pos"]), np.array(r["chrom"])))
        D = np.vstack(r["D"])[o] if r["D"] and len(r["D"]) == len(r["G"]) else None
        ms = MarkerSet(f"{k} ({flag}:{base})", k, flag, base, list(samples), np.array(r["chrom"])[o],
                       np.array(r["pos"], np.int64)[o], np.array(r["end"], np.int64)[o],
                       np.array(r["id"])[o], np.vstack(r["G"])[o], D)
        ms.cons = r.get("cons", {})
        if ms.cons:
            ms.notes.append(f"{len(ms.cons):,} records carry a split-read consensus (INFO CONSENSUS / CONSBP)")
        if n_multi:
            ms.notes.append(f"{n_multi} multi-allelic records of the file skipped")
        if n_ref:
            ms.notes.append(f"{n_ref} reference-only (gVCF) records of the file skipped")
        out.append(ms)
    if not out:
        raise InputError(f"no usable variant record with genotypes in {base}")
    return out


def read_consensus(path):
    """Only the split-read consensus of a (large) SV VCF: {(chrom, pos, id): (sequence, breakpoint index)}.
    Genotype columns are not parsed."""
    out = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or "CONSENSUS=" not in line:
                continue
            f = line.split("\t", 8)
            cs, cb = _INFO_CONS.search(f[7]), _INFO_CONSBP.search(f[7])
            if cs and cb:
                out[(f[0], int(f[1]), f[2])] = (cs.group(1).upper(), int(cb.group(1)))
    return out


def read_gt_matrix(path, flag):
    """Genotype matrix TSV: CHROM POS [END] [ID] [SVTYPE] then one column per sample (0/0, 0/1, 1/1 or 0/1/2)."""
    base = os.path.basename(path)
    d = pd.read_csv(path, sep="\t", dtype=str)
    d.columns = [c.lstrip("#") for c in d.columns]
    meta = [c for c in d.columns if c.upper() in ("CHROM", "POS", "END", "ID", "SVTYPE", "SVLEN", "REF", "ALT")]
    samples = [c for c in d.columns if c not in meta]
    up = {c.upper(): c for c in meta}
    raw = d[samples].to_numpy()
    G = np.full(raw.shape, -1, np.int8)
    for k, v in list(GT_CODE.items()) + [("0", 0), ("1", 1), ("2", 2)]:
        G[raw == k] = v
    pos = d[up["POS"]].astype(np.int64).to_numpy()
    end = d[up["END"]].astype(np.int64).to_numpy() if "END" in up else pos.copy()
    chrom = d[up["CHROM"]].to_numpy()
    ids = d[up["ID"]].to_numpy() if "ID" in up else np.array([f"{c}_{p}" for c, p in zip(chrom, pos)])
    if "SVTYPE" in up:
        kl = d[up["SVTYPE"]].str.upper().to_numpy()
    else:
        m = re.search(r"(DEL|DUP|INV|INS)", base.upper())
        kl = np.array([m.group(1) if m else ("SV" if "END" in up else "SNP")] * len(d))
    out = []
    for k in pd.unique(kl):
        sel = kl == k
        kk = k if k in SV_CLASSES + ("SNP",) else "SV"
        out.append(MarkerSet(f"{kk} ({flag}:{base})", kk, flag, base, list(samples), chrom[sel], pos[sel], end[sel], ids[sel], G[sel]))
    return out


def read_variants(path, flag, cfg, sample_names=None):
    kind = sniff(path)
    hard = cfg["genotypes"]["gl_hard_call"]
    if kind == "beagle":
        return read_beagle([path], sample_names, cfg["genotypes"]["bisnp_step"], hard, flag, cfg["genotypes"]["max_markers_snp"])
    if kind == "vcf":
        return read_vcf(path, flag, hard, want_lik=(flag == "bisnp"))
    if kind == "gt_matrix":
        return read_gt_matrix(path, flag)
    raise InputError(f"file type of {os.path.basename(path)} not recognised (expected beagle, VCF or a genotype matrix TSV)")


# --------------------------------------------------------------------------- haplotypes
def read_haplotypes(path, known_samples=None):
    """Structural-haplotype genotypes 0/1/2: wide (samples x intervals or intervals x samples) or long."""
    base = os.path.basename(path)
    d = pd.read_csv(path, sep="\t", dtype=str)
    low = [c.lower() for c in d.columns]
    if d.shape[1] == 3 and any(c in low for c in ("interval", "interval_id", "haplotype")) :
        ic = low.index([c for c in low if c in ("interval", "interval_id", "haplotype")][0])
        sc = [i for i, c in enumerate(low) if c in ("sample", "sample_id", "id", "individual")][0]
        gc = [i for i in range(3) if i not in (ic, sc)][0]
        w = d.pivot(index=d.columns[ic], columns=d.columns[sc], values=d.columns[gc])
    else:
        w = d.set_index(d.columns[0])
        ks = set(known_samples or [])
        if ks and len(ks & set(w.index)) >= len(ks & set(w.columns)):
            w = w.T                                    # rows were samples -> make rows intervals
        elif not ks and w.shape[0] > w.shape[1]:
            w = w.T
    raw = w.to_numpy()
    G = np.full(raw.shape, -1, np.int8)
    for k, v in list(GT_CODE.items()) + [("0", 0), ("1", 1), ("2", 2), ("0.0", 0), ("1.0", 1), ("2.0", 2)]:
        G[raw == k] = v
    ids = np.array([str(x) for x in w.index])
    chrom, pos, end = [], [], []
    for x in ids:
        m = re.match(r"(.+):(\d+)-(\d+)$", x)
        chrom.append(m.group(1) if m else "NA"); pos.append(int(m.group(2)) if m else 0); end.append(int(m.group(3)) if m else 0)
    return MarkerSet(f"HAPLOTYPE (haplotypes:{base})", "HAPLOTYPE", "haplotypes", base, [str(s) for s in w.columns],
                     np.array(chrom), np.array(pos, np.int64), np.array(end, np.int64), ids, G)


def read_intervals_bed(path):
    """--inversions: chrom, start, end, name, optional tier. No sample information is read."""
    rows = []
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[1].lower() in ("start", "chromstart"):
                continue
            name = f[3] if len(f) > 3 and f[3] else f"{f[0]}:{f[1]}-{f[2]}"
            rows.append({"interval": name, "chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                         "tier_supplied": f[4] if len(f) > 4 else ""})
    if not rows:
        raise InputError(f"no interval in {os.path.basename(path)}")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- relatedness
def _std_rel(df, src):
    df = df.copy()
    for c in ("k0", "k1", "k2"):
        if c not in df:
            df[c] = np.nan
    df["id1"] = df["id1"].astype(str); df["id2"] = df["id2"].astype(str)
    df["source"] = src
    # R0 (when the table has it, as NGSrelate output does) is kept for the opposite-homozygote panel
    return df[["id1", "id2", "theta", "k0", "k1", "k2", "source"] + (["R0"] if "R0" in df else [])]


def read_relatedness_file(path, sample_names=None):
    base = os.path.basename(path)
    d = pd.read_csv(path, sep=r"\s+" if path.endswith((".kin", ".kin0")) else "\t")
    cols = {c.lower(): c for c in d.columns}
    if "theta" in cols and ("j7" in cols or "ida" in cols or ("a" in cols and "b" in cols)):
        if "ida" in cols:
            id1, id2 = d[cols["ida"]], d[cols["idb"]]
        elif sample_names is not None:
            id1 = d[cols["a"]].map(lambda i: sample_names[int(i)]); id2 = d[cols["b"]].map(lambda i: sample_names[int(i)])
        else:
            id1, id2 = "Ind" + d[cols["a"]].astype(str), "Ind" + d[cols["b"]].astype(str)
        out = pd.DataFrame({"id1": id1, "id2": id2, "theta": d[cols["theta"]]})
        if "j7" in cols:
            out["k0"], out["k1"], out["k2"] = d[cols["j9"]], d[cols["j8"]], d[cols["j7"]]
        if "r0" in cols:
            out["R0"] = d[cols["r0"]]
        return _std_rel(out, "NGSrelate"), "NGSrelate .res"
    if "kinship" in cols and ("id1" in cols and "id2" in cols):
        out = pd.DataFrame({"id1": d[cols["id1"]], "id2": d[cols["id2"]], "theta": d[cols["kinship"]]})
        for k in ("k0", "k1", "k2"):
            if k in cols:
                out[k] = d[cols[k]]
        kind = "KING .kin0" if "fid2" in cols else ("KING .kin" if "fid" in cols or "fid1" in cols else "plain TSV (kinship)")
        return _std_rel(out, "KING" if kind.startswith("KING") else "table"), kind
    if "theta" in cols and "id1" in cols and "id2" in cols:
        out = pd.DataFrame({"id1": d[cols["id1"]], "id2": d[cols["id2"]], "theta": d[cols["theta"]]})
        for k in ("k0", "k1", "k2"):
            if k in cols:
                out[k] = d[cols[k]]
        if "r0" in cols:
            out["R0"] = d[cols["r0"]]
        return _std_rel(out, "table"), "plain TSV"
    raise InputError(f"relatedness file {base} not recognised (expected NGSrelate .res, KING .kin/.kin0 or a TSV with id1, id2, theta or kinship)")


def read_relatedness(path, sample_names=None):
    """A file or a folder of files (concatenated; sample order and the pair set are checked)."""
    files = expand(path)
    files = [f for f in files if not os.path.basename(f).startswith(".")]
    parts, kinds = [], []
    for f in files:
        r, kind = read_relatedness_file(f, sample_names)
        parts.append(r); kinds.append(kind)
    if len(set(kinds)) > 1:
        raise InputError(f"relatedness folder mixes formats: {sorted(set(kinds))}")
    R = pd.concat(parts, ignore_index=True)
    a = np.where(R.id1 < R.id2, R.id1, R.id2); b = np.where(R.id1 < R.id2, R.id2, R.id1)
    R["id1"], R["id2"] = a, b
    if (R.id1 == R.id2).any():
        R = R[R.id1 != R.id2]
    dup = R.duplicated(["id1", "id2"], keep=False)
    if dup.any():
        g = R[dup].groupby(["id1", "id2"]).theta.agg(lambda x: float(x.max() - x.min()))
        if (g > 1e-6).any():
            bad = g[g > 1e-6].index[:5].tolist()
            raise InputError(f"relatedness files are inconsistent: {int((g > 1e-6).sum())} pairs occur twice with different theta, e.g. {bad}")
        R = R.drop_duplicates(["id1", "id2"])
    ids = sorted(set(R.id1) | set(R.id2))
    expected = len(ids) * (len(ids) - 1) // 2
    info = {"format": kinds[0], "n_files": len(files), "n_samples": len(ids), "n_pairs": len(R), "n_pairs_expected": expected,
            "complete": len(R) == expected}
    if len(files) > 1 and len(R) != expected:
        raise InputError(f"relatedness folder does not hold a complete pair set: {len(R)} pairs for {len(ids)} samples (expected {expected}); "
                         "the files do not cover the same sample set")
    return R.reset_index(drop=True), info


# --------------------------------------------------------------------------- samples, depth, genome
def read_samples(path):
    d = pd.read_csv(path, sep="\t", dtype=str, header=None)
    first = [str(x).lower() for x in d.iloc[0]]
    if first[0] in ("id", "sample", "sample_id", "#id", "#sample", "fish_id", "individual"):
        d.columns = ["id"] + [c for c in first[1:]]
        d = d.iloc[1:].reset_index(drop=True)
    else:
        d.columns = ["id"] + [f"col{i}" for i in range(1, d.shape[1])]
    return d


def read_fai(path):
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"], dtype={"chrom": str})
    return dict(zip(d.chrom, d.length.astype(np.int64)))


def read_depth(path):
    """Per-sample depth: a table (sample, depth columns) or a folder of mosdepth summaries / one-number files."""
    rows = []
    if os.path.isdir(path):
        for f in sorted(glob.glob(os.path.join(path, "*"))):
            base = os.path.basename(f)
            sample = re.sub(r"(\.mosdepth)?(\.summary)?(\.txt|\.tsv)?$", "", base)
            try:
                t = pd.read_csv(f, sep="\t")
            except Exception:
                continue
            if {"chrom", "mean"} <= set(t.columns):
                tot = t[t.chrom == "total"]; reg = t[t.chrom == "total_region"]
                rows.append({"id": sample, "depth_assembly": float(tot["mean"].iloc[0]) if len(tot) else np.nan,
                             "depth_mask": float(reg["mean"].iloc[0]) if len(reg) else np.nan})
            elif t.shape[1] >= 3 and t.shape[0] > 0:               # samtools depth: chrom pos depth
                rows.append({"id": sample, "depth_assembly": float(pd.to_numeric(t.iloc[:, 2], errors="coerce").mean()), "depth_mask": np.nan})
        return pd.DataFrame(rows)
    t = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in t.columns}
    out = pd.DataFrame({"id": t.iloc[:, 0].astype(str)})
    mk = [c for c in cols if "mask" in c or "callable" in c]
    asm = [c for c in cols if ("assembly" in c or "genome" in c or c in ("depth", "mean_depth", "mean")) and c not in mk]
    out["depth_mask"] = pd.to_numeric(t[cols[mk[0]]], errors="coerce") if mk else np.nan
    out["depth_assembly"] = pd.to_numeric(t[cols[asm[0]]], errors="coerce") if asm else np.nan
    return out


def write_tsv(df, path, header_lines=()):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for h in header_lines:
            fh.write(f"# {h}\n")
        df.to_csv(fh, sep="\t", index=False, lineterminator="\n", float_format="%.6g")
