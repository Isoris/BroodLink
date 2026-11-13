"""Optional input --genes: gene models of the assembly, and the genes of each barcode interval.

`read_genes` reads a GFF3 (feature type "gene"; mRNA / transcript when the file has no gene
feature) or a BED file (chrom, start, end, name, score, strand). Coordinates are kept in the
convention of the interval tables of BroodLink (`io.read_intervals_bed`): 0-based start, end
exclusive, so a GFF3 start (1-based, closed) is lowered by one. When the GFF3 carries CDS
features, their intervals are kept per gene, so that "a variant hits coding sequence" can be
told from "a variant lies in the span of the gene".

`interval_genes` and `tag_hits` are extracted from the script of the published gene-level
barcode figure (HAT/figures/upgrades_2026-11-03/round_ed/ED_Fig_H6/ED_Fig_H6.py, function
`prepare`, and HAT/code/del_gene_loci_common.py, `load_hits` / `calls_table`):
  * a gene belongs to an interval when it overlaps it (gene end >= interval start and gene
    start <= interval end); the table says whether it lies inside or spans a boundary;
  * a "call" is an arrangement tag variant that hits coding sequence of a gene of the interval
    (the gene span when the gene file has no CDS); calls are numbered in genome order inside
    their interval (<interval number>.<k>).
The original took the hits from a cohort table built outside the pipeline (CDS overlap of
arrangement-private deletions); here they are computed from --genes and the tag variants of
step 6. Gene families and functional annotation of the original have no generic input.
"""
from __future__ import annotations
import os
import re

import numpy as np
import pandas as pd

from .io import InputError, _open

GENE_COLS = ["gene", "name", "chrom", "start", "end", "strand", "biotype"]
CDS_COLS = ["gene", "chrom", "start", "end"]
_BIOTYPE_KEYS = ("gene_biotype", "biotype", "gene_type")
_TRANSCRIPT_TYPES = ("mRNA", "transcript")


def _attrs(s):
    out = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1); out[k.strip()] = v.strip()
    return out


def _is_int(x):
    return bool(re.fullmatch(r"\d+", x.strip()))


def sniff_genes(path):
    """'gff3' or 'bed', from the first data line."""
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##gff-version"):
                return "gff3"
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and _is_int(f[3]) and _is_int(f[4]) and not _is_int(f[1]):
                return "gff3"
            return "bed"
    raise InputError(f"--genes: {os.path.basename(str(path))} is empty")


def _read_gff3(path):
    genes, tx, cds, parent = [], [], [], {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not (_is_int(f[3]) and _is_int(f[4])):
                continue
            typ = f[2]
            if typ == "gene" or typ in _TRANSCRIPT_TYPES:
                a = _attrs(f[8]); fid = a.get("ID", "")
                rec = (fid or a.get("Name", "") or f"{f[0]}:{f[3]}-{f[4]}", a.get("Name") or a.get("gene") or a.get("gene_name") or "", f[0], int(f[3]) - 1, int(f[4]),
                       f[6] if f[6] in "+-" and f[6] else ".", next((a[k] for k in _BIOTYPE_KEYS if k in a), ""))
                if typ == "gene":
                    genes.append(rec)
                else:
                    tx.append(rec)
                    if fid and a.get("Parent"):
                        parent[fid] = a["Parent"].split(",")[0]
            elif typ == "CDS":
                a = _attrs(f[8])
                if a.get("Parent"):
                    cds.append((a["Parent"].split(",")[0], f[0], int(f[3]) - 1, int(f[4])))
    level = "gene"
    if not genes:                                               # no gene feature: transcripts stand in for genes
        genes, level, parent = tx, "mRNA", {}
    G = pd.DataFrame(genes, columns=GENE_COLS)
    ids = set(G.gene)
    C = pd.DataFrame([(parent.get(p, p), c, s, e) for p, c, s, e in cds], columns=CDS_COLS)
    C = C[C.gene.isin(ids)].drop_duplicates().reset_index(drop=True)
    return G, C, level


def _read_bed(path):
    rows = []
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not (_is_int(f[1]) and _is_int(f[2])):
                continue                                        # a header line
            name = f[3] if len(f) > 3 and f[3] else f"{f[0]}:{f[1]}-{f[2]}"
            rows.append((name, name if len(f) > 3 and f[3] else "", f[0], int(f[1]), int(f[2]), f[5] if len(f) > 5 and f[5] in ("+", "-") else ".", ""))
    return pd.DataFrame(rows, columns=GENE_COLS), pd.DataFrame(columns=CDS_COLS), "BED line"


def read_genes(path, with_cds=False):
    """--genes: one row per gene (gene = identifier, name, chrom, start (0-based), end, strand,
    biotype when the file gives one). with_cds: also the CDS intervals per gene (empty for BED)
    and a dict describing what was read."""
    fmt = sniff_genes(path)
    G, C, level = _read_gff3(path) if fmt == "gff3" else _read_bed(path)
    if not len(G):
        raise InputError(f"--genes: no gene found in {os.path.basename(str(path))} (GFF3 with gene or mRNA features, or BED with chrom, start, end, name)")
    G = G[G.end > G.start].drop_duplicates(["gene", "chrom", "start"]).sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    info = {"format": "GFF3" if fmt == "gff3" else "BED", "feature": level, "n_genes": int(len(G)), "n_chrom": int(G.chrom.nunique()),
            "n_genes_with_cds": int(C.gene.nunique()) if len(C) else 0, "has_biotype": bool((G.biotype != "").any())}
    return (G, C, info) if with_cds else G


def drawn_genes(genes):
    """Genes shown and counted by the gene-level tables: the protein-coding genes when the file
    gives a biotype (as the original figure), every gene otherwise."""
    if len(genes) and (genes.biotype == "protein_coding").any():
        return genes[genes.biotype == "protein_coding"], "protein-coding genes (biotype of the gene file)"
    return genes, "all genes of the gene file (no biotype given)"


def chrom_key(c):
    """Natural order of sequence names (chr2 before chr10)."""
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", str(c))]


def barcode_intervals(hap):
    """Barcode intervals with coordinates, numbered in genome order (number = column of the barcode)."""
    rows = [{"interval": str(n), "chrom": str(c), "start": int(s), "end": int(e)} for n, c, s, e in zip(hap.ids, hap.chrom, hap.pos, hap.end)]
    I = pd.DataFrame(rows, columns=["interval", "chrom", "start", "end"])
    ok = (I.end > I.start).to_numpy() if len(I) else np.zeros(0, bool)
    A = I[ok].copy(); B = I[~ok].copy()
    if len(A):
        A = A.iloc[sorted(range(len(A)), key=lambda i: (chrom_key(A.chrom.iloc[i]), int(A.start.iloc[i])))]
    I = pd.concat([A, B], ignore_index=True)
    I["number"] = np.arange(1, len(I) + 1)
    I["has_coordinates"] = list(np.r_[np.ones(len(A), bool), np.zeros(len(B), bool)])
    return I


def tag_hits(genes, cds, intervals, tags):
    """Tag variants that hit a gene of their interval: one row per call x gene. `hit` says whether
    coding sequence is hit, or only the span of the gene could be tested (no CDS in the gene file)."""
    cols = ["call_id", "interval", "interval_number", "variant_id", "svtype", "chrom", "start", "end", "length_bp", "tagged_arrangement", "difference", "gene", "name", "hit"]
    if tags is None or not len(tags) or not len(genes):
        return pd.DataFrame(columns=cols)
    has_cds = len(cds) > 0
    T = (cds if has_cds else genes)[["gene", "chrom", "start", "end"]]
    name_of = dict(zip(genes.gene, genes.name)); ok_gene = set(genes.gene)
    by = {c: g.sort_values("start") for c, g in T.groupby("chrom")}
    rows = []
    for iv in intervals[intervals.has_coordinates].itertuples():
        tg = tags[tags.interval == iv.interval].sort_values(["start", "end"])
        g = by.get(iv.chrom)
        if g is None or not len(tg):
            continue
        gs, ge, gn = g.start.to_numpy(), g.end.to_numpy(), g.gene.to_numpy()
        k = 0
        for t in tg.itertuples():
            m = (ge >= t.start) & (gs <= t.end)                   # closed comparison, as the original
            hit = [x for x in dict.fromkeys(gn[m]) if x in ok_gene]
            if not hit:
                continue
            k += 1
            for x in hit:
                rows.append({"call_id": f"{iv.number}.{k}", "interval": iv.interval, "interval_number": iv.number, "variant_id": t.variant_id, "svtype": t.svtype, "chrom": t.chrom,
                             "start": int(t.start), "end": int(t.end), "length_bp": int(t.length_bp), "tagged_arrangement": t.tagged_arrangement, "difference": float(t.difference),
                             "gene": x, "name": name_of.get(x, ""), "hit": "coding sequence" if has_cds else "gene span (no CDS in the gene file)"})
    return pd.DataFrame(rows, columns=cols)


def interval_genes(genes, intervals, hits=None):
    """Genes per barcode interval: every gene that overlaps the interval, with its position
    relative to it (inside / spans the start / spans the end / spans both boundaries)."""
    cols = ["gene", "name", "chrom", "start", "end", "strand", "biotype", "interval", "interval_number", "position", "n_tag_variants_hitting"]
    nh = hits.groupby(["interval", "gene"]).variant_id.nunique().to_dict() if hits is not None and len(hits) else {}
    parts = []
    for iv in intervals[intervals.has_coordinates].itertuples():
        g = genes[(genes.chrom == iv.chrom) & (genes.end >= iv.start) & (genes.start <= iv.end)].copy()
        if not len(g):
            continue
        lo, hi = g.start < iv.start, g.end > iv.end
        g["interval"] = iv.interval; g["interval_number"] = iv.number
        g["position"] = np.where(lo & hi, "spans both boundaries", np.where(lo, "spans the start", np.where(hi, "spans the end", "inside")))
        g["n_tag_variants_hitting"] = [int(nh.get((iv.interval, x), 0)) for x in g.gene]
        parts.append(g)
    return pd.concat(parts, ignore_index=True)[cols] if parts else pd.DataFrame(columns=cols)


def barcode_tables(run):
    """Tables of step 7 that need --genes: (key, table, file name, description). Empty without
    --genes or without structural-haplotype genotypes; never needed by another step."""
    if getattr(run, "genes", None) is None:
        return []
    hap, _info = run.hap_set()
    if hap is None:
        return []
    I = barcode_intervals(hap)
    G, what = drawn_genes(run.genes)
    C = getattr(run, "genes_cds", None)
    C = C[C.gene.isin(set(G.gene))] if C is not None and len(C) else pd.DataFrame(columns=CDS_COLS)
    tags = run.t.get("haplotypes", {}).get("tags")
    H = tag_hits(G, C, I, tags)
    T = interval_genes(G, I, H)
    return [("genes", T, "barcode_interval_genes.tsv", f"Genes of every barcode interval ({what}): gene, name, chrom, start (0-based), end, strand, biotype, the interval and its number in genome order, position relative to the interval (inside / spans a boundary; a gene belongs to an interval when it overlaps it), and the number of arrangement tag variants that hit it."),
            ("gene_hits", H, "barcode_interval_tag_gene_hits.tsv", "Arrangement tag variants that hit a gene of their barcode interval (coding sequence when the gene file has CDS features, else the gene span): one row per call and gene; call_id = <interval number>.<k> in genome order.")]
