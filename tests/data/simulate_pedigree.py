#!/usr/bin/env python
"""Simulate the small synthetic cohort shipped with BroodLink (no real data involved).

    python tests/data/simulate_pedigree.py [--out tests/data] [--seed 7]

One generation of offspring from unrelated founders (the parents are not sampled):
  * full-sib families of varied size, some sharing ONE parent (candidate half-sib links),
    including a family linked through both of its parents and one parent with three families;
  * unassigned individuals (founders sampled directly);
  * one individual genotyped twice ("pair of samples with identical genotypes");
  * one low-coverage sample;
  * SNP genotype likelihoods at low coverage (Beagle format), DEL / DUP genotypes (VCF),
    large structural haplotypes (inversion-like intervals without recombination, with
    arrangement-associated SNPs and tag deletions), a callable mask, depth and sample tables,
    and a relatedness table.
Everything is written with fixed names; the truth tables are in truth/.
The relatedness table holds the realised identity-by-descent coefficients of the simulated
pedigree plus Gaussian noise; it stands in for the output of a likelihood-based estimator such
as NGSrelate and is not estimated from the simulated reads.
"""
from __future__ import annotations
import argparse
import gzip
import os

import numpy as np
import pandas as pd

CHROM_LEN = 7_500_000
N_CHROM = 8
SNPS_PER_CHROM = 750
N_DEL, N_DUP = 900, 260
XO_PER_MB = 0.2
INVERSIONS = [("chr1", 2_000_000, 2_800_000, 0.40), ("chr2", 4_500_000, 5_400_000, 0.35),
              ("chr3", 1_000_000, 1_900_000, 0.45), ("chr5", 5_000_000, 5_800_000, 0.30)]
# family name -> (parent 1, parent 2, number of offspring)
FAMILIES = {"TF01": ("P1", "Q1", 8), "TF02": ("P1", "Q2", 6),
            "TF03": ("P2", "Q3", 7), "TF04": ("P2", "Q4", 5),
            "TF05": ("P3", "P4", 6), "TF06": ("P3", "Q5", 5), "TF07": ("P4", "Q6", 4),
            "TF08": ("P5", "Q7", 5), "TF09": ("P5", "Q8", 4), "TF10": ("P5", "Q9", 4),
            "TF11": ("R1", "R2", 6), "TF12": ("R3", "R4", 4), "TF13": ("R5", "R6", 3),
            "TF14": ("R7", "R8", 2), "TF15": ("R9", "R10", 2)}
N_UNASSIGNED = 8
MEAN_DEPTH = 6.0
LOW_DEPTH = 1.2
SEQ_ERR = 0.01


def gz_write(path, text):
    with gzip.GzipFile(path, "wb", mtime=0) as fh:           # mtime=0: byte-identical files for a fixed seed
        fh.write(text.encode())


def write_genes(out, seed=7, genes_per_chrom=40):
    """genes.gff3: simulated gene models (gene, mRNA, exon, CDS) over the simulated chromosomes,
    for the optional input --genes. Uses its own random stream and only the constants of this
    file, so it can be called alone (`--genes-only`) without changing any other file."""
    rng = np.random.default_rng([seed, 20261114])
    inv = {}
    for c, s, e, _f in INVERSIONS:
        inv.setdefault(c, []).append((s, e))
    lines = ["##gff-version 3", "#!source BroodLink simulate_pedigree.py (synthetic gene models, no real annotation)"]
    k = 0
    for ci in range(N_CHROM):
        c = f"chr{ci + 1}"; lines.append(f"##sequence-region {c} 1 {CHROM_LEN}")
        pitch = CHROM_LEN // genes_per_chrom
        for j in range(genes_per_chrom):
            dense = any(s <= j * pitch < e for s, e in inv.get(c, []))      # one extra gene per slot inside the structural haplotypes
            for sub in range(2 if dense else 1):
                span = int(rng.integers(5_000, 40_000))
                lo = j * pitch + sub * (pitch // 2) + int(rng.integers(2_000, pitch // 2 - span - 2_000))
                hi = min(lo + span, CHROM_LEN - 1_000)
                k += 1
                strand = "+" if rng.random() < 0.5 else "-"
                coding = rng.random() < 0.9
                gid = f"gene-SIM{k:05d}"
                name = f"sim{k}{'ab'[int(rng.integers(2))]}" if rng.random() < 0.35 else f"LOC{900000 + k}"
                bt = "protein_coding" if coding else "lncRNA"
                lines.append(f"{c}\tsimulated\tgene\t{lo + 1}\t{hi}\t.\t{strand}\t.\tID={gid};Name={name};gene_biotype={bt}")
                tid = f"rna-SIM{k:05d}"
                lines.append(f"{c}\tsimulated\t{'mRNA' if coding else 'lnc_RNA'}\t{lo + 1}\t{hi}\t.\t{strand}\t.\tID={tid};Parent={gid};Name={name}")
                n_ex = int(rng.integers(2, 10))
                cuts = np.sort(rng.choice(np.arange(lo + 200, hi - 200), size=2 * n_ex - 2, replace=False))
                bounds = [lo] + [int(x) for x in cuts] + [hi]
                exons = [(bounds[2 * i], bounds[2 * i + 1]) for i in range(n_ex)]
                for i, (s, e) in enumerate(exons, 1):
                    lines.append(f"{c}\tsimulated\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tID=exon-SIM{k:05d}-{i};Parent={tid}")
                if coding:
                    for i, (s, e) in enumerate(exons, 1):
                        lines.append(f"{c}\tsimulated\tCDS\t{s + 1}\t{e}\t.\t{strand}\t0\tID=cds-SIM{k:05d};Parent={tid}")
    with open(os.path.join(out, "genes.gff3"), "w", newline="\n") as fo:
        fo.write("\n".join(lines) + "\n")
    return k


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--genes-only", action="store_true", help="write only genes.gff3 (simulated gene models); no other file is touched")
    a = ap.parse_args()
    if a.genes_only:
        print(f"wrote {write_genes(a.out, a.seed)} simulated genes to {os.path.join(a.out, 'genes.gff3')}"); return
    rng = np.random.default_rng(a.seed)
    out = a.out
    os.makedirs(os.path.join(out, "truth"), exist_ok=True)
    chroms = [f"chr{i + 1}" for i in range(N_CHROM)]

    # ---------------------------------------------------------------- callable mask
    mask = []
    for c in chroms:
        p = 0
        while p < CHROM_LEN:
            L = int(rng.integers(120_000, 400_000)); e = min(p + L, CHROM_LEN)
            mask.append((c, p, e)); p = e + int(rng.integers(30_000, 110_000))
    mask = pd.DataFrame(mask, columns=["chrom", "start", "end"])

    def callable_at(c, pos):
        m = mask[mask.chrom == c]
        i = np.searchsorted(m.start.to_numpy(), pos - 1, side="right") - 1
        return (i >= 0) & (pos - 1 < m.end.to_numpy()[np.clip(i, 0, None)])

    def draw_positions(c, n, share_callable):
        pos = set()
        while len(pos) < n:
            x = rng.integers(1000, CHROM_LEN - 6000, size=n * 3)
            ok = callable_at(c, x)
            keep = np.where(ok, rng.random(len(x)) < share_callable, rng.random(len(x)) < (1 - share_callable) * 3)
            for v in x[keep]:
                if len(pos) < n:
                    pos.add(int(v))
        return np.array(sorted(pos))

    # ---------------------------------------------------------------- loci
    loci = []                                                   # (type, chrom, pos, end, id)
    for c in chroms:
        for p in draw_positions(c, SNPS_PER_CHROM, 0.97):
            loci.append(("SNP", c, p, p, f"{c}_{p}"))
    inv_of = {}
    for k, n, share in (("DEL", N_DEL, 0.80), ("DUP", N_DUP, 0.80)):
        per = n // N_CHROM
        for c in chroms:
            for j, p in enumerate(draw_positions(c, per, share)):
                L = int(np.exp(rng.uniform(np.log(300), np.log(8000)))) if k == "DEL" else int(np.exp(rng.uniform(np.log(1000), np.log(20000))))
                loci.append((k, c, p, p + L, f"{k}_{c}_{p}"))
    for ii, (c, s, e, _f) in enumerate(INVERSIONS):             # tag deletions inside each interval
        for p in np.sort(rng.integers(s + 20_000, e - 20_000, size=20)):
            loci.append(("DEL", c, int(p), int(p) + int(rng.integers(400, 3000)), f"DEL_{c}_{int(p)}"))
    L = pd.DataFrame(loci, columns=["type", "chrom", "pos", "end", "id"]).drop_duplicates("id")
    L = L.sort_values(["chrom", "pos"]).reset_index(drop=True)
    L["inv"] = -1
    for ii, (c, s, e, _f) in enumerate(INVERSIONS):
        L.loc[(L.chrom == c) & (L.pos >= s) & (L.pos <= e), "inv"] = ii
    nL = len(L)

    # ---------------------------------------------------------------- founders
    founders = sorted({p for v in FAMILIES.values() for p in v[:2]}) + [f"U{i + 1}" for i in range(N_UNASSIGNED)]
    nF = len(founders); nH = 2 * nF
    arr = np.array([rng.random(nH) < f for (_c, _s, _e, f) in INVERSIONS]).astype(int)       # inversions x founder haplotypes
    for ii in range(len(INVERSIONS)):                            # each arrangement present on at least 8 founder haplotypes
        while arr[ii].sum() < 8 or arr[ii].sum() > nH - 8:
            arr[ii] = rng.random(nH) < INVERSIONS[ii][3]
    A = np.zeros((nH, nL), np.int8)                              # founder haplotype alleles
    typ = L.type.to_numpy(); inv = L.inv.to_numpy()
    for j in range(nL):
        if inv[j] >= 0 and ((typ[j] == "SNP" and rng.random() < 0.6) or (typ[j] == "DEL" and rng.random() < 0.75)):
            lo, hi = (rng.uniform(0.0, 0.05), rng.uniform(0.92, 1.0))
            if typ[j] == "DEL":
                lo, hi = 0.0, rng.uniform(0.9, 1.0)
            p0, p1 = (lo, hi) if rng.random() < 0.5 else (hi, lo)
            pr = np.where(arr[inv[j]] == 1, p1, p0)
            A[:, j] = rng.random(nH) < pr
        else:
            if typ[j] == "SNP":
                p = rng.uniform(0.08, 0.92)
            else:
                p = min(max(rng.beta(0.5, 1.6), 1.5 / nH), 0.9)
            A[:, j] = rng.random(nH) < p
    fh = {f: (2 * i, 2 * i + 1) for i, f in enumerate(founders)}

    # ---------------------------------------------------------------- meiosis
    cidx = {c: np.flatnonzero(L.chrom.to_numpy() == c) for c in chroms}
    cpos = {c: L.pos.to_numpy()[cidx[c]] for c in chroms}
    inv_iv = {c: [(s, e) for (cc, s, e, _f) in INVERSIONS if cc == c] for c in chroms}

    def gamete(h1, h2):
        g = np.empty(nL, np.int16)
        for c in chroms:
            n = rng.poisson(XO_PER_MB * CHROM_LEN / 1e6)
            xo = np.sort(rng.integers(1, CHROM_LEN, size=n))
            xo = [x for x in xo if not any(s <= x <= e for s, e in inv_iv[c])]      # no crossover inside an interval
            phase = (np.searchsorted(np.array(xo, dtype=np.int64), cpos[c]) + int(rng.integers(2))) % 2
            g[cidx[c]] = np.where(phase == 0, h1[cidx[c]], h2[cidx[c]])
        return g

    H = {}                                                       # individual -> (hap ids 1, hap ids 2)
    for f in founders:
        H[f] = (np.full(nL, fh[f][0], np.int16), np.full(nL, fh[f][1], np.int16))
    truth = []
    k = 0
    for fam, (p1, p2, n) in FAMILIES.items():
        for _ in range(n):
            k += 1; name = f"kid{k}"
            H[name] = (gamete(*H[p1]), gamete(*H[p2]))
            truth.append({"name": name, "family": fam, "parent_1": p1, "parent_2": p2, "note": ""})
    for i in range(N_UNASSIGNED):
        truth.append({"name": f"U{i + 1}", "family": "unassigned", "parent_1": "", "parent_2": "", "note": "founder sampled directly"})
    T = pd.DataFrame(truth)
    twin_src = T[T.family == "TF11"].name.iloc[0]
    T = pd.concat([T, pd.DataFrame([{"name": "twin", "family": "TF11", "parent_1": "R1", "parent_2": "R2",
                                     "note": "same individual genotyped twice (second library)"}])], ignore_index=True)
    H["twin"] = H[twin_src]
    order = rng.permutation(len(T))
    T = T.iloc[order].reset_index(drop=True)
    T.insert(0, "individual", [f"IND{i + 1:03d}" for i in range(len(T))])
    src_id = T.set_index("name").individual[twin_src]
    T.loc[T.name == "twin", "note"] = f"same individual as {src_id}, genotyped twice"
    T.loc[T.name == twin_src, "note"] = f"also genotyped as {T.set_index('name').individual['twin']}"
    low = T[T.family == "unassigned"].individual.iloc[0]
    T.loc[T.individual == low, "note"] = "founder sampled directly; low-coverage sample"
    ids = T.individual.tolist(); N = len(ids)
    hap1 = np.array([H[n][0] for n in T.name]); hap2 = np.array([H[n][1] for n in T.name])       # N x nL
    G = (A[hap1, np.arange(nL)[None, :]] + A[hap2, np.arange(nL)[None, :]]).astype(np.int8).T    # nL x N

    # ---------------------------------------------------------------- realised relatedness
    snp = np.flatnonzero(typ == "SNP")
    a1, a2 = hap1[:, snp], hap2[:, snp]
    rows = []
    for i in range(N):
        for j in range(i + 1, N):
            ibd = np.maximum((a1[i] == a1[j]).astype(int) + (a2[i] == a2[j]), (a1[i] == a2[j]).astype(int) + (a2[i] == a1[j]))
            k2, k1 = float((ibd == 2).mean()), float((ibd == 1).mean())
            k2n = min(max(k2 + rng.normal(0, 0.01), 0), 1); k1n = min(max(k1 + rng.normal(0, 0.015), 0), 1 - k2n)
            rows.append((ids[i], ids[j], k1n / 4 + k2n / 2, 1 - k1n - k2n, k1n, k2n))
    R = pd.DataFrame(rows, columns=["id1", "id2", "theta", "k0", "k1", "k2"])
    R.to_csv(os.path.join(out, "relatedness.tsv"), sep="\t", index=False, float_format="%.5f", lineterminator="\n")

    # ---------------------------------------------------------------- SNP likelihoods (Beagle)
    depth = MEAN_DEPTH * np.exp(rng.normal(0, 0.15, N))
    depth[ids.index(low)] = LOW_DEPTH
    Gs = G[snp]
    nreads = rng.poisson(depth[None, :], size=Gs.shape)
    pa = np.where(Gs == 0, SEQ_ERR, np.where(Gs == 1, 0.5, 1 - SEQ_ERR))
    nalt = rng.binomial(nreads, pa)
    lik = np.stack([(SEQ_ERR ** nalt) * ((1 - SEQ_ERR) ** (nreads - nalt)), 0.5 ** nreads * np.ones_like(pa),
                    ((1 - SEQ_ERR) ** nalt) * (SEQ_ERR ** (nreads - nalt))], 2)
    lik = lik / lik.sum(2, keepdims=True)
    lines = ["marker\tallele1\tallele2\t" + "\t".join(f"{s}\t{s}\t{s}" for s in ids)]
    txt = np.char.mod("%.4f", lik.reshape(len(snp), N * 3))
    for r, j in enumerate(snp):
        lines.append(f"{L.chrom[j]}_{L.pos[j]}\t0\t1\t" + "\t".join(txt[r]))
    gz_write(os.path.join(out, "snps.beagle.gz"), "\n".join(lines) + "\n")

    # ---------------------------------------------------------------- SV genotypes (VCF)
    sv = np.flatnonzero(typ != "SNP")
    Gv = G[sv].copy()
    err = rng.random(Gv.shape) < 0.01
    Gv[err] = rng.integers(0, 3, size=int(err.sum()))
    miss = rng.random(Gv.shape) < 0.02
    miss[:, ids.index(low)] |= rng.random(len(sv)) < 0.25
    keep = (np.where(miss, 0, Gv) > 0).any(1)
    code = np.array(["0/0", "0/1", "1/1"])
    v = ["##fileformat=VCFv4.2", "##source=BroodLink simulate_pedigree.py (synthetic data)",
         '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
         '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    v += [f"##contig=<ID={c},length={CHROM_LEN}>" for c in chroms]
    v.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids))
    for r, j in enumerate(sv):
        if not keep[r]:
            continue
        gt = np.where(miss[r], "./.", code[Gv[r]])
        v.append(f"{L.chrom[j]}\t{L.pos[j]}\t{L.id[j]}\tN\t<{typ[j]}>\t.\tPASS\tSVTYPE={typ[j]};END={L.end[j]}\tGT\t" + "\t".join(gt))
    gz_write(os.path.join(out, "sv.vcf.gz"), "\n".join(v) + "\n")

    # ---------------------------------------------------------------- structural haplotypes
    inv_names = [f"{c}:{s}-{e}" for (c, s, e, _f) in INVERSIONS]
    mid = [int(np.flatnonzero((L.inv.to_numpy() == ii))[0]) for ii in range(len(INVERSIONS))]
    K = np.array([arr[ii][hap1[:, mid[ii]]] + arr[ii][hap2[:, mid[ii]]] for ii in range(len(INVERSIONS))])     # inversions x N
    with open(os.path.join(out, "inversions.bed"), "w", newline="\n") as fo:
        for ii, (c, s, e, _f) in enumerate(INVERSIONS):
            fo.write(f"{c}\t{s}\t{e}\t{inv_names[ii]}\n")
    Kobs = K.copy()
    flip = rng.random(K.shape) < 0.01
    Kobs[flip] = (Kobs[flip] + 1) % 3
    hp = pd.DataFrame(Kobs.T, index=ids, columns=inv_names); hp.index.name = "individual"
    hp.to_csv(os.path.join(out, "haplotypes.tsv"), sep="\t", lineterminator="\n")

    # ---------------------------------------------------------------- mask, genome, samples, depth
    mask.to_csv(os.path.join(out, "callable.bed"), sep="\t", index=False, header=False, lineterminator="\n")
    with open(os.path.join(out, "genome.fai"), "w", newline="\n") as fo:
        for c in chroms:
            fo.write(f"{c}\t{CHROM_LEN}\t0\t60\t61\n")
    S = pd.DataFrame({"id": ids, "sex": rng.choice(["F", "M"], N), "batch": rng.choice(["B1", "B2"], N), "group": "cohort"})
    S.to_csv(os.path.join(out, "samples.tsv"), sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"sample": ids, "depth_mask": np.round(depth, 2), "depth_assembly": np.round(depth * rng.uniform(1.5, 1.9, N), 2)}
                 ).to_csv(os.path.join(out, "depth.tsv"), sep="\t", index=False, lineterminator="\n")

    # ---------------------------------------------------------------- truth
    tdir = os.path.join(out, "truth")
    T[["individual", "family", "parent_1", "parent_2", "note"]].to_csv(os.path.join(tdir, "truth_families.tsv"), sep="\t", index=False, lineterminator="\n")
    links = []
    fams = list(FAMILIES)
    for i, fa in enumerate(fams):
        for fb in fams[i + 1:]:
            sh = set(FAMILIES[fa][:2]) & set(FAMILIES[fb][:2])
            if sh:
                links.append({"family_a": fa, "family_b": fb, "shared_parent": sorted(sh)[0]})
    pd.DataFrame(links).to_csv(os.path.join(tdir, "truth_links.tsv"), sep="\t", index=False, lineterminator="\n")
    par = sorted({p for v_ in FAMILIES.values() for p in v_[:2]})
    pk = [{"parent": p, **{inv_names[ii]: int(arr[ii][fh[p][0]] + arr[ii][fh[p][1]]) for ii in range(len(INVERSIONS))}} for p in par]
    pd.DataFrame(pk).to_csv(os.path.join(tdir, "truth_parent_haplotypes.tsv"), sep="\t", index=False, lineterminator="\n")
    tk = pd.DataFrame(K.T, index=ids, columns=inv_names); tk.index.name = "individual"
    tk.to_csv(os.path.join(tdir, "truth_karyotypes.tsv"), sep="\t", lineterminator="\n")
    pd.DataFrame({"sample": [low, T.set_index("name").individual["twin"], src_id],
                  "what": ["low-coverage sample", "identical genotypes (pair)", "identical genotypes (pair)"]}
                 ).to_csv(os.path.join(tdir, "truth_qc.tsv"), sep="\t", index=False, lineterminator="\n")
    n_genes = write_genes(out, a.seed)
    print(f"wrote {N} individuals, {len(snp)} SNPs, {int(keep.sum())} SV calls, {len(INVERSIONS)} structural haplotypes, {n_genes} genes to {out}")


if __name__ == "__main__":
    main()
