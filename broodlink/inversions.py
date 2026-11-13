"""Step 6a: large structural haplotypes (inversion-like intervals) by local PCA.

Genotyping and tiers: extracted from inv41_unrelated_rerun.py (standardised dosage, SVD, 1-D
k-means with k = 3 on PC1, five criteria). Deletion-based classes: extracted from
arrangement_classes (relatedness_inversion_pipeline.py). Deletion scan: extracted from
step2_inversions.py. Sub-clusters and their deletion support: extracted from
compute_subclusters_v1.py and build_extdata_haplotype_subclusters_v1.py.
SNP scan: REIMPLEMENTED FROM METHODS (lostruct-style; the original scan script was not archived).

Genotype labels: 0 and 2 are the two homokaryotype clusters, 1 the middle cluster. 0 is the more
frequent homokaryotype. The labels carry no ancestral meaning ("STD"/"INV" are names only).
Everything found by the scan is "candidate, detected automatically, not validated".
"""
from __future__ import annotations
import time

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm

from .family_effect import eta_vec
from .io import allele_freq, log

SCAN_LABEL = "candidate, detected automatically, not validated"


# ------------------------------------------------------------------ building blocks
def kmeans1d(x, k=3, it=100):
    c = np.quantile(x, [0.1, 0.5, 0.9])
    lab = np.zeros(len(x), int)
    for _ in range(it):
        lab = np.argmin(np.abs(x[:, None] - c[None, :]), axis=1)
        new = np.array([x[lab == j].mean() if np.any(lab == j) else c[j] for j in range(k)])
        if np.allclose(new, c):
            break
        c = new
    order = np.argsort(c); remap = np.empty(k, int); remap[order] = np.arange(k)
    return remap[lab], np.sort(c)


def sep(a, b):
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sd = np.sqrt((a.var(ddof=1) * (len(a) - 1) + b.var(ddof=1) * (len(b) - 1)) / (len(a) + len(b) - 2))
    return float(abs(b.mean() - a.mean()) / sd) if sd > 0 else np.inf


def local_pca(X, fit):
    """X: sites x individuals (NaN missing). PCA fitted on the individuals `fit`, all projected.
    Returns pc1, pc2 (all individuals), share of variance of PC1, number of sites used."""
    Xf = X[:, fit]
    mu = np.nanmean(Xf, axis=1, keepdims=True); sd = np.nanstd(Xf, axis=1, keepdims=True)
    ok = (sd[:, 0] > 0) & np.isfinite(sd[:, 0])
    if ok.sum() < 3:
        return None
    Xc = np.nan_to_num((X[ok] - mu[ok]) / sd[ok])
    U, S, _ = np.linalg.svd(Xc[:, fit], full_matrices=False)
    pc = U[:, :2].T @ Xc
    return pc[0], (pc[1] if pc.shape[0] > 1 else np.zeros(X.shape[1])), float(S[0] ** 2 / (S ** 2).sum()), int(ok.sum())


def assign(pc1, fit, cfg):
    """k-means on the fitted individuals; every individual goes to the nearest centre;
    individuals with a class posterior below the threshold get a missing genotype."""
    lab_fit, cent = kmeans1d(pc1[fit])
    forced = np.argmin(np.abs(pc1[:, None] - cent[None, :]), axis=1)
    sds, w = [], []
    for j in range(3):
        v = pc1[fit][lab_fit == j]
        sds.append(v.std(ddof=1) if len(v) > 2 else np.nan); w.append(max(len(v), 1))
    pooled = np.sqrt(np.nansum([(n - 1) * s ** 2 for n, s in zip(w, sds) if not np.isnan(s)]) / max(sum(n - 1 for n, s in zip(w, sds) if not np.isnan(s)), 1))
    pooled = pooled if pooled > 0 else 1e-9
    dens = np.stack([w[j] * norm.pdf(pc1, cent[j], pooled) for j in range(3)], 1)
    post = dens.max(1) / np.maximum(dens.sum(1), 1e-300)
    lab = np.where(post >= cfg["inversions"]["ambiguous_posterior"], forced, -1)
    return forced, lab, cent, post


def _mw(a, b):
    return float(mannwhitneyu(a, b, alternative="greater").pvalue) if len(a) > 1 and len(b) > 1 else np.nan


def genotype_matrix(X, het, samples, unrel, fam_of, cfg, seed):
    """Karyotype one interval from a sites x individuals dosage matrix.
    Returns (per-individual frame, statistics dict) or (None, reason)."""
    c = cfg["inversions"]; t = c["tier"]
    n = X.shape[1]; allidx = np.arange(n)
    six = {s: i for i, s in enumerate(samples)}
    uidx = np.array([six[s] for s in unrel if s in six], int)
    use_unrel = len(uidx) >= c["min_unrelated"]
    full = local_pca(X, allidx)
    if full is None:
        return None, "too few variable sites"
    pc1_full, pc2_full, var_full, nsite = full
    forced_full, lab_full, _, _ = assign(pc1_full, allidx, cfg)
    if use_unrel:
        r = local_pca(X, uidx)
        if r is None:
            use_unrel = False
    if use_unrel:
        pc1, pc2, var1, nsite = r
        if np.corrcoef(pc1, pc1_full)[0, 1] < 0:
            pc1, pc2 = -pc1, pc2
        fit = uidx
    else:
        pc1, pc2, var1, fit = pc1_full, pc2_full, var_full, allidx
    forced, lab, cent, post = assign(pc1, fit, cfg)
    # orientation: 0 = the more frequent homokaryotype (labels carry no ancestral meaning)
    if (forced == 2).sum() > (forced == 0).sum():
        pc1 = -pc1; forced = 2 - forced; lab = np.where(lab >= 0, 2 - lab, -1)
    if np.corrcoef(pc1, pc1_full)[0, 1] < 0:
        pc1_full = -pc1_full; forced_full = 2 - forced_full
    ff = forced[fit]
    g = [pc1[fit][ff == j] for j in range(3)]; h = [het[fit][ff == j] for j in range(3)]
    s01, s12 = sep(g[0], g[1]), sep(g[1], g[2])
    p_het = np.nanmax([_mw(h[1], h[0]), _mw(h[1], h[2])]) if len(h[1]) > 1 else np.nan
    conc = float(np.mean(forced[fit] == forced_full[fit]))
    nn = [int((ff == j).sum()) for j in range(3)]
    minor = min(nn[0], nn[2]) * 2 + nn[1]
    # criterion 5: family eta squared of the arrangement above its permutation null
    cols = [i for i, s in enumerate(samples) if fam_of.get(s, "unassigned") != "unassigned" and lab[i] >= 0]
    eta, eta_q, eta_ok = np.nan, np.nan, False
    if len(cols) >= 8 and len({fam_of[samples[i]] for i in cols}) >= 2:
        Gk = lab[cols][None, :].astype(np.int8); fl = np.array([fam_of[samples[i]] for i in cols])
        eta = float(eta_vec(Gk, fl, 2)[0])
        rng = np.random.default_rng(seed)
        nul = np.array([eta_vec(Gk, rng.permutation(fl), 2)[0] for _ in range(t["eta2_permutations"])])
        eta_q = float(np.nanquantile(nul, t["eta2_null_quantile"])) if np.isfinite(nul).any() else np.nan
        eta_ok = bool(eta > eta_q) if np.isfinite(eta) and np.isfinite(eta_q) else False
    crit = {"separation": bool(np.isfinite(s01) and np.isfinite(s12) and s01 >= t["separation_sd"] and s12 >= t["separation_sd"]),
            "heterozygote_excess": bool(np.isfinite(p_het) and p_het < t["het_p"]),
            "concordance": bool(conc >= t["concordance"]) if use_unrel else None,
            "minor_copies": bool(minor >= t["min_minor_copies"]),
            "family_eta2": eta_ok}
    passed = [v for v in crit.values() if v is not None]
    called = lab >= 0
    stats = {"n_sites": nsite, "pc1_variance_share": var1, "clusters_fitted_on": f"unrelated subset ({len(uidx)} individuals), others projected" if use_unrel else "all individuals (no unrelated subset of the minimum size)",
             "n_fit": int(len(fit)), "n_0": int((lab == 0).sum()), "n_1": int((lab == 1).sum()), "n_2": int((lab == 2).sum()), "n_missing_ambiguous": int((~called).sum()),
             "minor_arrangement_frequency": float(np.where(called, lab, 0).sum() / (2.0 * max(called.sum(), 1))),
             "separation_0_1": s01, "separation_1_2": s12, "het_0": float(np.median(h[0])) if len(h[0]) else np.nan,
             "het_1": float(np.median(h[1])) if len(h[1]) else np.nan, "het_2": float(np.median(h[2])) if len(h[2]) else np.nan,
             "het_excess_p": p_het, "concordance_subset_vs_all": conc if use_unrel else np.nan, "minor_copies_in_fit": int(minor),
             "family_eta2": eta, "family_eta2_null_q": eta_q,
             **{f"pass_{k}": ("not evaluated" if v is None else ("pass" if v else "fail")) for k, v in crit.items()},
             "tier": "high_confidence" if all(passed) and len(passed) == 5 else "candidate"}
    ind = pd.DataFrame({"individual": samples, "genotype": lab, "genotype_forced": forced, "posterior": post, "pc1": pc1, "pc2": pc2, "het_inside": het,
                        "in_fit_subset": np.isin(np.arange(n), fit)})
    return ind, stats


def snp_block(ms, chrom, start, end):
    m = (ms.chrom == chrom) & (ms.pos >= start) & (ms.pos <= end)
    X = ms.dosage()[m]
    G = ms.G[m]
    with np.errstate(invalid="ignore", divide="ignore"):
        if ms.D is not None:
            het = np.nanmean(np.where(np.isnan(X), np.nan, ((X > 0.5) & (X < 1.5)).astype(float)), axis=0)
        else:
            het = np.where((G >= 0).sum(0) > 0, (G == 1).sum(0) / np.maximum((G >= 0).sum(0), 1), np.nan)
    return X, np.nan_to_num(het)


def del_block(ms, chrom, start, end, cfg):
    """Deletions of the interval co-segregating with the block (arrangement_classes)."""
    c = cfg["inversions"]
    af = allele_freq(ms.G)
    m = (ms.chrom == chrom) & (ms.pos >= start) & (ms.pos <= end) & (af > c["del_af_min"]) & (af < c["del_af_max"])
    if m.sum() < c["min_dels"]:
        return None, None, int(m.sum())
    blk = ms.G[m]; X = blk.astype(float); X[blk < 0] = np.nan
    mu = np.nanmean(X, axis=1, keepdims=True)
    Xi = np.where(np.isnan(X), mu, X)
    Xs = (Xi - Xi.mean(1, keepdims=True)) / (Xi.std(1, keepdims=True) + 1e-9)
    C = np.corrcoef(Xs); mac = (np.abs(C).sum(1) - 1) / max(len(C) - 1, 1)
    keep = mac > c["del_link_cut"]
    if keep.sum() < 5:
        keep = mac > np.percentile(mac, 50)
    Xk = X[keep]; Gk = blk[keep]
    het = np.where((Gk >= 0).sum(0) > 0, (Gk == 1).sum(0) / np.maximum((Gk >= 0).sum(0), 1), 0.0)
    return Xk, het, int(keep.sum())


def genotype_intervals(intervals, snp, dels, samples, unrel, fam_of, cfg, seed):
    """intervals: frame (interval, chrom, start, end, ...). Returns (summary, per-individual long table,
    genotype matrix dict interval -> array aligned to `samples`)."""
    c = cfg["inversions"]
    rows, inds, geno = [], [], {}
    for k, r in enumerate(intervals.itertuples()):
        base = {"interval": r.interval, "chrom": r.chrom, "start": int(r.start), "end": int(r.end), "span_bp": int(r.end - r.start),
                "source": getattr(r, "source", "supplied"), "tier_supplied": getattr(r, "tier_supplied", "")}
        res = {}
        if snp is not None:
            ss = snp.subset_samples(samples)
            X, het = snp_block(ss, r.chrom, r.start, r.end)
            if X.shape[0] >= c["min_snps"]:
                ind, st = genotype_matrix(X, het, ss.samples, unrel, fam_of, cfg, seed + k)
                if ind is not None:
                    res["SNP"] = (ind, st, ss.samples, X.shape[0])
        if dels is not None:
            dd = dels.subset_samples(samples)
            Xd, hetd, nd = del_block(dd, r.chrom, r.start, r.end, cfg)
            base["n_deletions_used"] = nd
            if Xd is not None:
                ind, st = genotype_matrix(Xd, hetd, dd.samples, unrel, fam_of, cfg, seed + k)
                if ind is not None:
                    res["DEL"] = (ind, st, dd.samples, nd)
        if not res:
            rows.append({**base, "genotyped_from": "not genotyped", "note": f"fewer than {c['min_snps']} SNPs" + (f" and fewer than {c['min_dels']} deletions" if dels is not None else "") + " in the interval",
                         "tier": "not genotyped"})
            continue
        prim = "SNP" if "SNP" in res else "DEL"
        ind, st, smp, nmark = res[prim]
        row = {**base, "genotyped_from": "SNPs" if prim == "SNP" else "from deletions only", "n_markers": nmark, **st, "note": ""}
        if "SNP" in res and "DEL" in res:
            a = res["SNP"][0].set_index("individual"); b = res["DEL"][0].set_index("individual").reindex(a.index)
            both = (a.genotype >= 0) & (b.genotype >= 0)
            if both.sum():
                agree = float((a.genotype[both] == b.genotype[both]).mean()); flip = float((a.genotype[both] == 2 - b.genotype[both]).mean())
                row["deletion_call_concordance"] = max(agree, flip)
                row["deletion_call_orientation"] = "same" if agree >= flip else "flipped (labels are arbitrary)"
                row["deletion_tier"] = res["DEL"][1]["tier"]
            dl = res["DEL"][0].copy(); dl.insert(0, "interval", r.interval); dl.insert(1, "markers", "deletions"); inds.append(dl)
        ind = ind.copy(); ind.insert(0, "interval", r.interval); ind.insert(1, "markers", "SNPs" if prim == "SNP" else "deletions")
        inds.append(ind)
        gmap = dict(zip(ind.individual, ind.genotype))
        geno[r.interval] = np.array([gmap.get(s, -1) for s in samples], np.int8)
        rows.append(row)
    S = pd.DataFrame(rows)
    I = pd.concat(inds, ignore_index=True) if inds else pd.DataFrame(columns=["interval", "markers", "individual", "genotype"])
    return S, I, geno


# ------------------------------------------------------------------ scans
def _mds(D2, k):
    n = D2.shape[0]
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ D2 @ J
    w, v = np.linalg.eigh(B)
    o = np.argsort(w)[::-1][:k]
    return v[:, o] * np.sqrt(np.maximum(w[o], 0))


def scan_snps(ms, cfg):
    """lostruct-style scan, REIMPLEMENTED FROM METHODS: windows of n SNPs, k local PCs, distance between
    window PCAs, multidimensional scaling, outlier windows at |z| >= threshold, merged into regions."""
    c = cfg["inversions"]["snp_scan"]; t0 = time.time()
    step = max(1, int(np.ceil(ms.n / c["thin_to"])))
    W, k = c["window_snps"], c["n_pcs"]
    wins, notes, cands = [], [], []
    chroms = list(dict.fromkeys(ms.chrom.tolist()))
    for ch in chroms:
        idx = np.flatnonzero(ms.chrom == ch)[::step]
        X = ms.subset_markers(idx).dosage()
        pos = ms.pos[idx]
        var = np.nanvar(X, axis=1); ok = var > 0
        X, pos = X[ok], pos[ok]
        nw = X.shape[0] // W
        if nw < c["min_windows_per_chrom"]:
            notes.append(f"{ch}: skipped by the SNP scan ({X.shape[0]} variable SNPs after thinning give {nw} windows; minimum {c['min_windows_per_chrom']})")
            continue
        V = np.zeros((nw, X.shape[1], k)); A = np.zeros((nw, k)); p1 = np.zeros(nw); hr = np.zeros(nw)
        for w in range(nw):
            Xw = X[w * W:(w + 1) * W]
            Z = np.nan_to_num(Xw - np.nanmean(Xw, axis=1, keepdims=True))
            C = Z.T @ Z / Z.shape[0]
            ev, evec = np.linalg.eigh(C)
            V[w] = evec[:, -k:][:, ::-1]; lam = np.maximum(ev[-k:][::-1], 0)
            A[w] = lam / max(lam.sum(), 1e-12); p1[w] = ev[-1] / max(ev.sum(), 1e-12)
            pc1 = evec[:, -1]; het = np.nanmean(((Xw > 0.5) & (Xw < 1.5)).astype(float), axis=0)
            q = np.percentile(pc1, [33, 67]); mid = (pc1 >= q[0]) & (pc1 < q[1])
            hr[w] = het[mid].mean() / max(het[~mid].mean(), 1e-9) if mid.any() and (~mid).any() else np.nan
        dots = np.einsum("ink,jnl->ijkl", V, V) ** 2
        cross = np.einsum("ik,jl,ijkl->ij", A, A, dots)
        selfn = (A ** 2).sum(1)
        D2 = np.maximum(selfn[:, None] + selfn[None, :] - 2 * cross, 0)
        M = _mds(D2, min(c["mds_coords"], nw - 1))
        z = (M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1)
        out = (np.abs(z[:, :c["test_coords"]]) >= c["z_outlier"]).any(1)
        for w in range(nw):
            wins.append({"chrom": ch, "window": w, "start": int(pos[w * W]), "end": int(pos[(w + 1) * W - 1]), "n_snps": W,
                         "mds1": M[w, 0], "mds2": M[w, 1] if M.shape[1] > 1 else np.nan, "z1": z[w, 0], "z2": z[w, 1] if z.shape[1] > 1 else np.nan,
                         "pc1_variance_share": p1[w], "het_ratio_middle_vs_outer": hr[w], "outlier": bool(out[w])})
        # merge outlier windows of the same sign on the same coordinate
        for co in range(min(c["test_coords"], z.shape[1])):
            for sign in (1, -1):
                hit = np.flatnonzero(sign * z[:, co] >= c["z_outlier"])
                if not len(hit):
                    continue
                runs = [[hit[0], hit[0]]]
                for w in hit[1:]:
                    if w - runs[-1][1] <= c["merge_gap_windows"] + 1:
                        runs[-1][1] = w
                    else:
                        runs.append([w, w])
                for a, b in runs:
                    nhit = int(((hit >= a) & (hit <= b)).sum())
                    s, e = int(pos[a * W]), int(pos[(b + 1) * W - 1])
                    if nhit >= c["min_windows"] and e - s >= c["min_size_bp"]:
                        cands.append({"chrom": ch, "start": s, "end": e, "n_outlier_windows": nhit, "mds_coordinate": co + 1,
                                      "max_abs_z": float(np.abs(z[a:b + 1, co]).max()), "mean_pc1_variance_share": float(p1[a:b + 1].mean()),
                                      "mean_het_ratio": float(np.nanmean(hr[a:b + 1])), "scan": "SNP local PCA"})
    C = _merge_candidates(pd.DataFrame(cands), "SNP")
    return C, pd.DataFrame(wins), notes, time.time() - t0, step


def scan_dels(ms, cfg):
    """Deletion scan (step2_inversions.py): fixed windows, LD = mean |r| of deletion dosages,
    windows at or above a quantile merged; a region needs a strong core window."""
    c = cfg["inversions"]["del_scan"]; ci = cfg["inversions"]; t0 = time.time()
    af = allele_freq(ms.G); poly = (af > ci["del_af_min"]) & (af < ci["del_af_max"])
    rows = []
    win = (ms.pos // c["window_bp"]) * c["window_bp"]
    for ch in dict.fromkeys(ms.chrom.tolist()):
        mc = (ms.chrom == ch) & poly
        for w in np.unique(win[mc]):
            pm = np.flatnonzero(mc & (win == w))
            if len(pm) < c["min_dels"]:
                continue
            D = ms.G[pm].astype(float); D[ms.G[pm] < 0] = np.nan
            C = pd.DataFrame(D.T).corr().to_numpy(); iu = np.triu_indices_from(C, 1)
            rows.append({"chrom": ch, "start": int(w), "end": int(w + c["window_bp"]), "n_deletions": len(pm), "LD_mean_abs_r": float(np.nanmean(np.abs(C[iu])))})
    T = pd.DataFrame(rows)
    notes, cands = [], []
    if len(T) < 10:
        notes.append(f"deletion scan: only {len(T)} windows with at least {c['min_dels']} deletions; scan skipped")
        return pd.DataFrame(), T, notes, time.time() - t0
    thr, core = T.LD_mean_abs_r.quantile(c["hit_quantile"]), T.LD_mean_abs_r.quantile(c["core_quantile"])
    T["outlier"] = T.LD_mean_abs_r >= thr
    for ch, g in T.groupby("chrom", sort=False):
        g = g.sort_values("start").reset_index(drop=True); run = None; runs = []
        for i, hval in enumerate(g.outlier):
            if hval:
                run = [i, i] if run is None else [run[0], i]
            elif run is not None and (g.start[i] - g.start[run[1]]) > c["merge_gap_bp"]:
                runs.append(run); run = None
        if run is not None:
            runs.append(run)
        for a, b in runs:
            seg = g.iloc[a:b + 1]
            if seg.LD_mean_abs_r.max() < core:
                continue
            m = (ms.chrom == ch) & (ms.pos >= seg.start.min()) & (ms.pos < seg.end.max()) & poly
            cands.append({"chrom": ch, "start": int(ms.pos[m].min()), "end": int(ms.pos[m].max()), "n_outlier_windows": int(seg.outlier.sum()),
                          "max_LD": float(seg.LD_mean_abs_r.max()), "n_deletions": int(seg.n_deletions.sum()), "scan": "deletion LD"})
    return _merge_candidates(pd.DataFrame(cands), "DEL"), T, notes, time.time() - t0


def _merge_candidates(C, tag):
    if not len(C):
        return pd.DataFrame(columns=["interval", "chrom", "start", "end", "source", "scan"])
    out = []
    for ch, g in C.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        cur = None
        for r in g.to_dict("records"):
            if cur is not None and r["start"] <= cur["end"]:
                cur["end"] = max(cur["end"], r["end"]); cur["n_outlier_windows"] = max(cur["n_outlier_windows"], r["n_outlier_windows"])
            else:
                if cur is not None:
                    out.append(cur)
                cur = dict(r)
        out.append(cur)
    C = pd.DataFrame(out)
    C.insert(0, "interval", [f"scan{tag}_{c}:{s}-{e}" for c, s, e in zip(C.chrom, C.start, C.end)])
    C["source"] = "scan"; C["label"] = SCAN_LABEL
    return C


def overlap_table(supplied, scanned):
    """Supplied intervals recovered by the scan, not recovered, and extra candidates."""
    rows = []
    used = set()
    for r in supplied.itertuples():
        best, bo = None, 0
        for s in scanned.itertuples():
            if s.chrom == r.chrom:
                ov = min(r.end, s.end) - max(r.start, s.start)
                if ov > bo:
                    best, bo = s, ov
        if best is not None:
            used.add(best.interval)
            rows.append({"category": "supplied interval recovered by the scan", "supplied_interval": r.interval, "scan_interval": best.interval, "chrom": r.chrom,
                         "overlap_bp": int(bo), "share_of_supplied": bo / max(r.end - r.start, 1), "share_of_scan": bo / max(best.end - best.start, 1)})
        else:
            rows.append({"category": "supplied interval not recovered", "supplied_interval": r.interval, "scan_interval": "", "chrom": r.chrom, "overlap_bp": 0,
                         "share_of_supplied": 0.0, "share_of_scan": np.nan})
    for s in scanned.itertuples():
        if s.interval not in used:
            rows.append({"category": "extra candidate interval found by the scan", "supplied_interval": "", "scan_interval": s.interval, "chrom": s.chrom,
                         "overlap_bp": 0, "share_of_supplied": np.nan, "share_of_scan": 0.0})
    return pd.DataFrame(rows, columns=["category", "supplied_interval", "scan_interval", "chrom", "overlap_bp", "share_of_supplied", "share_of_scan"])


# ------------------------------------------------------------------ arrangement-tag variants
def tag_variants(intervals, geno, samples, sv_sets, cfg):
    """Deletions / duplications whose carrier frequency differs between the two homokaryotypes."""
    c = cfg["inversions"]["tags"]; rows = []
    for ms in sv_sets:
        mm = ms.subset_samples(samples); six = {s: i for i, s in enumerate(samples)}
        col = np.array([six[s] for s in mm.samples])
        for r in intervals.itertuples():
            if r.interval not in geno:
                continue
            g = geno[r.interval][col]
            m = np.flatnonzero((mm.chrom == r.chrom) & (mm.pos >= r.start) & (mm.pos <= r.end))
            if not len(m) or (g == 0).sum() < c["min_called_per_class"] or (g == 2).sum() < c["min_called_per_class"]:
                continue
            G = mm.G[m]
            c0 = G[:, g == 0]; c2 = G[:, g == 2]
            n0 = (c0 >= 0).sum(1); n2 = (c2 >= 0).sum(1)
            f0 = (c0 > 0).sum(1) / np.maximum(n0, 1); f2 = (c2 > 0).sum(1) / np.maximum(n2, 1)
            ok = (n0 >= c["min_called_per_class"]) & (n2 >= c["min_called_per_class"])
            d = f2 - f0
            for i in np.flatnonzero(ok & (np.abs(d) >= c["min_delta"])):
                rows.append({"interval": r.interval, "chrom": r.chrom, "start": int(mm.pos[m[i]]), "end": int(max(mm.end[m[i]], mm.pos[m[i]])), "variant_id": mm.ids[m[i]],
                             "svtype": mm.klass, "length_bp": int(max(mm.end[m[i]] - mm.pos[m[i]], 0)), "carrier_frequency_class_0": float(f0[i]),
                             "carrier_frequency_class_2": float(f2[i]), "difference": float(abs(d[i])), "tagged_arrangement": "class 2" if d[i] > 0 else "class 0",
                             "n_called_class_0": int(n0[i]), "n_called_class_2": int(n2[i]),
                             "fully_private": bool((f0[i] == 0 and f2[i] == 1) or (f0[i] == 1 and f2[i] == 0))})
    cols = ["interval", "chrom", "start", "end", "variant_id", "svtype", "length_bp", "carrier_frequency_class_0", "carrier_frequency_class_2",
            "difference", "tagged_arrangement", "n_called_class_0", "n_called_class_2", "fully_private"]
    return pd.DataFrame(rows, columns=cols).sort_values(["chrom", "start", "end"]).reset_index(drop=True) if rows else pd.DataFrame(columns=cols)


# ------------------------------------------------------------------ sub-clusters inside arrangements
def _gmm_fit(X, k, rng, restarts):
    n, d = X.shape; best = None; floor = 1e-3 * np.var(X, axis=0).mean() + 1e-9
    for _ in range(restarts):
        c = [X[rng.integers(n)]]
        for _k in range(1, k):
            dd = np.min([((X - ci) ** 2).sum(1) for ci in c], axis=0)
            c.append(X[rng.choice(n, p=dd / dd.sum())] if dd.sum() > 0 else X[rng.integers(n)])
        mu = np.array(c); lab = np.argmin(((X[:, None] - mu[None]) ** 2).sum(2), 1)
        pi = np.array([max((lab == j).mean(), 1e-3) for j in range(k)]); pi /= pi.sum()
        cov = np.array([np.cov(X[lab == j].T) + floor * np.eye(d) if (lab == j).sum() > 2 else np.cov(X.T) + floor * np.eye(d) for j in range(k)])
        ll_old = -np.inf
        for _ in range(200):
            logp = np.zeros((n, k))
            for j in range(k):
                dif = X - mu[j]; inv = np.linalg.inv(cov[j]); _, ld = np.linalg.slogdet(cov[j])
                logp[:, j] = np.log(pi[j]) - .5 * (np.einsum("ij,jk,ik->i", dif, inv, dif) + ld + d * np.log(2 * np.pi))
            m = logp.max(1, keepdims=True); ll = float((m[:, 0] + np.log(np.exp(logp - m).sum(1))).sum())
            r = np.exp(logp - m); r /= r.sum(1, keepdims=True)
            if abs(ll - ll_old) < 1e-7:
                break
            ll_old = ll; nk = r.sum(0) + 1e-9; pi = nk / n; mu = (r.T @ X) / nk[:, None]
            cov = np.array([((r[:, j, None] * (X - mu[j])).T @ (X - mu[j])) / nk[j] + floor * np.eye(d) for j in range(k)])
        if best is None or ll > best[0]:
            best = (ll, r.argmax(1), mu)
    ll, lab, mu = best; p = k * (d + d * (d + 1) / 2) + (k - 1)
    return -2 * ll + p * np.log(n), lab


def _cluster_class(X, sc):
    n = len(X)
    if n < sc["min_class"]:
        return np.zeros(n, int)
    best = None
    for k in range(1, sc["k_max"] + 1):
        if n < sc["min_component"] * k:
            break
        b, lab = _gmm_fit(X, k, np.random.default_rng(11 + k), sc["restarts"])
        if np.bincount(lab, minlength=k).min() < sc["min_component"] and k > 1:
            continue
        if best is None or b < best[0]:
            best = (b, lab)
    lab = best[1].copy()
    while True:                                           # merge sub-clusters closer than merge_sd
        ids = sorted(set(lab))
        if len(ids) == 1:
            break
        worst = None
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                A, B = X[lab == ids[i]], X[lab == ids[j]]; dv = A.mean(0) - B.mean(0); w = dv / (np.linalg.norm(dv) + 1e-12)
                sd = np.sqrt((np.var(A @ w) + np.var(B @ w)) / 2 + 1e-12); sv = abs(dv @ w) / sd
                if worst is None or sv < worst[0]:
                    worst = (sv, ids[i], ids[j])
        if worst[0] >= sc["merge_sd"]:
            break
        lab[lab == worst[2]] = worst[1]
    ids = sorted(set(lab), key=lambda i: X[lab == i, 0].mean()); remap = {o: i for i, o in enumerate(ids)}
    return np.array([remap[v] for v in lab])


def diagnostic_deletions(Xi, sub, sc, rng):
    """Number of deletions with |difference of mean dosage| >= delta between two sub-clusters, with a
    label-permutation null. Xi: deletions x individuals of one class (NaN missing)."""
    with np.errstate(invalid="ignore"):
        af = np.nanmean(Xi, 1) / 2
    Xi = Xi[(af > sc["del_af_min"]) & (af < sc["del_af_max"])]
    ks = sorted(set(sub))

    def nd(lbl):
        with np.errstate(invalid="ignore"):
            mu = np.array([np.nanmean(Xi[:, lbl == k], 1) for k in ks])
        best = np.zeros(Xi.shape[0])
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                best = np.fmax(best, np.abs(mu[i] - mu[j]))
        return int((best >= sc["diagnostic_delta"]).sum())
    if Xi.shape[0] == 0:
        return 0, 0.0, np.nan, 0
    obs = nd(sub); null = np.array([nd(rng.permutation(sub)) for _ in range(sc["n_permutations"])])
    return obs, float(null.mean()), float((1 + (null >= obs).sum()) / (len(null) + 1)), int(Xi.shape[0])


def subclusters(intervals, per_ind, dels, cfg, seed):
    """Sub-clusters inside each karyotype class and whether deletion dosage supports each split."""
    sc = cfg["inversions"]["subclusters"]; rows, assign_rows = [], []
    rng = np.random.default_rng(seed + 3)
    dX = None
    if dels is not None:
        dX = dels.G.astype(float); dX[dels.G < 0] = np.nan
        dix = {s: i for i, s in enumerate(dels.samples)}
    for r in intervals.itertuples():
        a = per_ind[(per_ind.interval == r.interval)]
        a = a[a.markers == a.markers.iloc[-1]] if len(a) else a
        if not len(a):
            continue
        for g in (0, 1, 2):
            cl = a[a.genotype == g]
            if len(cl) < sc["min_class"]:
                rows.append({"interval": r.interval, "class": g, "n_individuals": len(cl), "n_subclusters": 1, "note": f"class below {sc['min_class']} individuals: one cluster"})
                continue
            lab = _cluster_class(cl[["pc1", "pc2"]].to_numpy(float), sc)
            row = {"interval": r.interval, "class": g, "n_individuals": len(cl), "n_subclusters": int(lab.max() + 1), "note": ""}
            for s, v in zip(cl.individual, lab):
                assign_rows.append({"interval": r.interval, "class": g, "individual": s, "subcluster": int(v)})
            if lab.max() >= 1 and dX is not None and g in (0, 2):
                keep = [(dix[s], v) for s, v in zip(cl.individual, lab) if s in dix]
                m = (dels.chrom == r.chrom) & (dels.pos >= r.start) & (dels.pos <= r.end)
                if len(keep) >= sc["min_class"] and m.any():
                    obs, nm, p, nt = diagnostic_deletions(dX[m][:, [k for k, _ in keep]], np.array([v for _, v in keep]), sc, rng)
                    row.update({"n_deletions_tested": nt, "n_diagnostic_deletions": obs, "expected_by_permutation": nm, "p_permutation": p,
                                "split_supported_by_deletions": bool(p < 0.05)})
            rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(assign_rows, columns=["interval", "class", "individual", "subcluster"])
