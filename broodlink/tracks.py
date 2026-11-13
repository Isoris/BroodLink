"""Genome-browser tracks. GFF3 is 1-based, closed; BED is 0-based, half-open.

Interval tables of step 06 follow the BED convention of --inversions (0-based start, end
exclusive); variant tables follow the VCF convention (1-based POS, END inclusive)."""
from __future__ import annotations
import os
from urllib.parse import quote

import pandas as pd

from . import __version__, families


def _attr(d):
    return ";".join(f"{k}={quote(str(v), safe=' ,.:()/+-_')}" for k, v in d.items() if v is not None and str(v) != "" and str(v) != "nan")


def write_tracks(run):
    d = os.path.join(run.out, "tracks"); os.makedirs(d, exist_ok=True)
    made = []
    S = run.t.get("haplotypes", {}).get("summary")
    head = f"##gff-version 3\n# BroodLink {__version__} | config {run.cfg_hash}\n# {families.HEADER}\n"
    if S is not None and len(S):
        S = S[S.tier != "not genotyped"] if "tier" in S else S
        with open(os.path.join(d, "structural_haplotypes.gff3"), "w", encoding="utf-8", newline="\n") as g, \
                open(os.path.join(d, "structural_haplotypes.bed"), "w", encoding="utf-8", newline="\n") as b:
            g.write(head + "# genotype labels 0 / 2 are arbitrary cluster labels; intervals from the scan are candidates, detected automatically, not validated\n")
            b.write('track name="structural_haplotypes" description="BroodLink structural haplotypes (candidates unless supplied)" itemRgb="On"\n')
            for r in S.itertuples():
                scan = str(r.source) == "scan"
                at = {"ID": r.interval, "Name": r.interval, "interval_source": r.source, "tier": r.tier, "genotyped_from": getattr(r, "genotyped_from", ""),
                      "n_markers": getattr(r, "n_markers", ""), "n_0": getattr(r, "n_0", ""), "n_1": getattr(r, "n_1", ""), "n_2": getattr(r, "n_2", ""),
                      "minor_arrangement_frequency": f"{getattr(r, 'minor_arrangement_frequency', float('nan')):.4f}",
                      "Note": "candidate, detected automatically, not validated" if scan else ""}
                g.write(f"{r.chrom}\tBroodLink\tstructural_haplotype\t{int(r.start) + 1}\t{int(r.end)}\t.\t.\t.\t{_attr(at)}\n")
                rgb = "150,150,150" if scan else ("0,114,178" if r.tier == "high_confidence" else "230,159,0")
                b.write(f"{r.chrom}\t{int(r.start)}\t{int(r.end)}\t{r.interval}\t0\t.\t{int(r.start)}\t{int(r.end)}\t{rgb}\n")
        made += [("structural_haplotypes.gff3", "intervals genotyped by local PCA; attributes: source, tier, markers, cluster sizes, minor-arrangement frequency"),
                 ("structural_haplotypes.bed", "the same intervals; blue = high confidence, orange = candidate (supplied interval), grey = found by the scan")]
    T = run.t.get("haplotypes", {}).get("tags")
    if T is not None and len(T):
        with open(os.path.join(d, "arrangement_tag_variants.gff3"), "w", encoding="utf-8", newline="\n") as g, \
                open(os.path.join(d, "arrangement_tag_variants.bed"), "w", encoding="utf-8", newline="\n") as b:
            g.write(head)
            b.write('track name="arrangement_tag_variants" description="BroodLink DEL / DUP tagging one arrangement class"\n')
            for i, r in enumerate(T.itertuples()):
                at = {"ID": f"tag{i + 1}", "Name": r.variant_id, "interval": r.interval, "svtype": r.svtype, "tagged_arrangement": r.tagged_arrangement,
                      "difference": f"{r.difference:.3f}", "fully_private": str(bool(r.fully_private)).lower()}
                g.write(f"{r.chrom}\tBroodLink\t{'deletion' if r.svtype == 'DEL' else 'copy_number_variation'}\t{int(r.start)}\t{int(r.end)}\t{r.difference:.3f}\t.\t.\t{_attr(at)}\n")
                b.write(f"{r.chrom}\t{int(r.start) - 1}\t{int(r.end)}\t{r.variant_id}\t{int(round(1000 * r.difference))}\t.\n")
        made += [("arrangement_tag_variants.gff3", "DEL / DUP whose carrier frequency differs between the two homokaryotype classes; score = difference"),
                 ("arrangement_tag_variants.bed", "the same variants; score = 1000 x difference")]
    P = run.t.get("family_effect", {}).get("per_marker")
    if P is not None and len(P):
        X = P[(P.eta2_obs > P.null_max) & (P.marker_class != "HAPLOTYPE")]
        if len(X):
            with open(os.path.join(d, "family_informative_variants.bed"), "w", encoding="utf-8", newline="\n") as b:
                b.write('track name="family_informative_variants" description="BroodLink variants with eta squared above every permutation"\n')
                for r in X.sort_values(["chrom", "pos"]).itertuples():
                    b.write(f"{r.chrom}\t{int(r.pos) - 1}\t{int(r.pos)}\t{r.marker_class}:{r.id}\t{int(round(1000 * r.eta2_obs))}\t.\n")
            made.append(("family_informative_variants.bed", "variants whose eta squared (variance explained by family) exceeds every permutation; position only (start of the variant); score = 1000 x eta squared"))
    with open(os.path.join(d, "README.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"BroodLink {__version__} - tracks\nGFF3: 1-based, closed coordinates. BED: 0-based, half-open coordinates.\n{families.HEADER}\n\n")
        for f, desc in made:
            fh.write(f"{f}\n  {desc}\n\n")
        if not made:
            fh.write("No track was written: no structural haplotype was genotyped and no family-informative variant was found.\n")
    return made
