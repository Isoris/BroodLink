"""Figures of a run (PNG + PDF), drawn from the tables the run produced.

`FIGURES` is the registry: (name, function, step it needs). A function takes (run, plt) and
returns None (nothing to draw: the input of the figure is absent), or
(figure, title, description[, note]), or a list of (suffix, figure, title, description, note)
when one figure is drawn several times (once per chain of linked families). File stems are
`figNN_<name>` with NN = position in the registry, so numbers do not depend on what was drawn.

Modules: `style` (house style of the published figures), `generic` (views without a published
counterpart), and one module per published figure family.
"""
from __future__ import annotations
import os

import pandas as pd

from ..io import detail, log, write_tsv
from . import generic, kinship, shared_parent, style

# (name, function, step whose tables it needs, "house" = style set by the figure itself)
FIGURES = [
    ("overview", kinship.fig_overview, "families"),
    ("qc_panels", generic.f_qc, "qc"),
    ("structural_haplotypes_pc1", generic.f_hap, "haplotypes"),
    ("structural_haplotype_scan", generic.f_scan, "haplotypes"),
    ("relatedness_theta_k2", generic.f_rel, "relatedness"),
    ("family_sizes_closure", generic.f_fam, "families"),
    ("arrangement_tag_variants", generic.f_tags, "haplotypes"),
    ("barcode", generic.f_barcode, "barcode"),
    ("shared_parent_nomination_audit_network", shared_parent.fig_shared_parent, "nomination"),
    ("nomination_matrix_all_units", shared_parent.fig_matrix_all_units, "nomination"),
    ("sv_genotypes_and_kinship", kinship.fig_sv_kinship, "concordance"),
    ("cross_marker_concordance", kinship.fig_concordance, "concordance"),
]
FIRST_NUMBER = 0


def stems():
    return [f"fig{FIRST_NUMBER + k:02d}_{name}" for k, (name, _f, _s) in enumerate(FIGURES)]


def _save(fig, d, stem, dpi):
    fig.canvas.draw()
    try:
        fig.set_layout_engine("none")                              # lay out once, then save twice
    except Exception:
        pass
    fig.savefig(os.path.join(d, stem + ".png"), dpi=dpi); fig.savefig(os.path.join(d, stem + ".pdf"), dpi=dpi)


def make_figures(run, only=None):
    """Draw every figure of the registry (or those whose name contains `only`); never stops a run."""
    d = os.path.join(run.out, "figures"); os.makedirs(d, exist_ok=True)
    dpi = run.cfg["figures"]["dpi"]; rows = []
    style.LED.rows = []
    log(f"figures: drawing up to {len(FIGURES)} figures (PNG {dpi} dpi + PDF) into {d}")
    for k, ((name, fn, step), stem) in enumerate(zip(FIGURES, stems())):
        if only and only not in name:
            continue
        detail(f"[{k + 1}/{len(FIGURES)}] {stem}")
        style.LED.figure = stem
        plt = generic._plt()                                       # resets the rc settings; ported figures then set the house style
        try:
            if step not in run.t or (step == "nomination" and "pairs" not in run.t["nomination"]):
                res = None
            else:
                res = fn(run, plt)
        except Exception as exc:                                   # a figure never stops a run
            log(f"WARNING: figure {stem} failed ({type(exc).__name__}: {exc})")
            rows.append({"file": stem, "title": "", "status": f"failed: {type(exc).__name__}: {exc}", "description": "", "note": ""}); plt.close("all"); continue
        if res is None:
            rows.append({"file": stem, "title": "", "status": "not drawn: the input of this figure is absent or the step produced nothing to draw", "description": "", "note": ""}); continue
        if isinstance(res, str):                                   # a reason for not drawing
            rows.append({"file": stem, "title": "", "status": "not drawn: " + res, "description": "", "note": ""}); continue
        items = res if isinstance(res, list) else [("",) + tuple(res)]
        for it in items:
            suffix, fig, title, desc = it[0], it[1], it[2], it[3]; note = it[4] if len(it) > 4 else ""
            st = stem + (f"_{suffix}" if suffix else "")
            try:
                _save(fig, d, st, dpi)
            except Exception as exc:                               # e.g. a write error: report it, keep the run
                log(f"WARNING: figure {st} could not be saved ({type(exc).__name__}: {exc})")
                rows.append({"file": st, "title": title, "status": f"failed: {type(exc).__name__}: {exc}", "description": desc, "note": note}); plt.close("all"); continue
            plt.close(fig)
            rows.append({"file": st, "title": title, "status": "drawn (PNG + PDF)", "description": desc, "note": note or ""})
        plt.close("all")
    M = pd.DataFrame(rows, columns=["file", "title", "status", "description", "note"])
    L = pd.DataFrame(style.LED.rows, columns=["figure", "panel", "printed_item", "value", "source"])
    if not only:
        write_tsv(L, os.path.join(d, "LEDGER.tsv"), run.header("Ledger of the ported figures: every number printed in a panel, with the table (and column) of the run it was computed from."))
    write_tsv(M, os.path.join(d, "MANIFEST.tsv"), run.header("Figures of the run: file stem, title, status, description, and a note on panels left out because their input is absent."))
    return M
