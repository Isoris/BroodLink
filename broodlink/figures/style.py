"""House style of the published figures (ported from the figure-style module of the original
study: fonts, type scale, palette, millimetre placement, panel letters, numbered banners).

Everything is a constant or a small helper. Figures are laid out in millimetres from the
top-left corner, as in the original scripts, so that panel positions can be carried over.
"""
from __future__ import annotations
import colorsys

import numpy as np

MM = 1 / 25.4
A4_H = 297.0                # mm: no figure is taller than an A4 page
WIDTH_FULL = 183.0          # mm: full text width of the published figures
# ------------------------------------------------------------------ colours
NAVY = "#17252D"            # all text, axes and ticks
NOTE_GREY = "#3D4650"       # notes and statistics inside panels
INK = "#3D4650"
BACKBONE = "#8A949C"
GREY_FAM = "#DDE3EB"        # "no family" in a colour band
RED = "#b91c1c"             # reference lines and the number they carry
GREEN = "#2F7D4F"           # Gate 1 pass
G1_FILL = "#DDEEDD"
ROH_BASE = {"navy": "#283B7C", "sky": "#95BADD", "teal": "#71A29C", "sand": "#F1E492", "brick": "#83251D"}
EXTRA = {"terra": "#D08C6E", "lav": "#B8A9CC", "grey": "#A7AEB3"}
PAL = [ROH_BASE["navy"], ROH_BASE["sand"], ROH_BASE["teal"], ROH_BASE["sky"], EXTRA["terra"], EXTRA["lav"], EXTRA["grey"], ROH_BASE["brick"]]
# marker classes: one colour per class in every figure
DEL_C, DUP_C, HAP_C, SNP_C, SNP_MASK_C, OTHER_C = "#E08E45", "#E0B72E", ROH_BASE["navy"], "#A7AEB3", "#4F5B66", "#71A29C"
# relatedness bands, lowest to highest (below third degree ... identical genotypes)
BAND_COL = ["#F4F6F9", "#C9D6E3", "#8CA9C9", "#283B7C", "#17252D"]
# barcode states (class 0 homokaryotype, heterokaryotype, class 2 homokaryotype, no information)
STATE0, STATE1, STATE2, NO_INFO = "#E8A33A", "#B98A92", "#8E6BBF", "#EEF1F4"
CONFLICT = "#2F7FC1"
# genotype heat maps of deletions
REF_C, HET_C, HOM_C, NA_C = "#F1F3F5", "#E3B193", ROH_BASE["brick"], "#FFFFFF"
# type scales (pt): main figures and Extended Data / Supplementary figures
SZ_MAIN = dict(letter=11.5, title=9.0, base=7.0, label=6.5, tick=6.0, small=5.0, note=4.8, in_panel=6.5)
SZ_ED = dict(letter=9.0, title=7.5, base=6.5, label=6.5, tick=6.0, small=5.5, note=5.5, in_panel=6.0)
LW = dict(axes=0.5, tick=0.5, data=0.7, thin=0.35, ref=0.7)
MIN_FONT = 5.5


def class_colour(klass):
    """Colour of a marker class (SNP, DEL, DUP, HAPLOTYPE, anything else)."""
    return {"SNP": SNP_C, "SNP_INDEL": SNP_C, "INDEL": SNP_C, "SNP_MASKED": SNP_MASK_C, "DEL": DEL_C, "DUP": DUP_C, "HAPLOTYPE": HAP_C}.get(str(klass), OTHER_C)


def class_label(klass, plural=True):
    """Short label of a marker class as printed in the figures."""
    lab = {"SNP": "biSNP", "SNP_INDEL": "SNP + INDEL", "INDEL": "INDEL", "DEL": "DEL", "DUP": "DUP", "HAPLOTYPE": "structural haplotype", "INV": "INV",
           "INS": "INS", "BND": "BND"}.get(str(klass), str(klass))
    if str(klass) == "SNP_MASKED":
        return "biSNPs, SV masked"
    return lab + ("s" if plural and not lab.endswith("INDEL") else "")


def ros_seq():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("ros_seq", ["#FFFFFF", "#F1E492", "#D08C6E", "#83251D"])


def nomination_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("hs", ["#FFFFFF", "#FBE3B8", "#E8A04A", "#B85A14", "#6B2A0E"])


def mute(h, sat=0.58, val_to=0.80):
    """Family hue kept, saturation reduced (same function as the original figure scripts)."""
    from matplotlib.colors import hsv_to_rgb, rgb_to_hsv, to_hex, to_rgb
    hh, ss, vv = rgb_to_hsv(to_rgb(h))
    vv = vv * 0.55 + val_to * 0.45 if vv > val_to else vv * 0.85 + 0.12
    return to_hex(hsv_to_rgb((hh, ss * sat, min(vv, 0.92))))


def family_colours(families):
    """One muted colour per family, deterministic from the order given (golden-angle hues, three
    lightness levels). The original figures read a stored colour file of one cohort; BroodLink
    derives the colours from the family list of the run."""
    from matplotlib.colors import to_hex
    out = {}
    for i, f in enumerate(families):
        h = (0.61 + i * 0.6180339887) % 1.0
        r, g, b = colorsys.hsv_to_rgb(h, (0.70, 0.55, 0.85)[i % 3], (0.72, 0.88, 0.60)[(i // 3) % 3])
        out[f] = mute(to_hex((r, g, b)))
    return out


def dark_text(colour, target_luminance=0.34):
    from matplotlib.colors import to_rgb
    r, g, b = to_rgb(colour); l = 0.299 * r + 0.587 * g + 0.114 * b; f = min(1.0, target_luminance / l) if l > 0 else 1.0
    return (r * f, g * f, b * f)


def mixw(c, f):
    """Tint (f > 0) or deepen (f < 0) a class colour."""
    from matplotlib.colors import to_rgb
    r = np.array(to_rgb(c)); return tuple(r + (1 - r) * f if f > 0 else r * (1 + f))


def house_style(plt, sz=None, thin=False):
    """Arial 7 pt (DejaVu Sans when Arial is not installed), navy text and axes, 0.5 pt axes, no
    top / right spines, no legend frame, text kept as text in PDF. Returns the type scale."""
    import matplotlib
    from matplotlib import font_manager
    sz = dict(SZ_MAIN if sz is None else sz)
    av = {f.name for f in font_manager.fontManager.ttflist}
    face = next((f for f in ("Arial", "TeX Gyre Heros", "Liberation Sans", "DejaVu Sans") if f in av), "DejaVu Sans")
    matplotlib.rcdefaults()
    lw = 0.35 if thin else LW["axes"]
    matplotlib.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42, "font.family": face, "font.size": sz["base"], "text.color": NAVY, "axes.labelcolor": NAVY,
        "axes.edgecolor": NAVY, "xtick.color": NAVY, "ytick.color": NAVY, "axes.linewidth": lw, "xtick.major.width": lw, "ytick.major.width": lw,
        "xtick.major.size": 1.6 if thin else 2.4, "ytick.major.size": 1.6 if thin else 2.4, "axes.grid": False, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "figure.facecolor": "white", "savefig.facecolor": "white", "figure.dpi": 100,
        "hatch.linewidth": 0.3, "mathtext.fontset": "custom", "mathtext.rm": face, "mathtext.it": face + ":italic", "mathtext.bf": face + ":bold",
        "mathtext.bfit": face + ":italic:bold"})
    if thin:
        matplotlib.rcParams.update({"xtick.major.pad": 1, "ytick.major.pad": 1})
    return sz


class Canvas:
    """A figure of exact size in millimetres. `ax(left, top, w, h)` places an axes from the top-left
    corner; `title`, `banner`, `text`, `swatch`, `line` draw in the same coordinates."""

    def __init__(self, plt, width_mm, height_mm, sz):
        self.plt, self.W, self.H, self.sz = plt, float(width_mm), float(height_mm), sz
        self.fig = plt.figure(figsize=(self.W * MM, self.H * MM))
        self._ov = None

    def ax(self, l, t, w, h):
        return self.fig.add_axes([l / self.W, 1 - (t + h) / self.H, w / self.W, h / self.H])

    @property
    def overlay(self):
        """Full-figure axes in millimetres (y down) for strips, swatches and marks outside panels."""
        if self._ov is None:
            self._ov = self.fig.add_axes([0, 0, 1, 1]); self._ov.set_xlim(0, self.W); self._ov.set_ylim(self.H, 0); self._ov.axis("off")
        return self._ov

    def text(self, x, y, s, fs=None, bold=False, color=NAVY, ha="left", va="center", **kw):
        return self.fig.text(x / self.W, 1 - y / self.H, s, fontsize=fs or self.sz["small"], fontweight="bold" if bold else "normal", color=color, ha=ha, va=va, **kw)

    def title(self, letter, text, top, left, dx=5.2):
        if letter:
            self.fig.text(left / self.W, 1 - top / self.H, letter, fontsize=self.sz["letter"], fontweight="bold", va="center")
        self.fig.text((left + (dx if letter else 0)) / self.W, 1 - top / self.H, text, fontsize=self.sz["title"], fontweight="bold", va="center")

    def banner(self, top, num, text, fs=8):
        from matplotlib.lines import Line2D
        self.fig.add_artist(Line2D([0, 1], [1 - top / self.H] * 2, color=NAVY, lw=.5, transform=self.fig.transFigure))
        self.fig.text(.5 / self.W, 1 - (top + 1.0) / self.H, num, fontsize=fs, fontweight="bold", color=NAVY, va="top")
        self.fig.text(5 / self.W, 1 - (top + 1.0) / self.H, text, fontsize=fs, fontweight="bold", color=NAVY, va="top")

    def swatch(self, x, y, col, lab, w=2.6, h=2.0, ec="none", fs=None, hatch=None, alpha=None):
        from matplotlib.patches import Rectangle
        self.fig.add_artist(Rectangle((x / self.W, 1 - (y + h / 2) / self.H), w / self.W, h / self.H, transform=self.fig.transFigure, fc=col, ec=ec, lw=.3, hatch=hatch, alpha=alpha))
        self.fig.text((x + w + .7) / self.W, 1 - y / self.H, lab, fontsize=fs or self.sz["small"], va="center", color=NAVY)

    def line(self, x0, x1, y, col, lw=.9, ls="-"):
        from matplotlib.lines import Line2D
        self.fig.add_artist(Line2D([x0 / self.W, x1 / self.W], [1 - y / self.H] * 2, color=col, lw=lw, ls=ls, transform=self.fig.transFigure))

    def text_width(self, s, fs, bold=False):
        """Width in mm of a string at a font size (for keys packed left to right)."""
        t = self.fig.text(0, 0, s, fontsize=fs, fontweight="bold" if bold else "normal")
        w = t.get_window_extent(self.fig.canvas.get_renderer()).width / self.fig.dpi * 25.4
        t.remove(); return w

    def set_height(self, height_mm):
        """Shrink or grow the page after the content is known is not possible with absolute axes;
        figures compute their height first. Kept for clarity: asserts the A4 limit."""
        assert height_mm <= A4_H + 1e-6, f"figure taller than A4 ({height_mm:.0f} mm)"


def canvas(plt, width_mm, height_mm, sz=None, thin=False):
    sz = house_style(plt, sz, thin)
    return Canvas(plt, width_mm, min(float(height_mm), A4_H), sz)


def frame(ax, lw=.35):
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(lw)


def fam_key(f):
    """Sort key of unit names: families (F12) by number, then unassigned individuals (S_...)."""
    s = str(f); d = "".join(ch for ch in s if ch.isdigit())
    return (1 if s.startswith("S_") else 0, int(d) if d else 0, s)


def fmt_p(p, n_null=None):
    """P value as printed in panels: exact to 4 decimals, or '< 10^-4'."""
    if p < 2e-4:
        return "$\\mathit{P}$ < 10$^{-4}$"
    return f"$\\mathit{{P}}$ = {p:.4f}"


class Ledger:
    """Every number printed in a ported panel, with the table and column it came from.
    `L(panel, item, value, source)` records the entry and returns the value, so it can be used
    inside the expression that prints it. make_figures writes the entries to figures/LEDGER.tsv."""

    def __init__(self):
        self.rows, self.figure = [], ""

    def __call__(self, panel, item, value, source):
        self.rows.append({"figure": self.figure, "panel": panel, "printed_item": item, "value": value if isinstance(value, str) else (f"{value:.6g}" if isinstance(value, float) else str(value)), "source": source})
        return value


LED = Ledger()
