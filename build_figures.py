#!/usr/bin/env python3
"""Render every figure for the BrainVar GRCh38-vs-T2T eQTL report.

Each figure is emitted twice -- <name>.light.svg and <name>.dark.svg -- from the
same code path, with only the palette swapped.  The page shows one or the other
via prefers-color-scheme, so dark mode is a selected set of steps rather than an
automatic inversion of the light one.

Colours are the validated categorical slots (validate_palette.js, all checks pass
in both modes).  Two light-mode slots sit below 3:1 contrast, so every figure
carries direct value labels and every figure is mirrored by a data table in the
page -- the relief the validator requires.
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)
DATA = json.loads(
    (Path("/mnt/ssd/lalli/nf_stage/brainvar_t2t_neurodevelopment_gain_loss")
     / "scratch/report_data.json").read_text()
)

CELLS = ["linear_grch38_dv", "graph_grch38_dv", "linear_t2t_dv", "graph_t2t_dv"]
LABEL = {
    "linear_grch38_dv": "linear · GRCh38",
    "graph_grch38_dv": "graph · GRCh38",
    "linear_t2t_dv": "linear · T2T",
    "graph_t2t_dv": "graph · T2T",
}

LIGHT = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
             grid="#e1e0d9", axis="#c3c2b7",
             s1="#2a78d6", s2="#eb6834", s3="#1baf7a", s4="#eda100")
DARK = dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
            grid="#2c2c2a", axis="#383835",
            s1="#3987e5", s2="#d95926", s3="#199e70", s4="#c98500")

THOUSANDS = FuncFormatter(lambda v, _: f"{int(v):,}")


def frame(ax, c, *, ygrid=True):
    """Recessive chrome: hairline grid, no top/right spine, muted ticks."""
    ax.set_facecolor(c["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(c["axis"])
        ax.spines[side].set_linewidth(1.0)
    if ygrid:
        ax.yaxis.grid(True, color=c["grid"], linewidth=1.0)
        ax.set_axisbelow(True)
    ax.tick_params(colors=c["muted"], labelsize=9, length=0)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(c["ink2"])


def finish(fig, axes, c, name):
    fig.patch.set_facecolor(c["surface"])
    for ax in (axes if isinstance(axes, (list, tuple)) else [axes]):
        ax.set_facecolor(c["surface"])
    fig.savefig(FIG / name, format="svg", bbox_inches="tight",
                facecolor=c["surface"], transparent=False)
    plt.close(fig)


def title(ax, c, text, sub=None, width=96):
    """Title above subtitle, with the title's pad sized to the wrapped subtitle.

    matplotlib gives the title a pad measured from the axes top, while an
    annotation grows upward from its own offset, so the two collide unless the
    pad is computed from the subtitle's line count.
    """
    pad = 9
    if sub:
        lines = textwrap.wrap(sub, width=width)
        pad = 15 + 13 * len(lines)
        ax.annotate("\n".join(lines), xy=(0, 1), xytext=(0, 7),
                    xycoords="axes fraction", textcoords="offset points",
                    color=c["ink2"], fontsize=9.5, ha="left", va="bottom",
                    linespacing=1.35)
    ax.set_title(text, color=c["ink"], fontsize=12.5, fontweight="600",
                 loc="left", pad=pad)


def legend(ax, c, **kw):
    lg = ax.legend(frameon=False, fontsize=9.5, labelcolor=c["ink2"], **kw)
    return lg


# ---------------------------------------------------------------- figure 1
def fig_missing_burden(c, mode):
    raw = DATA["nan_whole_sample"]
    arms = [a for a in CELLS if a in raw]
    pct = [100 * raw[a]["source_variant_rows_with_nan"] / raw[a]["variants"] for a in arms]
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    y = range(len(arms))
    ax.barh(list(y), pct, height=0.58, color=c["s1"], zorder=3)
    for i, (p, a) in enumerate(zip(pct, arms)):
        n = raw[a]["source_variant_rows_with_nan"]
        ax.annotate(f"{p:.1f}%   ({n:,} variants)", xy=(p, i), xytext=(8, 0),
                    textcoords="offset points", va="center", ha="left",
                    color=c["ink2"], fontsize=9.5)
    ax.set_yticks(list(y), [LABEL[a] for a in arms])
    ax.invert_yaxis()
    ax.set_xlim(0, max(pct) * 1.45)
    ax.set_xlabel("share of cis-tested variants carrying at least one missing call",
                  color=c["ink2"], fontsize=9.5)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.xaxis.grid(True, color=c["grid"], linewidth=1.0)
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)
    frame(ax, c, ygrid=False)
    title(ax, c, "A quarter of variants carried a missing call",
          "Each was silently dropped before the fix. Afterwards, NaNs reaching TensorQTL: "
          "zero, in every analysis point.")
    finish(fig, ax, c, f"missing-burden.{mode}.svg")


# ---------------------------------------------------------------- figure 2
def fig_pc_sweep(c, mode):
    curve = DATA["expression_pc_curve"]
    ks = sorted(int(k) for k in curve[CELLS[0]])
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for a, key in zip(CELLS, ("s1", "s2", "s3", "s4")):
        ys = [curve[a][str(k)] if str(k) in curve[a] else curve[a][k] for k in ks]
        ax.plot(ks, ys, color=c[key], linewidth=2.0, marker="o", markersize=4.5,
                markeredgecolor=c["surface"], markeredgewidth=1.2,
                label=LABEL[a], zorder=3)
    # No end-of-line labels: all four arms converge within ~70 eGenes at k=50 and
    # the labels overplot each other.  The accompanying table carries the values.
    ax.axvline(35, color=c["axis"], linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
    # Placed at mid-height right of the rule: the curves occupy the top of the
    # panel and the legend the bottom-right corner.
    ax.annotate("k = 35 in force", xy=(35, 0.42), xytext=(7, 0),
                xycoords=("data", "axes fraction"), textcoords="offset points",
                color=c["ink2"], fontsize=9, va="center")
    ax.set_xticks(ks)
    ax.set_xlim(min(ks) - 2, max(ks) + 2)
    ax.set_xlabel("expression principal components", color=c["ink2"], fontsize=9.5)
    ax.set_ylabel("eGenes", color=c["ink2"], fontsize=9.5)
    ax.yaxis.set_major_formatter(THOUSANDS)
    frame(ax, c)
    title(ax, c, "Expression-PC sweep on corrected genotypes",
          "Three arms peak at 45 and linear · T2T is still rising at 50, so 35 is not an "
          "argmax.")
    legend(ax, c, loc="lower right", ncol=2)
    finish(fig, ax, c, f"pc-sweep.{mode}.svg")


# ---------------------------------------------------------------- figure 3
def fig_four_arm(c, mode):
    a_int = DATA["four_arm_int_k35"]
    a_log = DATA["four_arm_logcpm_k35"]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    x = range(len(CELLS))
    w = 0.36
    b1 = [a_int[a]["egenes_q05"] for a in CELLS]
    b2 = [a_log[a]["egenes_q05"] for a in CELLS]
    ax.bar([i - w / 2 - 0.01 for i in x], b1, w, color=c["s1"], label="INT scale", zorder=3)
    ax.bar([i + w / 2 + 0.01 for i in x], b2, w, color=c["s2"], label="log2(CPM+1) scale", zorder=3)
    for i, (v1, v2) in enumerate(zip(b1, b2)):
        ax.annotate(f"{v1:,}", xy=(i - w / 2 - 0.01, v1), xytext=(0, 4),
                    textcoords="offset points", ha="center", color=c["ink2"], fontsize=9)
        ax.annotate(f"{v2:,}", xy=(i + w / 2 + 0.01, v2), xytext=(0, 4),
                    textcoords="offset points", ha="center", color=c["ink2"], fontsize=9)
    ax.set_xticks(list(x), [LABEL[a] for a in CELLS])
    ax.set_ylim(0, max(b2) * 1.18)
    ax.set_ylabel("eGenes at BH 5%", color=c["ink2"], fontsize=9.5)
    ax.yaxis.set_major_formatter(THOUSANDS)
    frame(ax, c)
    title(ax, c, "Four-arm cis-eQTL yield at k = 35",
          f"{sum(b1):,} eGenes on the INT scale, {sum(b2):,} on log-CPM, across the four arms.")
    legend(ax, c, loc="upper right", ncol=2)
    finish(fig, ax, c, f"four-arm-yield.{mode}.svg")


# ---------------------------------------------------------------- figure 4
def fig_stratified(c, mode):
    strat = DATA["stratified_maps"]
    panels = ["linear_grch38_dv", "linear_t2t_dv"]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0), sharey=True)
    for ax, arm in zip(axes, panels):
        for stratum, key in (("XX", "s1"), ("XY", "s2")):
            d = strat[f"{arm}.{stratum}"]
            ks = sorted(int(k) for k in d)
            ys = [d[str(k)] if str(k) in d else d[k] for k in ks]
            ax.plot(ks, ys, color=c[key], linewidth=2.0, marker="o", markersize=4.5,
                    markeredgecolor=c["surface"], markeredgewidth=1.2,
                    label=stratum, zorder=3)
            peak = max(range(len(ys)), key=lambda i: ys[i])
            ax.annotate(f"{ys[peak]:,}", xy=(ks[peak], ys[peak]), xytext=(0, 7),
                        textcoords="offset points", ha="center",
                        color=c["ink2"], fontsize=9)
        ax.set_xticks(ks)
        ax.set_xlabel("expression PCs", color=c["ink2"], fontsize=9.5)
        ax.set_title(LABEL[arm], color=c["ink2"], fontsize=10, loc="left", pad=6)
        frame(ax, c)
    axes[0].set_ylabel("eGenes at BH 5%", color=c["ink2"], fontsize=9.5)
    axes[0].yaxis.set_major_formatter(THOUSANDS)
    axes[0].set_ylim(0, 1600)
    legend(axes[1], c, loc="center right")
    fig.text(0.045, 1.105, "Sex-stratified maps: XY yields about 2.2x XX at every k",
             color=c["ink"], fontsize=12.5, fontweight="600", ha="left", va="top")
    fig.text(0.045, 1.035, "XX n = 92, XY n = 133. The gap tracks sample size, and the two "
             "strata peak at different k.",
             color=c["ink2"], fontsize=9.5, ha="left", va="top")
    finish(fig, list(axes), c, f"stratified-maps.{mode}.svg")


# ---------------------------------------------------------------- figure 5
def fig_interaction(c, mode):
    inter = DATA["logcpm_interaction"]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    x = range(len(CELLS))
    w = 0.26
    series = [("INT, k = 35", [0] * 4, "s1"),
              ("log-CPM, k = 10", [inter["e10"][a] for a in CELLS], "s2"),
              ("log-CPM, k = 35", [inter["e35"][a] for a in CELLS], "s3")]
    for j, (name, vals, key) in enumerate(series):
        off = (j - 1) * (w + 0.02)
        ax.bar([i + off for i in x], vals, w, color=c[key], label=name, zorder=3)
        for i, v in enumerate(vals):
            ax.annotate(f"{v}", xy=(i + off, v), xytext=(0, 4),
                        textcoords="offset points", ha="center",
                        color=c["ink2"], fontsize=9, fontweight="600" if v == 0 else "normal")
    ax.set_xticks(list(x), [LABEL[a] for a in CELLS])
    ax.set_ylim(0, 7.6)
    ax.set_ylabel("interaction eGenes at BH 5%", color=c["ink2"], fontsize=9.5)
    frame(ax, c)
    title(ax, c, "Genotype×sex interaction appears only on the log-CPM scale",
          "INT is a rank transform, so it discards the effect magnitude an interaction term "
          "estimates.")
    legend(ax, c, loc="upper right", ncol=3)
    finish(fig, ax, c, f"interaction-scale.{mode}.svg")


# ---------------------------------------------------------------- figure 6
def fig_calibration(c, mode):
    ct = DATA["sex_contrast"]
    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    x = range(len(CELLS))
    w = 0.36
    ivw = [ct[a]["ivw"]["simulated_null_rate"] for a in CELLS]
    eit = [ct[a]["either"]["simulated_null_rate"] for a in CELLS]
    ax.bar([i - w / 2 - 0.01 for i in x], ivw, w, color=c["s1"],
           label="ivw selection (valid)", zorder=3)
    ax.bar([i + w / 2 + 0.01 for i in x], eit, w, color=c["s2"],
           label="either-sex selection (invalid)", zorder=3)
    for i, (v1, v2) in enumerate(zip(ivw, eit)):
        ax.annotate(f"{v1:.4f}", xy=(i - w / 2 - 0.01, v1), xytext=(0, 4),
                    textcoords="offset points", ha="center", color=c["ink2"], fontsize=9)
        ax.annotate(f"{v2:.4f}", xy=(i + w / 2 + 0.01, v2), xytext=(0, 4),
                    textcoords="offset points", ha="center", color=c["ink2"], fontsize=9)
    ax.axhline(0.05, color=c["axis"], linewidth=1.4, linestyle=(0, (4, 3)), zorder=4)
    # Reserve a strip of empty axis to the right of the last bar group so the
    # reference-line label never overplots a bar or a per-bar value label.
    ax.set_xlim(-0.62, len(CELLS) - 1 + 0.95)
    ax.annotate("nominal 0.05", xy=(len(CELLS) - 1 + 0.42, 0.05), xytext=(0, 6),
                textcoords="offset points", color=c["ink2"], fontsize=9, ha="left")
    ax.set_xticks(list(x), [LABEL[a] for a in CELLS])
    ax.set_ylim(0, 0.21)
    ax.set_ylabel("simulated null rate", color=c["ink2"], fontsize=9.5)
    frame(ax, c)
    title(ax, c, "Only one selection rule is FDR-valid",
          "Selecting on the inverse-variance-weighted mean is orthogonal to the contrast. "
          "Selecting on either sex is not.")
    legend(ax, c, loc="upper left", ncol=2)
    finish(fig, ax, c, f"selection-calibration.{mode}.svg")


# ---------------------------------------------------------------- figure 7
def fig_direction(c, mode):
    nc = DATA["direction_null_check"]
    arms = [a for a in CELLS if a in nc]
    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    y = list(range(len(arms)))
    for i, a in enumerate(arms):
        o, n = nc[a]["observed"], nc[a]["null"]
        ax.plot([n, o], [i, i], color=c["axis"], linewidth=2.0, zorder=2,
                solid_capstyle="round")
        ax.plot(n, i, "o", color=c["s2"], markersize=9,
                markeredgecolor=c["surface"], markeredgewidth=1.5, zorder=3)
        ax.plot(o, i, "o", color=c["s1"], markersize=9,
                markeredgecolor=c["surface"], markeredgewidth=1.5, zorder=3)
        if i == 0:
            # Direct labels on the first row instead of a legend box: a legend
            # placed inside these axes renders its keys at plot coordinates and
            # reads as two extra data points.
            ax.annotate("null", xy=(n, i), xytext=(0, -16), ha="center",
                        textcoords="offset points", color=c["s2"], fontsize=9.5,
                        fontweight="600")
            ax.annotate("observed", xy=(o, i), xytext=(0, -16), ha="center",
                        textcoords="offset points", color=c["s1"], fontsize=9.5,
                        fontweight="600")
        ax.annotate(f"+{o - n:.3f}", xy=(o, i), xytext=(10, 0),
                    textcoords="offset points", va="center",
                    color=c["ink2"], fontsize=9.5)
        ax.annotate(f"{n:.3f}", xy=(n, i), xytext=(-8, 0), textcoords="offset points",
                    va="center", ha="right", color=c["ink2"], fontsize=9)
    ax.axvline(0.5, color=c["grid"], linewidth=1.2, zorder=1)
    ax.annotate("0.50", xy=(0.5, len(arms) - 0.35), xytext=(3, 0),
                textcoords="offset points", color=c["muted"], fontsize=8.5)
    ax.set_yticks(y, [LABEL[a] for a in arms])
    ax.invert_yaxis()
    ax.set_xlim(0.44, 0.83)
    ax.set_xlabel("fraction of lead genes with the stronger effect in XX",
                  color=c["ink2"], fontsize=9.5)
    ax.xaxis.grid(True, color=c["grid"], linewidth=1.0)
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)
    frame(ax, c, ygrid=False)
    title(ax, c, "XX-stronger excess survives its null",
          "The null sits just above 0.50; that offset is the sample-size artifact, measured "
          "rather than assumed.")
    finish(fig, ax, c, f"direction-asymmetry.{mode}.svg")


# ---------------------------------------------------------------- figure 8
def fig_contrast_pairs(c, mode):
    ct = DATA["sex_contrast"]
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    vals = [ct[a]["bh05_within_ivw_family"] for a in CELLS]
    x = range(len(CELLS))
    ax.bar(list(x), vals, 0.5, color=c["s1"], zorder=3)
    for i, v in enumerate(vals):
        ax.annotate(f"{v:,}", xy=(i, v), xytext=(0, 4), textcoords="offset points",
                    ha="center", color=c["ink2"], fontsize=9.5)
    ax.set_xticks(list(x), [LABEL[a] for a in CELLS])
    ax.set_ylim(0, max(vals) * 1.2)
    ax.set_ylabel("pairs at BH 5% (ivw family)", color=c["ink2"], fontsize=9.5)
    ax.yaxis.set_major_formatter(THOUSANDS)
    frame(ax, c)
    title(ax, c, "Variant-gene pairs with a significant between-sex effect difference",
          "Counted inside the ivw selection family, the only one where BH applies.")
    finish(fig, ax, c, f"contrast-pairs.{mode}.svg")


# ---------------------------------------------------------------- figure 9
CONTRAST_ORDER = ["t2t_minus_grch38_linear", "t2t_minus_grch38_graph",
                  "graph_minus_linear_grch38", "graph_minus_linear_t2t"]
CONTRAST_LABEL = {
    "t2t_minus_grch38_linear": "T2T − GRCh38  ·  linear",
    "t2t_minus_grch38_graph": "T2T − GRCh38  ·  graph",
    "graph_minus_linear_grch38": "graph − linear  ·  GRCh38",
    "graph_minus_linear_t2t": "graph − linear  ·  T2T",
}


def mb(interval: str) -> str:
    """chr15:82200001-82300000 -> chr15:82.2 Mb"""
    chrom, span = interval.split(":")
    return f"{chrom}:{int(span.split('-')[0]) / 1e6:.1f} Mb"


def fig_reference_vs_aligner(c, mode):
    bc = DATA["hotspots"]["by_contrast"]
    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    for i, key in enumerate(CONTRAST_ORDER):
        q = bc[key]["quantiles_abs_delta_z"]
        colour = c["s1"] if bc[key]["dimension"] == "reference" else c["s2"]
        bx = ax.bxp(
            [{"med": q["0.5"], "q1": q["0.25"], "q3": q["0.75"],
              "whislo": q["0.05"], "whishi": q["0.95"], "fliers": []}],
            positions=[i], widths=0.52, vert=False, patch_artist=True,
            showfliers=False, zorder=3)
        for box in bx["boxes"]:
            box.set(facecolor=colour, edgecolor=colour, linewidth=0)
        for part in ("whiskers", "caps"):
            for a in bx[part]:
                a.set(color=colour, linewidth=1.4)
        for med in bx["medians"]:
            med.set(color=c["surface"], linewidth=1.8)
        ax.annotate(f'{q["0.5"]:.3f}', xy=(q["0.95"], i), xytext=(9, 0),
                    textcoords="offset points", va="center", ha="left",
                    color=c["ink2"], fontsize=9.5)
    ax.set_xscale("log")
    ax.set_yticks(range(len(CONTRAST_ORDER)),
                  [CONTRAST_LABEL[k] for k in CONTRAST_ORDER])
    ax.invert_yaxis()
    ax.set_xlim(0.006, 1.35)
    ax.set_xlabel("per-window mean |Δ Z| across all eligible 100-kb windows (log scale)",
                  color=c["ink2"], fontsize=9.5)
    ax.xaxis.grid(True, color=c["grid"], linewidth=1.0)
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)
    frame(ax, c, ygrid=False)
    title(ax, c, "Changing the reference moves effect estimates ten times more than "
                 "changing the aligner",
          "Box spans the interquartile range, whiskers the 5th to 95th percentile. The two "
          "reference contrasts sit almost entirely to the right of the two aligner "
          "contrasts; labelled values are medians.")
    finish(fig, ax, c, f"reference-vs-aligner.{mode}.svg")


# ---------------------------------------------------------------- figure 10
def fig_hotspot_context(c, mode):
    hs = DATA["hotspots"]
    wins = hs["windows"]
    split = hs["context_split"]
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    groups = [
        (2, "in both workflows", c["s1"], "o", split["replicating"]),
        (1, "one workflow only", c["s2"], "s", split["single_workflow"]),
    ]
    for repl, name, colour, marker, med in groups:
        pts = [w for w in wins if w["replication"] == repl]
        ax.scatter([w["segdup_bp_fraction"] for w in pts],
                   [w["mappability_bp_fraction"] for w in pts],
                   s=68, marker=marker, color=colour, zorder=4,
                   edgecolor=c["surface"], linewidth=1.2,
                   label=f'{name}  (n = {len(pts)})')
        ax.axvline(med["segdup_bp_fraction"], color=colour, linewidth=1.2,
                   linestyle=(0, (4, 3)), zorder=2, alpha=0.85)
        ax.annotate(f'median {med["segdup_bp_fraction"]:.3f}',
                    xy=(med["segdup_bp_fraction"], 1.055), xytext=(0, 0),
                    textcoords="offset points", color=colour, fontsize=9,
                    va="bottom", ha="center")
    # Two windows carry the argument and are named on the plot: the extreme
    # single-workflow case, and the one replicating window that is segdup-rich.
    for w in wins:
        if w["segdup_bp_fraction"] > 0.6:
            ax.annotate(mb(w["interval"]),
                        xy=(w["segdup_bp_fraction"], w["mappability_bp_fraction"]),
                        xytext=(0, -16), textcoords="offset points",
                        color=c["ink2"], fontsize=9, ha="center")
    # Square root, not log: the low-segdup cluster needs spreading and zero has to
    # stay representable.
    ax.set_xscale("function", functions=(lambda v: np.sqrt(np.clip(v, 0, None)),
                                         lambda v: np.square(v)))
    ax.set_xticks([0, 0.01, 0.05, 0.1, 0.2, 0.35, 0.55, 0.85])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlim(-0.004, 0.95)
    ax.set_ylim(0.25, 1.06)
    ax.set_xlabel("segmental-duplication fraction of the window",
                  color=c["ink2"], fontsize=9.5)
    ax.set_ylabel("100-mer mappability fraction", color=c["ink2"], fontsize=9.5)
    ax.xaxis.grid(True, color=c["grid"], linewidth=1.0)
    frame(ax, c)
    title(ax, c, "Windows that need a particular aligner are the hard-sequence ones",
          "Each point is one of the 13 top-ranked reference-discordance windows. Those "
          "recovered by both the linear and the graph workflow sit at low segmental-"
          "duplication content and high mappability.")
    legend(ax, c, loc="lower left")
    finish(fig, ax, c, f"hotspot-context.{mode}.svg")


FIGURES = [fig_missing_burden, fig_pc_sweep, fig_four_arm, fig_stratified,
           fig_interaction, fig_calibration, fig_direction, fig_contrast_pairs,
           fig_reference_vs_aligner, fig_hotspot_context]

if __name__ == "__main__":
    for mode, palette in (("light", LIGHT), ("dark", DARK)):
        for fn in FIGURES:
            fn(palette, mode)
    made = sorted(p.name for p in FIG.glob("*.svg"))
    print(f"wrote {len(made)} figures to {FIG}")
    for m in made:
        print("  ", m)
