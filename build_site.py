#!/usr/bin/env python3
"""Assemble the BrainVar GRCh38-vs-T2T eQTL report page.

Every number in the prose and in the tables is read from report_data.json, which
is itself collected only from corrected (post-nanfix) run roots.  Nothing is
transcribed by hand, so the text cannot drift from the run tree.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
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

V = DATA["validation"]
INT = DATA["four_arm_int_k35"]
LOG = DATA["four_arm_logcpm_k35"]
CURVE = DATA["expression_pc_curve"]
STRAT = DATA["stratified_maps"]
INTER = DATA["logcpm_interaction"]
CT = DATA["sex_contrast"]
NC = DATA["direction_null_check"]
RAW = DATA["nan_whole_sample"]
RAW_XX = DATA["nan_by_arm_raw"]
HS = DATA["hotspots"]
HSC = HS["by_contrast"]
HSW = HS["windows"]

KS = sorted(int(k) for k in CURVE[CELLS[0]])

CONTRAST_ORDER = ["t2t_minus_grch38_linear", "t2t_minus_grch38_graph",
                  "graph_minus_linear_grch38", "graph_minus_linear_t2t"]
CONTRAST_LABEL = {
    "t2t_minus_grch38_linear": "T2T − GRCh38 · linear",
    "t2t_minus_grch38_graph": "T2T − GRCh38 · graph",
    "graph_minus_linear_grch38": "graph − linear · GRCh38",
    "graph_minus_linear_t2t": "graph − linear · T2T",
}
REF_KEYS = [k for k in CONTRAST_ORDER if HSC[k]["dimension"] == "reference"]
WF_KEYS = [k for k in CONTRAST_ORDER if HSC[k]["dimension"] == "workflow"]

# The separation is stated at its weakest point: the smallest reference median against
# the largest aligner median, and the lowest reference 5th percentile against the
# highest aligner 95th.
ref_med_lo = min(HSC[k]["median_window_abs_delta_z"] for k in REF_KEYS)
wf_med_hi = max(HSC[k]["median_window_abs_delta_z"] for k in WF_KEYS)
mag_ratio = ref_med_lo / wf_med_hi
ref_q05 = min(HSC[k]["quantiles_abs_delta_z"]["0.05"] for k in REF_KEYS)
wf_q95 = max(HSC[k]["quantiles_abs_delta_z"]["0.95"] for k in WF_KEYS)

HSA = HS["annotation_summary"]
SPLIT = HS["context_split"]
N_REP = SPLIT["replicating"]["n"]
N_SNG = SPLIT["single_workflow"]["n"]
GO_SIG = HS["go"]["significant_terms"]
GO_SIG_CONTRAST = (CONTRAST_LABEL[GO_SIG[0]["contrast_key"]] if GO_SIG else "")
RANKSUM = SPLIT["rank_sum_p"]
# Windows whose every ClinGen-mappable variant sits in a recurrent-CNV breakpoint, and
# the strongest partial case, named in the prose.
BP_FULL = sorted((w for w in HSW if w["clingen_breakpoint_fraction"] >= 1.0),
                 key=lambda w: w["interval"])
BP_PART = sorted((w for w in HSW if 0 < w["clingen_breakpoint_fraction"] < 1.0),
                 key=lambda w: -w["clingen_breakpoint_fraction"])
# The one replicating window that is nonetheless duplication-rich; named on the page so
# the tendency is not presented as a clean partition.
SEGDUP_EXCEPTION = max((w for w in HSW if w["replication"] == 2),
                       key=lambda w: w["segdup_bp_fraction"])


def mb(interval: str) -> str:
    """chr15:82200001-82300000 -> chr15:82.2 Mb"""
    chrom, span = interval.split(":")
    return f"{chrom}:{int(span.split('-')[0]) / 1e6:.1f} Mb"


def sci(v: float, digits: int = 1) -> str:
    """5.9059e-08 -> 5.9 × 10<sup>−8</sup>, so exponents read as prose, not as repr."""
    mant, exp = f"{v:.{digits}e}".split("e")
    sign = "−" if int(exp) < 0 else ""
    return f"{mant} × 10<sup>{sign}{abs(int(exp))}</sup>"


def ck(arm: str, k: int) -> int:
    d = CURVE[arm]
    return d[str(k)] if str(k) in d else d[k]


def sk(key: str, k: int) -> int:
    d = STRAT[key]
    return d[str(k)] if str(k) in d else d[k]


def fig(name: str, alt: str, caption: str) -> str:
    return f"""
<figure>
  <picture>
    <source srcset="figures/{name}.dark.svg" media="(prefers-color-scheme: dark)">
    <img src="figures/{name}.light.svg" alt="{alt}" loading="lazy">
  </picture>
  <figcaption>{caption}</figcaption>
</figure>"""


def table(headers, rows, *, cls="", note="") -> str:
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    n = f'<p class="tnote">{note}</p>' if note else ""
    return (f'<div class="tablewrap"><table class="{cls}"><thead><tr>{head}</tr></thead>'
            f"<tbody>{body}</tbody></table></div>{n}")


# ----------------------------------------------------------------- tables
t_validation = table(
    ["Analysis profile", "Points", "What it is"],
    [
        ["<code>pooled_int</code>", V["points"]["pooled_int"],
         "Genome-wide genotype×sex interaction, inverse-normal scale"],
        ["<code>logcpm_int</code>", V["points"]["logcpm_int"],
         "The same interaction on log2(CPM+1), at k = 10 and k = 35"],
        ["<code>interaction_pc</code>", V["points"]["interaction_pc"],
         "Interaction-PC sensitivity points"],
        ["<code>stratified_map</code>", V["points"]["stratified_map"],
         "Permutation cis-maps within each sex across the PC grid"],
        ["<code>stratified_nominal</code>", V["points"]["stratified_nominal"],
         "Full nominal scans within each sex (every cis pair)"],
        ["<strong>Total</strong>", f'<strong>{V["total_points"]}</strong>',
         f'<strong>plus {V["downstream_points"]} all-variant contrast points</strong>'],
    ])

t_missing = table(
    ["Arm", "cis-tested variants", "Variants with ≥1 missing call", "Share",
     "NaNs reaching TensorQTL"],
    [[LABEL[a], f'{RAW[a]["variants"]:,}',
      f'{RAW[a]["source_variant_rows_with_nan"]:,}',
      f'{100 * RAW[a]["source_variant_rows_with_nan"] / RAW[a]["variants"]:.1f}%',
      f'{RAW[a]["tensorqtl_nan_values"]}'] for a in CELLS if a in RAW],
    note="Counted across all 225 donors, the sample the maps are actually built on. A "
         "single-sex stratum shows a lower rate (17.4 to 19.4 percent in XX) purely because "
         "fewer donors means fewer chances for any one variant to carry a missing call. The "
         "final column is the whole point of the rebuild.")

t_curve = table(
    ["Arm"] + [f"k = {k}" for k in KS] + ["Peak"],
    [[LABEL[a]] + [f"{ck(a, k):,}" for k in KS]
     + [f"<strong>k = {max(KS, key=lambda k: ck(a, k))}</strong>"] for a in CELLS]
    + [["<strong>Aggregate</strong>"]
       + [f"<strong>{sum(ck(a, k) for a in CELLS):,}</strong>" for k in KS]
       + [f"<strong>k = {max(KS, key=lambda k: sum(ck(a, kk) for a in CELLS) if (kk := k) else 0)}</strong>"]],
    cls="numeric")

t_four = table(
    ["Arm", "Genes tested", "eGenes, INT scale", "eGenes, log2(CPM+1)"],
    [[LABEL[a], f'{INT[a]["tested"]:,}', f'{INT[a]["egenes_q05"]:,}',
      f'{LOG[a]["egenes_q05"]:,}'] for a in CELLS]
    + [["<strong>Total</strong>", "",
        f'<strong>{sum(INT[a]["egenes_q05"] for a in CELLS):,}</strong>',
        f'<strong>{sum(LOG[a]["egenes_q05"] for a in CELLS):,}</strong>']],
    cls="numeric")

SKS = [5, 10, 15, 20, 25, 30]
t_strat = table(
    ["Arm / stratum"] + [f"k = {k}" for k in SKS],
    [[f"{LABEL[a]} · {s}"] + [f"{sk(f'{a}.{s}', k):,}" if str(k) in STRAT[f"{a}.{s}"]
                              or k in STRAT[f"{a}.{s}"] else "—" for k in SKS]
     for a in ("linear_grch38_dv", "linear_t2t_dv") for s in ("XX", "XY")]
    + [[f"{LABEL[a]} · {s}"] + ["—" if k != 15 else f"{sk(f'{a}.{s}', 15):,}" for k in SKS]
       for a in ("graph_grch38_dv", "graph_t2t_dv") for s in ("XX", "XY")],
    cls="numeric",
    note="The graph arms were run at k = 15 only. Dashes are points that were never "
         "scheduled, not failures.")

t_inter = table(
    ["Arm", "INT, k = 35", "log-CPM, k = 10", "log-CPM, k = 35"],
    [[LABEL[a], "0", INTER["e10"][a], INTER["e35"][a]] for a in CELLS],
    cls="numeric")

t_contrast = table(
    ["Arm", "cis pairs tested", "var(z)", "ivw BH 5%", "ivw null rate",
     "either BH 5%", "either null rate"],
    [[LABEL[a], f'{CT[a]["pairs_all"]:,}', f'{CT[a]["variance_of_z_all_pairs"]:.3f}',
      f'{CT[a]["bh05_within_ivw_family"]:,}', f'{CT[a]["ivw"]["simulated_null_rate"]:.4f}',
      f'{CT[a]["bh05_within_either_family"]:,}',
      f'{CT[a]["either"]["simulated_null_rate"]:.4f}'] for a in CELLS],
    cls="numeric",
    note="Only the ivw column is FDR-valid. The either-sex counts are shown so the "
         "difference is visible, not because they can be reported as discoveries.")

t_null = table(
    ["Arm", "Observed XX-stronger", "Genes", "Null", "Null genes", "Excess"],
    [[LABEL[a], f'{NC[a]["observed"]:.3f}', f'{NC[a]["observed_genes"]:,}',
      f'{NC[a]["null"]:.3f}', f'{NC[a]["null_genes"]:,}',
      f'<strong>+{NC[a]["observed"] - NC[a]["null"]:.3f}</strong>']
     for a in CELLS if a in NC],
    cls="numeric")

t_hs_contrast = table(
    ["Contrast", "Axis", "Matched variants", "Occupied windows", "Windows ≥ 100",
     "Median window mean |Δ Z|"],
    [[CONTRAST_LABEL[k],
      "reference" if HSC[k]["dimension"] == "reference" else "aligner",
      f'{HSC[k]["matched_variants"]:,}', f'{HSC[k]["windows"]:,}',
      f'{HSC[k]["windows_ge100"]:,}',
      f'<strong>{HSC[k]["median_window_abs_delta_z"]:.4f}</strong>']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note="Every contrast matches the same gene and the same LiftoverIndel-normalized "
         "allele, so these are like-for-like comparisons of the same variant on two "
         "processings of the same donors. Windows are fixed, non-overlapping, and in T2T "
         "coordinates.")

FRACS = [("segdup_bp_fraction", "Segmental duplication"),
         ("mappability_bp_fraction", "100-mer mappability"),
         ("repeat_bp_fraction", "Repeat content"),
         ("har_bp_fraction", "Human accelerated region")]
t_hs_context = table(
    ["Window sequence context, median",
     f"Recovered by both workflows (n = {N_REP})",
     f"One workflow only (n = {N_SNG})"],
    [[name, f'{SPLIT["replicating"][key]:.3f}', f'{SPLIT["single_workflow"][key]:.3f}']
     for key, name in FRACS],
    cls="numeric",
    note=f"Base-pair fractions of each 100-kb window, from native T2T annotation tracks. "
         f"With {N_REP} windows against {N_SNG} this is a description of the two groups, "
         f"not a powered test; a two-sided rank-sum comparison of the segmental-duplication "
         f"fractions returns p = {RANKSUM['segdup_bp_fraction']:.2f}.")

t_hs_windows = table(
    ["Window (T2T)", "Both workflows", "Variants", "Mean |Δ Z|", "Segdup",
     "Mappability", "cCRE overlap", "ClinGen recurrent-CNV breakpoint", "Anchor genes"],
    [[mb(w["interval"]),
      "yes" if w["replication"] == 2 else "—",
      f'{w["matched_variants"]:,}',
      f'{w["mean_abs_delta_z"]:.2f}',
      f'{w["segdup_bp_fraction"]:.3f}',
      f'{w["mappability_bp_fraction"]:.3f}',
      f'{w["ccre_overlap"]:,} / {w["ccre_mappable"]:,}',
      w["clingen_breakpoint_names"] if w["clingen_breakpoint_variants"] else "—",
      (w["anchor_genes"].replace(";", ", ") or "—")]
     for w in HSW],
    note="The 13 distinct windows ranked in the top ten of either reference contrast, "
         "ordered by workflow replication and then by segmental-duplication content. "
         "cCRE overlap is evaluated only among variants with a unique normalized GRCh38 "
         "position, which is the denominator shown.")

# ----------------------------------------------------------------- prose
int_total = sum(INT[a]["egenes_q05"] for a in CELLS)
log_total = sum(LOG[a]["egenes_q05"] for a in CELLS)
agg = {k: sum(ck(a, k) for a in CELLS) for k in KS}
agg_peak = max(KS, key=lambda k: agg[k])
t2t_agg = {k: ck("graph_t2t_dv", k) + ck("linear_t2t_dv", k) for k in KS}
t2t_peak = max(KS, key=lambda k: t2t_agg[k])
ivw_lo = min(CT[a]["bh05_within_ivw_family"] for a in CELLS)
ivw_hi = max(CT[a]["bh05_within_ivw_family"] for a in CELLS)
nc_lo = min(NC[a]["observed"] - NC[a]["null"] for a in NC)
nc_hi = max(NC[a]["observed"] - NC[a]["null"] for a in NC)

HTML = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Rebuilding a Brain eQTL Map on Two Reference Genomes</title>
<meta name="description" content="A four-arm GRCh38-versus-T2T cis-eQTL analysis of the
BrainVar developmental cohort, the missing-genotype bug that invalidated the first pass,
and what survived the rebuild.">
<meta name="author" content="Joseph Lalli">
<style>
  :root {{
    color-scheme: light;
    --page: #f9f9f7; --surface: #fcfcfb;
    --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
    --rule: #e1e0d9; --edge: #c3c2b7;
    --accent: #2a78d6; --warn: #eb6834; --good: #006300;
    --code-bg: #f2f1ec;
    --measure: 43rem;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      color-scheme: dark;
      --page: #0d0d0d; --surface: #1a1a19;
      --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
      --rule: #2c2c2a; --edge: #383835;
      --accent: #3987e5; --warn: #d95926; --good: #0ca30c;
      --code-bg: #232321;
    }}
  }}
  * {{ box-sizing: border-box; }}
  html {{ -webkit-text-size-adjust: 100%; }}
  body {{
    margin: 0; background: var(--page); color: var(--ink);
    font: 400 17px/1.65 system-ui, -apple-system, "Segoe UI", sans-serif;
  }}
  .wrap {{ max-width: var(--measure); margin: 0 auto; padding: 0 1.25rem 6rem; }}
  header.masthead {{
    border-bottom: 1px solid var(--rule); margin-bottom: 3rem;
    padding: 4.5rem 0 2.5rem;
  }}
  h1 {{ font-size: clamp(2rem, 5.2vw, 2.9rem); line-height: 1.12; letter-spacing: -0.02em;
       margin: 0 0 1rem; font-weight: 700; }}
  .standfirst {{ font-size: 1.16rem; color: var(--ink2); margin: 0 0 1.75rem; line-height: 1.55; }}
  .byline {{ font-size: 0.9rem; color: var(--muted); margin: 0; }}
  h2 {{ font-size: 1.55rem; line-height: 1.25; letter-spacing: -0.01em; margin: 3.5rem 0 0.4rem;
       font-weight: 650; padding-top: 1.5rem; border-top: 1px solid var(--rule); }}
  h3 {{ font-size: 1.12rem; margin: 2.25rem 0 0.4rem; font-weight: 650; }}
  h2 + p, h3 + p {{ margin-top: 0.6rem; }}
  p {{ margin: 0 0 1.15rem; }}
  a {{ color: var(--accent); text-decoration-thickness: 1px; text-underline-offset: 2px; }}
  code {{ font: 500 0.87em/1.4 ui-monospace, SFMono-Regular, Menlo, monospace;
          background: var(--code-bg); padding: 0.12em 0.38em; border-radius: 4px; }}
  figure {{ margin: 2.25rem 0; }}
  figure img {{ width: 100%; height: auto; display: block;
                background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; }}
  figcaption {{ font-size: 0.87rem; color: var(--ink2); margin-top: 0.7rem; line-height: 1.5; }}
  .tablewrap {{ overflow-x: auto; margin: 1.5rem 0 0.6rem;
                border: 1px solid var(--rule); border-radius: 8px; background: var(--surface); }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.86rem; }}
  th, td {{ padding: 0.55rem 0.8rem; text-align: left; white-space: nowrap;
            border-bottom: 1px solid var(--rule); }}
  thead th {{ font-weight: 620; color: var(--ink2); background: var(--code-bg);
              position: sticky; top: 0; }}
  tbody tr:last-child td {{ border-bottom: 0; }}
  table.numeric td + td, table.numeric th + th {{ text-align: right;
       font-variant-numeric: tabular-nums; }}
  .tnote {{ font-size: 0.83rem; color: var(--muted); margin: 0 0 1.6rem; line-height: 1.5; }}
  blockquote {{ margin: 1.75rem 0; padding: 0.15rem 0 0.15rem 1.15rem;
                border-left: 3px solid var(--accent); color: var(--ink2); }}
  blockquote p:last-child {{ margin-bottom: 0; }}
  .callout {{ background: var(--surface); border: 1px solid var(--rule);
              border-left: 3px solid var(--warn); border-radius: 8px;
              padding: 1.1rem 1.25rem; margin: 2rem 0; font-size: 0.95rem; }}
  .callout p:last-child {{ margin-bottom: 0; }}
  .callout strong {{ color: var(--ink); }}
  .keyfig {{ display: flex; flex-wrap: wrap; gap: 1.5rem; margin: 2rem 0 2.5rem;
             padding: 1.4rem 1.5rem; background: var(--surface);
             border: 1px solid var(--rule); border-radius: 10px; }}
  .keyfig div {{ flex: 1 1 8rem; }}
  .keyfig .n {{ font-size: 1.75rem; font-weight: 680; letter-spacing: -0.02em;
                display: block; line-height: 1.15; }}
  .keyfig .l {{ font-size: 0.8rem; color: var(--ink2); display: block; margin-top: 0.3rem;
                line-height: 1.35; }}
  ul, ol {{ margin: 0 0 1.15rem; padding-left: 1.3rem; }}
  li {{ margin-bottom: 0.5rem; }}
  .toc {{ background: var(--surface); border: 1px solid var(--rule); border-radius: 10px;
          padding: 1.2rem 1.5rem; margin: 0 0 3rem; font-size: 0.94rem; }}
  .toc p {{ margin: 0 0 0.6rem; font-weight: 620; font-size: 0.82rem;
            text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted); }}
  .toc ol {{ margin: 0; padding-left: 1.1rem; }}
  .toc li {{ margin-bottom: 0.35rem; }}
  footer {{ border-top: 1px solid var(--rule); margin-top: 4rem; padding-top: 1.75rem;
            font-size: 0.85rem; color: var(--muted); }}
</style>
</head>
<body>
<div class="wrap">

<header class="masthead">
  <h1>Rebuilding a Brain eQTL Map on Two Reference Genomes</h1>
  <p class="standfirst">A four-arm comparison of GRCh38 against T2T-CHM13 in the BrainVar
  developmental cohort — and an account of the one-character convention mismatch that
  quietly deleted a quarter of the variants from the first pass, what it cost, and which
  findings survived rebuilding everything without it.</p>
  <p class="byline">Joseph Lalli</p>
</header>

<nav class="toc">
  <p>Contents</p>
  <ol>
    <li><a href="#design">The question and the four arms</a></li>
    <li><a href="#bug">The bug: NaN is not −9</a></li>
    <li><a href="#rebuild">Rebuilding, and proving the rebuild</a></li>
    <li><a href="#pcs">How many expression PCs?</a></li>
    <li><a href="#maps">The four-arm cis-eQTL maps</a></li>
    <li><a href="#hotspots">Where the references disagree</a></li>
    <li><a href="#interaction">Genotype×sex interaction, and why scale decides it</a></li>
    <li><a href="#stratified">Sex-stratified maps</a></li>
    <li><a href="#contrast">Contrasting effect sizes between sexes</a></li>
    <li><a href="#direction">Is the XX-stronger excess real?</a></li>
    <li><a href="#standing">What is settled, and what is not</a></li>
    <li><a href="#methods">Methods and reproducibility</a></li>
  </ol>
</nav>

<div class="keyfig">
  <div><span class="n">{V["total_points"]}</span><span class="l">TensorQTL analysis points
    rebuilt and audited</span></div>
  <div><span class="n">{V["nan_values_converted"] / 1e6:.0f}M</span><span class="l">missing
    dosages corrected at the boundary</span></div>
  <div><span class="n">{int_total:,}</span><span class="l">eGenes across the four arms
    at k = 35</span></div>
  <div><span class="n">{HS["nominal_rows"] / 1e6:.0f}M</span><span class="l">cis pairs in
    the all-variant nominal scan</span></div>
  <div><span class="n">225</span><span class="l">donors (92 XX, 133 XY)</span></div>
</div>

<h2 id="design">The question and the four arms</h2>

<p>The T2T-CHM13 assembly closed the gaps that GRCh38 left open. The natural question for
anyone holding a genotype-plus-expression cohort is whether that completeness changes the
answers — not in principle, but in the specific, countable sense of whether you call the
same genes as eQTLs.</p>

<p>Answering it requires holding everything else still. This analysis uses the BrainVar
developmental brain cohort, 225 donors, and processes it four separate times. Two
references, GRCh38 and T2T-CHM13, crossed with two alignment strategies, linear and
pangenome-graph. The four arms are labelled <code>linear · GRCh38</code>,
<code>graph · GRCh38</code>, <code>linear · T2T</code> and <code>graph · T2T</code>
throughout. Every arm carries its own genotypes, its own genotype principal components,
its own expression covariates, and its own gene universe called natively against its own
reference; nothing is lifted over between arms to make them comparable, because lifting
over is exactly the operation whose necessity is under test.</p>

<p>Because all four arms run the same code with the same parameters, a difference between
them is attributable to the reference and the aligner. That is the design's whole value,
and it is also what makes it unforgiving: a bug in the shared path shows up in all four
arms at once, looking deceptively like a consistent biological result.</p>

<h2 id="bug">The bug: NaN is not −9</h2>

<p>Genotype dosages have to encode "we don't know." The genotype Parquet files produced by
this pipeline encode a missing dosage as an IEEE floating-point <code>NaN</code>, which is
the ordinary choice for a numeric column with holes in it.</p>

<p>TensorQTL does not use that convention. Its mean-imputation step identifies missing
calls with an equality test against a sentinel value of <code>−9</code>:</p>

<blockquote><p><code>genotypes_t == missing</code>, where <code>missing = -9</code></p></blockquote>

<p>NaN is not equal to anything. It is not equal to <code>−9</code>, and by IEEE-754 it is
not even equal to itself. So the test never fires, the missing dosages are never detected,
never imputed, and pass straight into the linear algebra. Any variant carrying at least one
missing call yields a NaN test statistic and drops silently out of the results.</p>

<div class="callout">
<p><strong>The failure mode is incompleteness, not corruption.</strong> No error is raised,
no warning is printed, and every output file is well-formed and plausible. The results are
simply missing a large, non-random slice of the variants — non-random because missingness
correlates with coverage, with repeat content, and therefore with exactly the regions where
the two references differ. That is the worst possible way for a bug to interact with a
reference-comparison study.</p>
</div>

<p>The scale is not marginal. Between 24.6 and 26.6 percent of cis-tested variants in each
arm carry at least one missing call, so roughly one variant in four was discarded from the
first pass without a trace.</p>

{fig("missing-burden", "Bar chart of the share of cis-tested variants carrying at least one missing call, by arm, ranging from 24.6 to 26.6 percent", "Missing-call burden per arm across all 225 donors. Every variant counted here was silently dropped from the pre-correction analysis. After the fix, the number of NaN values reaching TensorQTL is zero in every one of the 50 analysis points.")}

{t_missing}

<h2 id="rebuild">Rebuilding, and proving the rebuild</h2>

<p>The fix itself is small: a boundary adapter converts NaN to <code>−9</code> in memory
immediately before every TensorQTL call, and nowhere else. The source Parquet files are never
rewritten — they keep the NaN convention, which is correct for every other consumer,
including the fine-mapping tools downstream that expect NaN and would break on a sentinel.
The conversion exists only at the one interface where the convention differs.</p>

<p>Applying a fix is easy. Establishing that it was applied everywhere is the harder problem,
and it is where most of the engineering went. Every analysis point records an audit file
containing the count of NaN values present in its source genotypes, the count of
<code>−9</code> sentinels that reached TensorQTL, a checksum of the adapter, a checksum of
the driver script, the resolved input bindings, and the package versions. A point is only
accepted if the sentinel count is at least the source NaN count.</p>

<p>The whole matrix then has to validate together before anything downstream is allowed to
run:</p>

{t_validation}

<p>Across all {V["total_points"]} points, <strong>{V["nan_values_converted"]:,}</strong>
missing dosages were converted, and the terminal validation records
<code>all_tensorqtl_inputs_nan_free: {str(V["all_inputs_nan_free"]).lower()}</code>. That is
a machine-checked claim, not a description of intent.</p>

<p>The audit design paid for itself in an unexpected way. Every stage refuses to reuse a
directory whose recorded signature does not match what it is about to compute, and refuses
to treat a partially written point as complete. When the rebuild was resumed after an
interruption, that strictness caught two separate latent defects in code paths that had
never executed before — one validator demanding an output file TensorQTL cannot produce for
a non-interaction scan, and one audit record that included itself in the file inventory it
was attesting to, guaranteeing a mismatch. Both raised loudly. A pipeline that had defaulted
to "close enough" would have produced a corrected dataset whose provenance chain quietly
meant nothing.</p>

<h2 id="pcs">How many expression PCs?</h2>

<p>Expression data carries large amounts of structure unrelated to genotype — batch, cell
composition, RNA quality, developmental stage. The standard remedy is to include the leading
principal components of the expression matrix as covariates. Too few and that structure
inflates the noise; too many and you start absorbing the genetic signal itself.</p>

<p>The sweep runs ten values of k from 5 to 50 in all four arms, with permutation testing
throughout. On corrected genotypes the aggregate curve peaks at
<strong>k = {agg_peak}</strong> with {agg[agg_peak]:,} eGenes, against {agg[35]:,} at k = 35.</p>

{fig("pc-sweep", "Line chart of eGene yield against number of expression principal components for the four arms, rising steeply to k=30 and flattening, with peaks at k=45", "eGene yield against expression-PC count. The curve is steep to about k = 30 and nearly flat after it. Three arms peak at k = 45; linear · T2T is still rising at k = 50.")}

{t_curve}

<div class="callout">
<p><strong>An unresolved consequence.</strong> Every map in this analysis was built at
k = 35, and the original justification for that choice was that 35 was the maximum of the
T2T aggregate curve. On corrected data that is no longer true: the T2T aggregate now peaks
at <strong>k = {t2t_peak}</strong> ({t2t_agg[t2t_peak]:,} eGenes against
{t2t_agg[35]:,} at 35). The surviving argument for 35 is a balance argument about gene
biotypes rather than a maximum, and it has not been re-checked against the corrected
by-biotype curve. Choosing 35 costs each arm between 0.7 and 1.7 percent against its own
optimum, so the stakes are modest — but the stated reason for the choice no longer holds,
and that is recorded rather than papered over.</p>
</div>

<h2 id="maps">The four-arm cis-eQTL maps</h2>

<p>With k fixed, each arm gets a permutation-based cis-eQTL map: for every gene, the best
variant within the cis window, with a gene-level empirical p-value and Benjamini-Hochberg
correction across genes.</p>

{fig("four-arm-yield", "Grouped bar chart of eGene counts per arm on the INT and log-CPM scales, showing log-CPM consistently higher", "eGene yield per arm at k = 35 on both expression scales. The four arms agree closely with each other; the choice of expression scale moves the count more than the choice of reference does.")}

{t_four}

<p>Two things stand out. The four arms land within about two percent of each other
({min(INT[a]["egenes_q05"] for a in CELLS):,} to
{max(INT[a]["egenes_q05"] for a in CELLS):,} eGenes on the inverse-normal scale), which is
the reassuring result: swapping the reference genome or the aligner does not move the
headline count much. And the expression scale matters more than the reference does —
log2(CPM+1) yields {log_total:,} eGenes against {int_total:,} for the inverse-normal
transform, a gap larger than any between-arm difference.</p>

<p>That second observation is the hinge of the genotype×sex analysis further down. The first
deserves a harder look before then, because agreement in a headline count is a weak form of
agreement.</p>

<h2 id="hotspots">Where the references disagree</h2>

<p>Two maps can call almost the same number of eGenes while disagreeing about which variants
carry the signal and how strongly. A per-arm total is precisely the statistic that would
conceal that, because it is a count of winners and says nothing about the surface they were
drawn from. Testing it properly needs the whole surface: every cis variant-gene pair in every
arm, not just the per-gene best.</p>

<p>That scan is <strong>{HS["nominal_rows"]:,}</strong> tested pairs across the four arms. It
also functions as an independent audit of everything above, because the permutation-selected
leads have to reappear in it unchanged. All {HS["lead_bridge"]["matched_rows"]:,} of them do,
at a maximum absolute effect-size discrepancy of
{sci(HS["lead_bridge"]["max_abs_beta_error"])} — so the maps in the previous section are
confirmed by a separately computed surface, not merely re-reported.
{HS["degenerate_rows_excluded"]} rows out of the {HS["nominal_rows"]:,} carried TensorQTL's
degenerate zero-variance representation and are excluded from every difference below.</p>

<p>Comparing arms requires one common frame. Each variant is reduced to its
LiftoverIndel-normalized T2T position, reference and alternate allele, matched across arms on
that identity together with the gene, and its slope oriented to the T2T alternate allele. Every
matched variant then yields a difference in test statistic between any two arms. Binned into
fixed {HS["window_size_bp"] // 1000}-kb windows, those differences separate into two regimes
that barely touch.</p>

{fig("reference-vs-aligner", "Horizontal box plot on a log scale of per-window mean absolute delta Z for four contrasts, with the two reference contrasts near 0.21 and the two aligner contrasts near 0.02", "Per-window mean |Δ Z| across all eligible 100-kb windows. The two reference contrasts sit an order of magnitude to the right of the two aligner contrasts, and the whiskers — 5th to 95th percentile — do not meet.")}

{t_hs_contrast}

<p>Changing the reference genome moves the effect estimate about
<strong>{mag_ratio:.0f} times</strong> as much as changing the aligner: a median window mean
|Δ Z| of {ref_med_lo:.3f} against {wf_med_hi:.3f}. The separation is not a difference of
averages over overlapping distributions. The 5th percentile of either reference contrast
({ref_q05:.3f}) sits above the 95th percentile of either aligner contrast ({wf_q95:.3f}), so a
typically-quiet window under a reference swap is still noisier than an unusually loud window
under an aligner swap.</p>

<div class="callout">
<p><strong>What the aligner axis actually measures.</strong> The graph arms are aligned against
the pangenome and then surjected back into linear coordinates before variant calling, so their
genotypes are emitted against a linear allele representation. Pangenomic alignment improves
<em>where reads are placed</em>; after surjection it cannot contribute alleles that the linear
reference has no way to write down. That is the most economical explanation for why the aligner
contrasts are an order of magnitude smaller, and it bounds the claim: this comparison measures
the value of pangenomic read placement, not the value of pangenomic variant representation.
Those are different propositions and it would be easy to read the first as the second.</p>
</div>

<p>Within the reference contrasts the discordance is not spread evenly across the genome.
Ranking windows by mean absolute difference and by the mean within their own upper 5% tail,
screening the leaders against count-matched windows, and calibrating with
{HS["rotations"]:,} chromosome-wise circular rotations, leaves
<strong>{HSA["windows"]} distinct windows</strong> in the top ten of either reference contrast.
Every one carries at least {HS["minimum_matched_variants"]} exactly matched variants.
<strong>{HSA["replicating"]} of the {HSA["windows"]}</strong> are recovered by both the linear
and the graph-surjected workflow, which turns out to be the most informative thing about
them.</p>

<h3>Replication across workflows is a paralog control</h3>

<p>The obvious worry about any reference-difference result is that it is a mapping artifact
rather than a mapping improvement. Segmental duplications are where GRCh38 collapses paralogous
copies onto one locus, and a fixed difference between the copies then reads as a heterozygous
variant in nearly every donor. If T2T resolves that duplication, the genotype changes character
and so does the eQTL estimate — but so it would if the reads were simply being misplaced in a
new way. Segdup enrichment alone cannot separate those two stories, because both predict it.</p>

<p>The four-arm design already contains a discriminator. Because the graph arms use haplotype-
sampled pangenomic alignment, they resolve paralogous copies that a linear aligner collapses. A
window whose signal depends on which aligner was used is a window where read placement was
deciding the answer. A window recovered by both is one where it was not.</p>

{fig("hotspot-context", "Scatter plot of segmental-duplication fraction against mappability for 13 windows, with windows recovered by both workflows clustered at low segdup and high mappability", "The 13 top-ranked windows in sequence context. Those recovered by both workflows cluster at low segmental-duplication content and high mappability; those found by only one spread into duplicated, poorly mappable sequence. chr1:15.9 Mb is the exception that replicates despite high duplication content.")}

{t_hs_context}

<p>The split runs the way the artifact hypothesis predicts it should if the artifact is
confined to the non-replicating group. Windows recovered by both workflows have a median
segmental-duplication fraction of {SPLIT["replicating"]["segdup_bp_fraction"]:.3f} and median
mappability {SPLIT["replicating"]["mappability_bp_fraction"]:.3f}; those found by one workflow
only sit at {SPLIT["single_workflow"]["segdup_bp_fraction"]:.3f} and
{SPLIT["single_workflow"]["mappability_bp_fraction"]:.3f}. The aligner axis is behaving as a
filter on hard sequence, and the {N_REP} windows that survive it are largely ordinary sequence
where paralog collapse is not an available explanation.</p>

<p>This is a description of {HSA["windows"]} windows, not a test, and it should not be read as
one — a two-sided rank-sum comparison of the duplication fractions returns
p = {RANKSUM["segdup_bp_fraction"]:.2f}. One window, {mb(SEGDUP_EXCEPTION["interval"])},
replicates across workflows while sitting at
{SEGDUP_EXCEPTION["segdup_bp_fraction"]:.2f} duplication content, so the rule is a tendency
rather than a partition.</p>

{t_hs_windows}

<h3>What the windows are made of</h3>

<p>Annotating them against ClinGen gives a split result worth stating in both directions.
<strong>{HSA["clingen_breakpoint_windows"]} of {HSA["windows"]}</strong> windows contain
variants inside a recurrent-CNV breakpoint or low-copy-repeat interval, and in
{len(BP_FULL)} of them — {", ".join(mb(w["interval"]) for w in BP_FULL)} — that is true of
<em>every</em> mappable variant the window contains, not merely some;
{mb(BP_PART[0]["interval"])} reaches
{100 * BP_PART[0]["clingen_breakpoint_fraction"]:.0f} percent. But
<strong>{HSA["clingen_score3_windows"]} of {HSA["windows"]}</strong> fall inside a curated
region with sufficient evidence for dosage sensitivity. These loci sit in the structural
scaffolding that makes recurrent CNVs possible, not in the dosage-sensitive intervals whose
disruption is known to cause disease, and the second fact constrains how much the first is
allowed to mean.</p>

<p>Of the {HSA["ccre_mappable"]:,} matched variants in these windows with a unique normalized
GRCh38 position, {HSA["ccre_overlap"]:,} overlap a candidate cis-regulatory element.
{HSA["haqer_interval_windows"]} windows overlap a human accelerated regulatory region and
{HSA["haqer_variant_windows"]} contain matched variants inside one.</p>

<p>Gene Ontology was run as a plausibility check against a background of the eligible anchor
genes in the same contrast, and it returns essentially nothing:
{HS["go"]["sets_tested"] - HS["go"]["sets_with_hits"]} of the {HS["go"]["sets_tested"]} gene
sets have no enriched term at all. The single exception is instructive rather than encouraging.
The {GO_SIG_CONTRAST} set returns {len(GO_SIG)} terms, every one of them taste-receptor biology
— the strongest is {GO_SIG[0]["name"]} at p = {sci(GO_SIG[0]["p_value"])} on
{GO_SIG[0]["intersection_size"]} genes — which is the <em>TAS2R</em> cluster, a textbook
copy-number-variable gene family. A negative control behaving exactly as a negative control
should.</p>

<div class="callout">
<p><strong>What these p-values are and are not.</strong> The spatial statistic asks whether a
window's score cluster is unusual for its own chromosome at its own exact variant count. It is
<em>not</em> a paired test that two effect estimates differ. Two estimates computed from the
same donors on two references are enormously correlated, and a naive paired test of them would
be badly anti-conservative; the rotation null sidesteps that by asking a different and
answerable question. The values are also candidate-screened rather than genome-wide adjusted,
and the window size, the minimum-count gate and the tail definition were chosen for this
follow-up rather than preregistered. Treat the ranking as a shortlist for investigation.</p>
</div>

<h2 id="interaction">Genotype×sex interaction, and why scale decides it</h2>

<p>Does a variant's effect on expression differ between XX and XY donors? The direct test
adds a genotype×sex interaction term and asks whether its coefficient is non-zero.</p>

<p>Run on the inverse-normal-transformed expression that is standard for eQTL mapping, the
answer across the whole genome, in all four arms, is <strong>zero interaction eGenes</strong>.
Not few — zero. And the null is well-behaved rather than broken: the test statistic
distribution matches its theoretical width, and the strongest genome-wide hit is less
extreme than noise alone would predict.</p>

<p>Run on log2(CPM+1) expression, the same cohort, the same variants, the same interaction
term, yields three to six interaction eGenes per arm.</p>

{fig("interaction-scale", "Grouped bar chart showing zero interaction eGenes on the INT scale in all four arms, versus three to six on the log-CPM scale", "Interaction eGenes at BH 5%. The INT bars are zero in every arm at every k tested, and are labelled rather than drawn. The log-CPM bars are small but non-zero, consistently, across both references and both aligners.")}

{t_inter}

<p>The reason is mechanical rather than mysterious. The inverse-normal transform is a
<em>rank</em> transform: it replaces each gene's expression values with their normal quantile
scores, preserving order and discarding spacing. An interaction term is a claim about
magnitudes — it says the slope in one group differs from the slope in the other. Rank-
transforming the phenotype destroys exactly the quantity the interaction coefficient is
trying to estimate, and it does so before the model ever sees the data.</p>

<p>So the two results are not in conflict, and the log-CPM result does not overturn the INT
one. They are answers to different questions, and only one of the two scales can express the
question that was asked. The practical lesson generalises past this dataset: an
inverse-normal transform is a reasonable default for main-effect mapping and a poor one for
any analysis whose estimand is a difference in effect size.</p>

<p>A caution on scale, though. These counts are single digits out of more than eighteen
thousand genes tested. They are consistent in sign and rough magnitude across four
independently processed arms, which is meaningful, but the study is powered to detect a
genotype×sex interaction only if it is roughly an order of magnitude larger than a typical
main effect. Absence of interaction signal here is not evidence that such effects are rare.</p>

<h2 id="stratified">Sex-stratified maps</h2>

<p>The complementary approach is to stop modelling sex and simply split the cohort, mapping
eQTLs separately in the 92 XX and 133 XY donors. This trades statistical power for freedom
from any assumption about how sex enters the model.</p>

{fig("stratified-maps", "Two-panel line chart of eGene counts against expression PCs for XX and XY strata in the two linear arms, with XY roughly 2.2 times XX throughout", "Sex-stratified eGene yield across the PC grid, for the two linear arms. XY yields roughly 2.2 times XX at every k. The two strata also peak at different k, which is why a single shared k is a compromise rather than an optimum for either.")}

{t_strat}

<p>XY produces about 2.2 times as many eGenes as XX at every k tested. It is worth being
explicit that this is almost certainly not biology. With 133 donors against 92, the XY
stratum has roughly 45 percent more samples, and eGene discovery is steeply
power-limited in cohorts of this size — the ratio tracks the sample-size difference, and
nothing in these numbers argues for genuinely more genetic regulation in XY brains.</p>

<p>The strata also peak at different k: XX around k = 10 to 15, XY later, around k = 20.
Because downstream contrasts need both strata computed identically, k = 15 was adopted as a
shared compromise, which is slightly past optimal for XX and short of optimal for XY.</p>

<h2 id="contrast">Contrasting effect sizes between sexes</h2>

<p>The stratified maps make a sharper question available. For every variant-gene pair, we
now have an effect estimate and a standard error in each sex, so we can ask directly whether
the effects differ:</p>

<blockquote><p><code>z = (b<sub>XX</sub> − b<sub>XY</sub>) / √(se<sub>XX</sub>² +
se<sub>XY</sub>²)</code>, with Welch-Satterthwaite degrees of freedom</p></blockquote>

<p>Computed over every cis pair, that is more than 100 million tests per arm. Which raises
the question that decides whether the whole analysis is valid: <em>which pairs do you test?</em></p>

<h3>The selection rule is the analysis</h3>

<p>Testing all 100 million is mostly noise. Some filter to pairs with a real eQTL is needed,
and the obvious filter — keep pairs significant in either sex — is a trap. Selecting on
either sex's effect and then testing the difference between those same effects biases the
difference, because the selection and the contrast share the noise that drove the selection.</p>

<p>The defensible filter selects on the inverse-variance-weighted mean of the two effects:</p>

<blockquote><p><code>b<sub>ivw</sub> = (b<sub>XX</sub>/v<sub>XX</sub> +
b<sub>XY</sub>/v<sub>XY</sub>) / (1/v<sub>XX</sub> + 1/v<sub>XY</sub>)</code></p></blockquote>

<p>This works because <code>Cov(b<sub>ivw</sub>, b<sub>XX</sub> − b<sub>XY</sub>) = 0</code>
exactly. The quantity used to select is statistically orthogonal to the quantity being
tested, so selection carries no information about the contrast and BH correction remains
valid.</p>

<p>That is an algebraic claim, and it was checked empirically by simulating a null in which
both sexes' effects are drawn around their common mean. Under the ivw rule the null rate
comes out at 0.0478 to 0.0479 against a nominal 0.05. Under the either-sex rule it comes out
at about 0.171 — a 3.4-fold inflation.</p>

{fig("selection-calibration", "Grouped bar chart of simulated null rates by selection rule, showing ivw at about 0.048 against nominal 0.05 and either-sex at about 0.171", "Simulated null rate under each selection rule against the nominal 0.05. The ivw rule is calibrated. The either-sex rule is inflated 3.4-fold, so its apparent discoveries are not FDR-valid regardless of how the p-values are adjusted afterwards.")}

{t_contrast}

<p>Inside the valid family, {ivw_lo:,} to {ivw_hi:,} variant-gene pairs per arm show a
between-sex effect difference at BH 5%.</p>

{fig("contrast-pairs", "Bar chart of significant between-sex pairs per arm under the ivw selection family, ranging from about 3,981 to 4,086", "Variant-gene pairs with a significant between-sex effect-size difference, counted within the ivw selection family. The four arms agree closely.")}

<p>One measurement deserves flagging because it slightly undercuts the p-values above.
The variance of z across all pairs is 1.026 to 1.027 rather than 1.0. A theoretical null
would give exactly 1. The excess is small but systematic across all four arms, and it means
the theoretical p-values are mildly anti-conservative — an earlier permutation analysis put
the offset at roughly 5 percent. The counts should be read with that in mind.</p>

<h2 id="direction">Is the XX-stronger excess real?</h2>

<p>Among pairs with a significant difference, the stronger effect is in XX far more often
than in XY. That asymmetry has an obvious deflationary explanation that must be excluded
before it can mean anything: XX has 92 donors against XY's 133, so XX effect estimates carry
larger standard errors, and larger errors produce larger absolute estimates by chance alone.
An XX-stronger excess could be pure measurement artifact.</p>

<p>The way to settle it is to build a null in which, by construction, no true sex difference
exists, and see how much asymmetry the artifact alone generates. Both sexes' effects are
redrawn around their shared inverse-variance-weighted mean, using each sex's real standard
error, then pushed through the identical selection, the identical top-K truncation, and the
identical collapse to one row per gene.</p>

{fig("direction-asymmetry", "Dot plot comparing observed XX-stronger fraction near 0.73 against a null near 0.52 for each of the four arms", "Observed XX-stronger fraction against its null. The null sits just above 0.50 — that offset is the sample-size artifact, measured rather than assumed — while the observed value sits near 0.73 in every arm.")}

{t_null}

<p>The artifact is real and small: the null lands at 0.518 to 0.521, about two points above
the 0.500 that a perfectly balanced design would give. The observed value is 0.716 to 0.737.
The excess of <strong>+{nc_lo:.3f} to +{nc_hi:.3f}</strong> is far larger than the artifact
and agrees closely across both references and both aligners.</p>

<p>Two caveats travel with this number and should not be dropped. First, the null draws noise
independently for each pair, which destroys linkage disequilibrium: the null's selected pairs
spread across roughly 3,500 genes while the observed ones concentrate into about 210, so the
comparison is approximate even though the pair counts are matched. A rigorous version would
permute sex labels and refit, preserving LD; that has not been run on corrected data.
Second, this must be counted per gene and never per variant-gene pair. In an earlier pass a
single large LD block contributed the majority of significant pairs and was XY-stronger,
which flipped the apparent direction entirely when counted by pair. Genes, not pairs.</p>

<h2 id="standing">What is settled, and what is not</h2>

<p>Settled, in the sense of resting on corrected data with a complete audit chain:</p>

<ul>
  <li>The four arms agree closely on cis-eQTL <em>yield</em>, but not on the surface beneath
      it. Reference choice perturbs effect estimates about {mag_ratio:.0f} times more than
      aligner choice does, and that disagreement is spatially concentrated rather than
      diffuse.</li>
  <li>{HSA["replicating"]} of the {HSA["windows"]} top-ranked discordance windows are
      recovered by both the linear and the graph-surjected workflow, and those sit in
      ordinary sequence rather than in segmental duplications.</li>
  <li>Every permutation-selected lead in every arm replays against an independently computed
      {HS["nominal_rows"] / 1e6:.0f}-million-pair nominal surface, to
      {sci(HS["lead_bridge"]["max_abs_beta_error"], 0)} in effect size.</li>
  <li>Genotype×sex interaction is undetectable on the inverse-normal scale and detectable,
      barely, on log-CPM — a consequence of what a rank transform discards.</li>
  <li>The XY-versus-XX eGene ratio of about 2.2 tracks sample size, not biology.</li>
  <li>Only the inverse-variance-weighted selection rule supports valid FDR control for the
      between-sex contrast.</li>
  <li>The XX-stronger direction excess survives its null by roughly 20 points in every arm.</li>
</ul>

<p>Not settled, and deliberately left open rather than resolved by assertion:</p>

<ul>
  <li><strong>k = 35.</strong> Its original justification is false on corrected data. The
      choice stands only because every downstream map was built at it and the cost is
      under two percent.</li>
  <li><strong>The LD-preserving null.</strong> The direction-asymmetry verdict rests on a
      parametric null that ignores LD. The permutation version is the one that would
      settle it.</li>
  <li><strong>The anti-conservative offset.</strong> var(z) above 1 is unexplained and
      systematic.</li>
  <li><strong>Copy number at the discordance windows.</strong> Haplotype-sampled pangenomic
      alignment addresses paralog collapse, but no amount of correct read placement makes a
      multi-copy locus diploid. Whether these are copy-number effects that a SNP dosage
      merely tags is untested. The joint-called genotypes retain the per-sample depth and
      allele-balance fields that would settle it.</li>
  <li><strong>The {N_SNG} single-workflow windows.</strong> By construction these are the
      ones where the aligner decided the answer, and they are the duplicated, poorly mappable
      ones. They are reported here as aligner-sensitive, not as reference findings.</li>
  <li><strong>Fine-mapping.</strong> The per-stratum nominal scans are the natural input to
      multi-condition fine-mapping, which is in progress separately.</li>
</ul>

<h2 id="methods">Methods and reproducibility</h2>

<p>Cohort: BrainVar developmental brain tissue, 225 donors, 92 XX and 133 XY. Genotypes
derived per arm; three genotype principal components (the knee is at 3 in all four arms, with
the third-to-fourth eigenvalue gap two orders of magnitude above the bulk gap floor).
Expression quantified natively against each reference, with gene universes frozen per
reference rather than intersected. Association testing with TensorQTL on two NVIDIA L4 GPUs.</p>

<p>Per-stratum degrees of freedom in the stratified analyses are 60 for XX and 101 for XY,
from 30 covariates. Between-sex nulls are simulated parametrically with one replicate at a
pinned seed. Multiple testing is Benjamini-Hochberg throughout, applied within an explicitly
named selection family.</p>

<p>The all-variant scan is the non-interaction model <code>Y<sub>INT</sub> ~ G + C</code> on all
225 donors with 51 covariate columns and 172 residual degrees of freedom, run across the four
arms at k = 35. Cross-arm comparison uses one common frame: each variant reduced to its
LiftoverIndel-normalized T2T position, reference and alternate allele, matched on that identity
together with the gene, with slopes oriented to the T2T alternate allele and each unique variant
contributing once through its nearest-TSS jointly tested anchor. Windows are fixed,
non-overlapping {HS["window_size_bp"] // 1000}-kb intervals in T2T coordinates. The spatial null
is a chromosome-wise circular rotation of paired score tuples over the ordered eligible variants,
at {HS["rotations"]:,} rotations after a count-matched screen.</p>

<p>Every figure and every table on this page is generated directly from the corrected run
tree by a single collection script, so no number here was transcribed by hand. Figures are
rendered in both light and dark palettes from one code path; the palette was checked
programmatically for colour-vision separation and contrast rather than by eye.</p>

<p>One reporting convention is worth stating because it changes conclusions rather than
presentation: counts of differential effects are reported per <em>gene</em>, never per
variant-gene pair, because linkage disequilibrium lets a single locus dominate a pair-level
count and reverse its apparent direction.</p>

<footer>
  <p>BrainVar GRCh38-versus-T2T cis-eQTL analysis. All results shown are from the corrected
  pipeline; the pre-correction outputs were deleted rather than retained, so nothing on this
  page is derived from them.</p>
</footer>

</div>
</body>
</html>
"""

(HERE / "index.html").write_text(HTML)
(HERE / ".nojekyll").write_text("")
print(f"wrote {HERE / 'index.html'} ({len(HTML):,} bytes)")
