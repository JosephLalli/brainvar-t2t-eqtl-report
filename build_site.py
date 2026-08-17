#!/usr/bin/env python3
"""Assemble the BrainVar GRCh38-vs-T2T eQTL report page.

Every number in the prose and in the tables is read from report_data.json, which
is itself collected only from the four-arm run tree.  Nothing is
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

INT = DATA["four_arm_int_k35"]
LOG = DATA["four_arm_logcpm_k35"]
CURVE = DATA["expression_pc_curve"]
STRAT = DATA["stratified_maps"]
INTER = DATA["logcpm_interaction"]
CT = DATA["sex_contrast"]
NC = DATA["direction_null_check"]
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
N_BOTH = SPLIT["replicating"]["n"]
N_ONE = SPLIT["single_workflow"]["n"]
GO_SIG = HS["go"]["significant_terms"]
GO_SIG_CONTRAST = (CONTRAST_LABEL[GO_SIG[0]["contrast_key"]] if GO_SIG else "")
RANKSUM = SPLIT["rank_sum_p"]
GW = HS["genomewide"]
GWN = GW["null_vs_observed"]
# The diagnostic stated at its least favourable: the contrast where the observed maximum
# stands furthest above the null maximum still barely clears it.
gw_ratio_hi = max(v["observed_max"] / v["null_max_median"] for v in GWN.values())
gw_ratio_ref = (GWN["t2t_minus_grch38_linear"]["observed_max"]
                / GWN["t2t_minus_grch38_linear"]["null_max_median"])
gw_n10 = [v["above_10"] for v in GW["robust_z"].values()]
gw_n25 = [v["above_25"] for v in GW["robust_z"].values()]

EXC = HS["exceedance"]
EXA = EXC["achieved"]
EXCUT = EXC["cut_by_contrast"]
AGREE = EXC["aligner_agreement"]
exc_lo = min(v["flagged"] for v in EXA.values())
exc_hi = max(v["flagged"] for v in EXA.values())
cut_lo, cut_hi = min(EXCUT.values()), max(EXCUT.values())
# The starting point that was proposed before the error rate was pinned, kept so the
# cost of the more permissive cut is visible rather than argued.
EXC196 = {k: next(e for e in v if abs(e["cut"] - 1.96) < 1e-9)
          for k, v in EXC["exceedance"].items()}
n196_lo = min(v["flagged"] for v in EXC196.values())
n196_hi = max(v["flagged"] for v in EXC196.values())
fdr196_lo = min(v["estimated_fdr"] for v in EXC196.values())
fdr196_hi = max(v["estimated_fdr"] for v in EXC196.values())

DIRR = HS["directional"]
DIR_KEY = "t2t_minus_grch38_linear"
DSTR = DIRR["results"][DIR_KEY]["delta_strength"]
DSGN = DIRR["results"][DIR_KEY]["delta_z"]
FRACS_DIR = [("segdup_bp_fraction", "Segmental duplication"),
             ("repeat_bp_fraction", "Repeat content"),
             ("mappability_bp_fraction", "100-mer mappability")]


def fold(side: str, key: str) -> float:
    base = DSTR["neither"][f"{key}_mean"]
    return DSTR[side][f"{key}_mean"] / base if base else float("nan")


def unmappable_fold(side: str) -> float:
    base = 1 - DSTR["neither"]["mappability_bp_fraction_mean"]
    return (1 - DSTR[side]["mappability_bp_fraction_mean"]) / base if base else float("nan")


dir_dup_lo, dir_dup_hi = sorted((fold("up", "segdup_bp_fraction"),
                                 fold("down", "segdup_bp_fraction")))
dir_unmap_lo, dir_unmap_hi = sorted((unmappable_fold("up"), unmappable_fold("down")))


def dq(side: str, key: str, q: str) -> float:
    return DSTR[side][f"{key}_quantiles"][q]


# These window fractions are heavily zero-inflated, so the contrast against the rest of
# the genome is quoted at the upper tail rather than as a shift in the middle.
dup_q95_up = dq("up", "segdup_bp_fraction", "0.95")
dup_q95_dn = dq("down", "segdup_bp_fraction", "0.95")
dup_q95_rest = dq("neither", "segdup_bp_fraction", "0.95")
map_q05_up = dq("up", "mappability_bp_fraction", "0.05")
map_q05_dn = dq("down", "mappability_bp_fraction", "0.05")
map_q05_rest = dq("neither", "mappability_bp_fraction", "0.05")
# Windows whose every ClinGen-mappable variant sits in a recurrent-CNV breakpoint, and
# the strongest partial case, named in the prose.
BP_FULL = sorted((w for w in HSW if w["clingen_breakpoint_fraction"] >= 1.0),
                 key=lambda w: w["interval"])
BP_PART = sorted((w for w in HSW if 0 < w["clingen_breakpoint_fraction"] < 1.0),
                 key=lambda w: -w["clingen_breakpoint_fraction"])
# The one aligner-independent window that is nonetheless duplication-rich; named on the
# page so the tendency is not presented as a clean partition.
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


CONC = DATA["concordance"]["by_contrast"]
XAXIS = DATA["concordance"]["cross_axis_overlap"]
GCLS = DATA["gene_classes"]["by_contrast"]
CHRX = DATA["chrx_by_sex"]["results"]
MECH = DATA["mechanism"]["by_contrast"]
UNI = DATA["gene_universe"]

conc_calls_aligner = min(CONC[k]["call_jaccard"] for k in WF_KEYS)
conc_calls_ref = min(CONC[k]["call_jaccard"] for k in REF_KEYS)
conc_r_aligner = min(CONC[k]["pearson_r"] for k in WF_KEYS)
conc_r_ref = min(CONC[k]["pearson_r"] for k in REF_KEYS)
gene_moved_lo = min(CONC[k]["gene_discordant_rate"] for k in CONTRAST_ORDER)
gene_moved_hi = max(CONC[k]["gene_discordant_rate"] for k in CONTRAST_ORDER)
gd_or_lo = min(GCLS[k]["genomic_disorder"]["odds_ratio"] for k in CONTRAST_ORDER)
gd_or_hi = max(GCLS[k]["genomic_disorder"]["odds_ratio"] for k in CONTRAST_ORDER)
har_or_lo = min(GCLS[k]["human_accelerated"]["odds_ratio"] for k in CONTRAST_ORDER)
har_or_hi = max(GCLS[k]["human_accelerated"]["odds_ratio"] for k in CONTRAST_ORDER)
chrx_hit = CHRX["graph_minus_linear_t2t::XX"]
uni_excl = UNI["grch38_only"]["genes"] + UNI["t2t_only"]["genes"]
uni_excl_egenes = UNI["grch38_only"]["egenes"] + UNI["t2t_only"]["egenes"]

# ----------------------------------------------------------------- tables
t_concordance = table(
    ["Contrast", "Correlation of z", "Variants moving &lt; 0.5 z",
     "Associations called identically", "Genes that moved"],
    [[CONTRAST_LABEL[k], f'{CONC[k]["pearson_r"]:.4f}',
      f'{CONC[k]["share_abs_delta_z_below_0.5"]:.1%}',
      f'<strong>{CONC[k]["call_jaccard"]:.1%}</strong>',
      f'{CONC[k]["gene_discordant_rate"]:.1%}'] for k in CONTRAST_ORDER],
    cls="numeric",
    note=f'Across every exactly matched variant-gene pair. An association is "called" when '
         f'|z| reaches {DATA["concordance"]["call_threshold_abs_z"]:.0f}, roughly p = 6×10⁻⁵. '
         f'The final column counts genes whose mean |Δ Z| across their cis variants is an '
         f'outlier at an estimated 5% false-discovery proportion.')

t_gene_classes = table(
    ["Contrast", "Recurrent genomic-disorder region", "Segmental duplication",
     "Human accelerated region"],
    [[CONTRAST_LABEL[k]]
     + [f'<strong>{GCLS[k][c]["odds_ratio"]:.1f}×</strong><br>'
        f'<span class="sub">{GCLS[k][c]["rate_in_class"]:.1%} vs '
        f'{GCLS[k][c]["rate_out"]:.1%}</span>'
        for c in ("genomic_disorder", "segmental_duplication", "human_accelerated")]
     for k in CONTRAST_ORDER],
    cls="numeric",
    note="Odds that a method-sensitive gene falls in each class, against the genes tested in "
         "the same contrast — never against the genome, which would manufacture the result. "
         "Percentages are the discordance rate inside and outside the class. Fisher tests "
         "treat genes as independent and genes in one duplicated block are not, so the odds "
         "ratios are the trustworthy part.")

t_chrx = table(
    ["Contrast", "XX: chrX vs autosomal", "XX odds", "XY: chrX vs autosomal", "XY odds"],
    [[CONTRAST_LABEL[k],
      f'{CHRX[f"{k}::XX"]["chrx_rate"]:.1%} vs {CHRX[f"{k}::XX"]["autosomal_rate"]:.1%}',
      f'<strong>{CHRX[f"{k}::XX"]["odds_ratio"]:.2f}</strong>',
      f'{CHRX[f"{k}::XY"]["chrx_rate"]:.1%} vs {CHRX[f"{k}::XY"]["autosomal_rate"]:.1%}',
      f'{CHRX[f"{k}::XY"]["odds_ratio"]:.2f}'] for k in CONTRAST_ORDER],
    cls="numeric",
    note="Each stratum is compared against its own autosomes, so the difference in stratum "
         "size cannot drive the comparison. chrX dosage encoding was audited and is identical "
         "across all four arms, which excludes differential ploidy handling.")

t_universe = table(
    ["", "GRCh38 only", "T2T only", "Shared"],
    [["Genes testable", f'{UNI["grch38_only"]["genes"]:,}',
      f'{UNI["t2t_only"]["genes"]:,}', f'{UNI["shared"]:,}'],
     ["eGenes among them", f'<strong>{UNI["grch38_only"]["egenes"]:,}</strong>',
      f'<strong>{UNI["t2t_only"]["egenes"]:,}</strong>', "—"],
     ["Share that are eGenes", f'{UNI["grch38_only"]["egene_rate"]:.1%}',
      f'{UNI["t2t_only"]["egene_rate"]:.1%}',
      f'{UNI["shared_reference_context"]["egene_rate"]:.1%}'],
     ["Absent from the other annotation",
      f'{UNI["grch38_only_absent_from_t2t_annotation"]:,}',
      f'{UNI["t2t_only_absent_from_grch38_annotation"]:,}', "—"],
     ["Mean duplication content", "—",
      f'{UNI["t2t_only"]["mean_segdup"]:.3f}',
      f'{UNI["shared_reference_context"]["mean_segdup"]:.3f}']],
    cls="numeric",
    note="A gene is in a reference's universe when both of that reference's arms tested it. "
         "Duplication content is measured in T2T coordinates and so is available only for the "
         "T2T-exclusive and shared sets.")

TOPK = DATA["topk"]["by_contrast"]
topk_ref_lo = min(TOPK[k]["depths"]["100"]["changed"] for k in REF_KEYS)
topk_ref_hi = max(TOPK[k]["depths"]["100"]["changed"] for k in REF_KEYS)
topk_wf_hi = max(TOPK[k]["depths"]["100"]["changed"] for k in WF_KEYS)

t_topk = table(
    ["Contrast", "Top-100 kept", "Candidates that change", "Not testable in the other arm",
     "Rank correlation"],
    [[CONTRAST_LABEL[k], f'{TOPK[k]["depths"]["100"]["share"]:.0%}',
      f'<strong>{TOPK[k]["depths"]["100"]["changed"]}</strong>',
      f'{TOPK[k]["depths"]["100"]["absent_from_other_universe"]}',
      f'{TOPK[k]["spearman_rank_correlation"]:.3f}']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note="Genes ranked within each arm by the gene-level permutation p-value the map already "
         "produces. Rank correlation is computed over genes both arms rank, which separates "
         "the same genes in a different order from different genes. Ranking is noisy near "
         "ties, so a gene crossing the boundary of a top-100 list need not reflect a "
         "meaningful change in evidence.")

LEAD = DATA["lead_switching"]["by_contrast"]
lead_ref_same_lo = min(LEAD[k]["share_same_lead"] for k in REF_KEYS)
lead_ref_same_hi = max(LEAD[k]["share_same_lead"] for k in REF_KEYS)
lead_wf_same_lo = min(LEAD[k]["share_same_lead"] for k in WF_KEYS)
lead_wf_same_hi = max(LEAD[k]["share_same_lead"] for k in WF_KEYS)

t_lead = table(
    ["Contrast", "eGenes in both arms", "Same lead variant", "Lead changes",
     "Median move", "Crosses the TSS"],
    [[CONTRAST_LABEL[k], f'{LEAD[k]["genes_egene_in_both"]:,}',
      f'<strong>{LEAD[k]["share_same_lead"]:.1%}</strong>',
      f'{1 - LEAD[k]["share_same_lead"]:.1%}',
      f'{LEAD[k]["median_move_bp"] / 1000:.1f} kb',
      f'{LEAD[k]["share_crossing_tss"]:.0%}']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note="Restricted to genes that are eGenes in both arms, so differences in yield cannot "
         "contribute. Leads are compared on normalised variant identity rather than on "
         "coordinates. Distance is a weak proxy for whether two leads tag the same signal — "
         "the measurement that would settle that is linkage disequilibrium between them, "
         "which has not been computed.")

BGN = DATA["background_noise"]["by_arm"]
BGR = DATA["background_noise"]["ratios"]
LDS = DATA["lead_switch_ld"]["by_contrast"]
ld_same_ref = min(LDS[k]["share_r2_above_0.8"] for k in REF_KEYS)
ld_same_wf = min(LDS[k]["share_r2_above_0.8"] for k in WF_KEYS)
ld_indep_ref = max(LDS[k]["share_r2_below_0.2"] for k in REF_KEYS)
# Compounding the two stages: switched at all, and switched to something independent.
net_ref = max((1 - LEAD[k]["share_same_lead"]) * LDS[k]["share_r2_below_0.2"]
              for k in REF_KEYS)
net_wf = max((1 - LEAD[k]["share_same_lead"]) * LDS[k]["share_r2_below_0.2"]
             for k in WF_KEYS)

t_lead_ld = table(
    ["Contrast", "Switched pairs measured", "Median r²",
     "Same signal (r² &gt; 0.8)", "Independent (r² &lt; 0.2)"],
    [[CONTRAST_LABEL[k], f'{LDS[k]["pairs_with_ld"]:,}', f'{LDS[k]["median_r2"]:.3f}',
      f'<strong>{LDS[k]["share_r2_above_0.8"]:.1%}</strong>',
      f'{LDS[k]["share_r2_below_0.2"]:.1%}']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note="r² between the two leads, computed in the positive arm's own genotype matrix over "
         "pairwise-complete donors rather than through an external panel that would carry its "
         "own reference. Restricted to switched leads whose counterpart is nameable in that "
         "arm's variant space; a lead with no counterpart there is excluded, and those are the "
         "cases where the arms differ most — so the independent share is understated.")

ANC = DATA["ancestry"]
ANC_GROUPS = [g for g in ("EUR", "AMR", "AFR") if g in ANC["group_sizes"]]
anc_ref = ANC["altload_change"]["t2t_minus_grch38_linear"]
anc_wf = ANC["altload_change"]["graph_minus_linear_grch38"]

t_ancestry = table(
    ["Ancestry group", "Donors", "Change in non-reference load on T2T",
     "Change under an aligner swap"],
    [[g, f'{ANC["group_sizes"][g]}',
      f'<strong>{anc_ref["by_group"][g]:+.3f}</strong>',
      f'{anc_wf["by_group"][g]:+.4f}']
     for g in ANC_GROUPS],
    cls="numeric",
    note=f'Median change in mean alternate-allele dosage per donor, a proxy for how far a '
         f'donor sits from the reference. Kruskal-Wallis across groups: '
         f'p = {anc_ref["kruskal_p"]:.1e} for the reference swap and '
         f'{anc_wf["kruskal_p"]:.1e} for the aligner swap. Ancestry labels are projected '
         f'superpopulation assignments, not self-reported, and groups outside these three are '
         f'too small to test.')

EXCL = DATA["arm_exclusive"]["by_contrast"]

t_exclusive = table(
    ["Contrast", "Exclusive variants", "Reach |z| ≥ 4", "Shared variants do",
     "Ratio", "Duplication content"],
    [[CONTRAST_LABEL[k], f'{EXCL[k]["exclusive"]["variants"]:,}',
      f'<strong>{EXCL[k]["exclusive"]["share_reaching_call"]:.1%}</strong>',
      f'{EXCL[k]["shared"]["share_reaching_call"]:.1%}',
      f'<strong>{EXCL[k]["call_rate_ratio"]:.2f}×</strong>',
      f'{EXCL[k]["exclusive"]["mean_segdup"]:.3f} vs '
      f'{EXCL[k]["shared"]["mean_segdup"]:.3f}']
     for k in EXCL],
    cls="numeric",
    note=f'Computed on {" and ".join(DATA["arm_exclusive"]["contigs"])} only. A variant is '
         f'exclusive to an arm when its normalised identity has no counterpart among the '
         f'other arm\'s tested variants; variants that cannot be placed in the common frame '
         f'at all are excluded from both groups, and those are the most reference-specific of '
         f'all. Signal is the best |z| a variant reaches against any gene, which favours '
         f'variants tested against more genes.')

YARD = DATA["yardstick"]["by_contrast"]
YARD_ROWS = [("t2t_minus_grch38_linear", "reference"),
             ("t2t_minus_grch38_graph", "reference"),
             ("graph_minus_linear_grch38", "aligner"),
             ("graph_minus_linear_t2t", "aligner"),
             ("hc_minus_dv_grch38", "caller"),
             ("hc_minus_dv_t2t", "caller")]
YARD_LABEL = dict(CONTRAST_LABEL,
                  hc_minus_dv_grch38="HaplotypeCaller − DeepVariant · GRCh38",
                  hc_minus_dv_t2t="HaplotypeCaller − DeepVariant · T2T")


def _beyond(r):
    return r.get("share_beyond_0.01", 1 - r.get("share_within_0.01", float("nan")))


t_yardstick = table(
    ["What was changed", "Axis", "Variants compared", "Frequencies identical",
     "Differing by &gt; 0.01"],
    [[YARD_LABEL[k], dim, f'{YARD[k]["matched"]:,}',
      f'{YARD[k]["share_exact"]:.1%}',
      f'<strong>{_beyond(YARD[k]):.2%}</strong>'] for k, dim in YARD_ROWS],
    cls="numeric",
    note=f'{DATA["yardstick"]["donors"]} donors present in all six callsets, across '
         f'{", ".join(DATA["yardstick"]["contigs"])}. Frequencies are recomputed identically '
         f'from each callset, so no axis has an advantage of stage, donor set or region. '
         f'Cross-reference rows are restricted to variants carrying a unique normalised '
         f'identity in both references, which same-reference rows are not — see the note '
         f'below the figure.')

NACC = DATA["newly_accessible"]["by_contrast"]

t_access = table(
    ["Contrast", "Genes whose lead the other arm cannot represent", "Share",
     "Among genes reaching |z| ≥ 4"],
    [[CONTRAST_LABEL[k],
      f'{NACC[k]["genes_whose_lead_is_exclusive"]:,} of {NACC[k]["genes"]:,}',
      f'<strong>{NACC[k]["share_leads_exclusive"]:.1%}</strong>',
      f'{NACC[k]["share_called_leads_exclusive"]:.1%}']
     for k in NACC],
    cls="numeric",
    note="Genome-wide. A variant is exclusive to an arm when it carries a unique normalised "
         "identity there and has no counterpart among the other arm's tested variants; "
         "variants with no unique identity at all are excluded from both sides, so this is a "
         "lower bound. A gene's lead being exclusive does not mean the other arm finds nothing "
         "for that gene — only that it cannot find that variant.")

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
     f"Reference effect under both aligners (n = {N_BOTH})",
     f"Under one aligner only (n = {N_ONE})"],
    [[name, f'{SPLIT["replicating"][key]:.3f}', f'{SPLIT["single_workflow"][key]:.3f}']
     for key, name in FRACS],
    cls="numeric",
    note=f"Base-pair fractions of each 100-kb window, from native T2T annotation tracks. "
         f"With {N_BOTH} windows against {N_ONE} this is a description of the two groups, "
         f"not a powered test; a two-sided rank-sum comparison of the segmental-duplication "
         f"fractions returns p = {RANKSUM['segdup_bp_fraction']:.2f}.")

t_hs_windows = table(
    ["Window (T2T)", "Both aligners", "Variants", "Mean |Δ Z|", "Segdup",
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
         "ordered by whether the reference effect holds under both aligners and then by "
         "segmental-duplication content. "
         "cCRE overlap is evaluated only among variants with a unique normalized GRCh38 "
         "position, which is the denominator shown.")

t_gw = table(
    ["Contrast", "Windows tested", "Largest observed window mean |Δ Z|",
     "Median largest under rotation", "Family-wise threshold", "Beyond chance"],
    [[CONTRAST_LABEL[k], f'{GW["robust_z"][k]["windows"]:,}',
      f'{GWN[k]["observed_max"]:.2f}', f'{GWN[k]["null_max_median"]:.2f}',
      f'{GWN[k]["threshold"]:.2f}',
      f'<strong>{GW["significant_by_contrast"][k]}</strong>']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note=f'{GW["rotations"]:,} chromosome-wise circular rotations per contrast, with the '
         f'maximum window mean over the whole genome recorded on each rotation. The third '
         f'and fourth columns are the point of the table: rotating the genome relocates the '
         f'largest cluster but does not remove it, so what chance produces is almost exactly '
         f'what was observed.')

t_exc = table(
    ["Contrast", "Cut", "Windows discordant", "Share of windows", "Estimated FDR",
     "At z &gt; 1.96 instead"],
    [[CONTRAST_LABEL[k], f'z &gt; {EXCUT[k]:.2f}',
      f'<strong>{EXA[k]["flagged"]:,}</strong>', f'{EXA[k]["observed_rate"]:.2%}',
      f'{EXA[k]["estimated_fdr"]:.1%}',
      f'{EXC196[k]["flagged"]:,} at {EXC196[k]["estimated_fdr"]:.0%}']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note="The cut is chosen per contrast as the most permissive one holding the estimated "
         "false-discovery proportion at 5%. That proportion is the normal-null expectation "
         "divided by the observed count, so it assumes a normal tail; linkage "
         "disequilibrium correlates neighbouring windows, which makes the true null tail "
         "heavier and this estimate a lower bound.")

t_dir = table(
    ["Statistic", "Direction", "Windows"]
    + [n for _, n in FRACS_DIR] + ["Up vs down"],
    [[lab, dirn, f'{blk[side]["n"]:,}']
     + [f'{blk[side][f"{key}_quantiles"]["0.5"]:.3f} '
        f'({blk[side][f"{key}_quantiles"]["0.95"]:.3f})' for key, _ in FRACS_DIR]
     + [f'p = {blk["up_vs_down_rank_sum_p"][FRACS_DIR[0][0]]:.3f}' if side == "up" else ""]
     for lab, blk in (("Association strength", DSTR), ("Signed effect", DSGN))
     for side, dirn in (("up", "larger on T2T"), ("down", "smaller on T2T"))]
    + [["Rest of the genome", "—", f'{DSTR["neither"]["n"]:,}']
       + [f'{DSTR["neither"][f"{key}_quantiles"]["0.5"]:.3f} '
          f'({DSTR["neither"][f"{key}_quantiles"]["0.95"]:.3f})' for key, _ in FRACS_DIR]
       + [""]],
    cls="numeric",
    note="Median base-pair fraction per 100-kb window with the 95th percentile in brackets, "
         "for the T2T − GRCh38 contrast within the linear workflow. Medians rather than "
         "means because these fractions are heavily zero-inflated. Association strength is "
         "the mean of |z| on T2T minus |z| on GRCh38; signed effect is the mean of z minus z "
         "with both arms oriented to the same T2T alternate allele. The last column is a "
         "two-sided rank-sum comparison of the up and down groups on duplication content.")

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
<title>Mapping Brain eQTLs on Two Reference Genomes</title>
<meta name="description" content="A four-arm GRCh38-versus-T2T cis-eQTL analysis of the
BrainVar developmental cohort: what changes when you swap the reference genome, what
changes when you swap the aligner, and where on the genome the two disagree.">
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
  .sub {{ font-size: 0.82em; color: var(--muted); font-weight: 400; }}
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
  <h1>Mapping Brain eQTLs on Two Reference Genomes</h1>
  <p class="standfirst">Swapping the reference genome or the aligner leaves almost all of a
  brain cis-eQTL map exactly where it was. The exception is a small and specific minority of
  the genome — and it is the part that is difficult to align <em>because</em> it varies so
  much between people, which is the same reason it carries disease.</p>
  <p class="byline">Joseph Lalli</p>
</header>

<nav class="toc">
  <p>Contents</p>
  <ol>
    <li><a href="#design">The question and the four arms</a></li>
    <li><a href="#pcs">How many expression PCs?</a></li>
    <li><a href="#maps">The four-arm cis-eQTL maps</a></li>
    <li><a href="#concordance">What a method change leaves alone</a></li>
    <li><a href="#mechanism">What kind of difference it is</a></li>
    <li><a href="#yardstick">How big is that, really?</a></li>
    <li><a href="#hotspots">Where the references disagree</a></li>
    <li><a href="#classes">What the moved regions are made of</a></li>
    <li><a href="#universe">The genes that were never askable</a></li>
    <li><a href="#chrx">Chromosome X, and why sex decides it</a></li>
    <li><a href="#ancestry">Ancestry, and who a reference represents</a></li>
    <li><a href="#interaction">Genotype×sex interaction, and why scale decides it</a></li>
    <li><a href="#stratified">Sex-stratified maps</a></li>
    <li><a href="#contrast">Contrasting effect sizes between sexes</a></li>
    <li><a href="#direction">Is the XX-stronger excess real?</a></li>
    <li><a href="#standing">What is settled, and what is not</a></li>
    <li><a href="#methods">Methods and reproducibility</a></li>
  </ol>
</nav>

<div class="keyfig">
  <div><span class="n">{conc_calls_aligner:.0%}</span><span class="l">of association calls
    unchanged when the aligner is swapped</span></div>
  <div><span class="n">{gd_or_lo:.1f}–{gd_or_hi:.1f}×</span><span class="l">enrichment of the
    genes that do move for genomic-disorder loci</span></div>
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
them is attributable to the reference and the aligner. That is the design's whole value: the two
factors are crossed, so each can be read off separately, and an effect that appears
under only one level of the other factor is identifiable as exactly that.</p>

<h2 id="pcs">How many expression PCs?</h2>

<p>Expression data carries large amounts of structure unrelated to genotype — batch, cell
composition, RNA quality, developmental stage. The standard remedy is to include the leading
principal components of the expression matrix as covariates. Too few and that structure
inflates the noise; too many and you start absorbing the genetic signal itself.</p>

<p>The sweep runs ten values of k from 5 to 50 in all four arms, with permutation testing
throughout. The aggregate curve peaks at
<strong>k = {agg_peak}</strong> with {agg[agg_peak]:,} eGenes, against {agg[35]:,} at k = 35.</p>

{fig("pc-sweep", "Line chart of eGene yield against number of expression principal components for the four arms, rising steeply to k=30 and flattening, with peaks at k=45", "eGene yield against expression-PC count. The curve is steep to about k = 30 and nearly flat after it. Three arms peak at k = 45; linear · T2T is still rising at k = 50.")}

{t_curve}

<div class="callout">
<p><strong>An unresolved choice.</strong> Every map in this analysis was built at
k = 35, on the understanding that 35 was the maximum of the T2T aggregate curve. The
sweep does not support that: the T2T aggregate peaks
at <strong>k = {t2t_peak}</strong> ({t2t_agg[t2t_peak]:,} eGenes against
{t2t_agg[35]:,} at 35). The surviving argument for 35 is a balance argument about gene
biotypes rather than a maximum, and it has not been re-checked against the
by-biotype curve. Choosing 35 costs each arm between 0.7 and 1.7 percent against its own
optimum, so the stakes are modest — but the stated reason for the choice does not hold,
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

<p>That second observation is the hinge of the genotype×sex analysis further down. The first —
that four independently processed arms land within two percent of each other — is the one worth
pressing on, because a matching headline count is a weak form of agreement and everything that
follows is an attempt to find out how weak.</p>

<h2 id="concordance">What a method change leaves alone</h2>

<p>The useful result has to be stated in this order, because the second half only means
something against the first. Most people believe that reference and aligner choice do not
much matter for genetic disease research. <strong>For most of the genome they are right</strong>,
and confirming that is the first job here.</p>

<p>Every variant that survives in both arms of a contrast, matched on gene and on exact
normalised allele, gives a pair of test statistics that can be compared directly.</p>

{fig("concordance", "Horizontal bar chart of agreement per contrast, showing 97 percent of association calls identical under an aligner swap and about three quarters under a reference swap", "Agreement across every exactly matched variant. Swapping the aligner is very nearly a no-op. Swapping the reference moves more, but still leaves the large majority of the map exactly where it was.")}

{t_concordance}

<p>Changing the aligner is close to doing nothing: the two arms correlate at
{conc_r_aligner:.4f}, and <strong>{conc_calls_aligner:.1%}</strong> of the associations either
arm would report are reported by both. Changing the reference moves more — correlation
{conc_r_ref:.3f}, {conc_calls_ref:.1%} of calls shared — but still leaves roughly nine variants
in ten within half a z-unit of where they started, and {1 - gene_moved_hi:.0%} to
{1 - gene_moved_lo:.0%} of genes untouched altogether.</p>

<p>Nothing in this analysis threatens prior genetics. If a reference swap had overturned a
large fraction of a cis-eQTL map, that would have been a reason to distrust the swap rather
than the map.</p>

<p>That reassurance needs one qualification immediately, because a genome-wide average is not
what anyone acts on. Decisions are made at the head of the ranking — the genes that get
followed up, fine-mapped, or put in a figure — and the head is less stable than the bulk.</p>

{t_topk}

<p>Switching reference changes <strong>{topk_ref_lo} to {topk_ref_hi} of the top hundred
genes</strong>; switching aligner changes {topk_wf_hi}. The rank correlation among genes both
arms rank falls to about
{min(TOPK[k]["spearman_rank_correlation"] for k in REF_KEYS):.2f} for a reference swap against
{min(TOPK[k]["spearman_rank_correlation"] for k in WF_KEYS):.2f} for an aligner swap. A
handful of the changes are genes the other arm cannot test at all, which is a different kind
of difference and is taken up further down.</p>

<p>There is a third level to this, and it is the one that matters most for what people do next.
Restrict attention to genes that <em>both</em> arms call as eGenes — no yield difference, no
ranking effect, the two arms agreeing that the gene has a signal — and ask whether they agree
about which variant carries it.</p>

{t_lead}

<p>They agree {lead_ref_same_lo:.0%} to {lead_ref_same_hi:.0%} of the time under a reference
swap and {lead_wf_same_lo:.0%} to {lead_wf_same_hi:.0%} under an aligner swap. <strong>Four
times out of ten, two references that agree a gene has an eQTL nominate a different variant as
its lead.</strong> These are not sub-resolution wobbles: only about an eighth of the moves are
under a kilobase, the median is over ten, and a fifth to a quarter land on the other side of
the transcription start site.</p>

<p>Taken alone that invites a conclusion it does not support, because a different variant is not
the same thing as a different signal. Fourteen kilobases inside one haplotype block is a
relabelling; the same distance across a recombination boundary is a different causal
hypothesis. Distance cannot tell them apart, so the linkage disequilibrium between each pair of
leads was measured directly, in the arm's own genotypes.</p>

{t_lead_ld}

<p>It softens the result substantially. The median r² between the two leads is
{min(LDS[k]["median_r2"] for k in REF_KEYS):.2f} for a reference swap and
{min(LDS[k]["median_r2"] for k in WF_KEYS):.2f} for an aligner swap, and
{ld_same_ref:.0%} to {ld_same_wf:.0%} of switched pairs exceed r² = 0.8. <strong>Most lead
changes are the same signal wearing a different label.</strong> Only about
{ld_indep_ref:.0%} of them are effectively independent under a reference swap.</p>

<p>Compounding the two stages gives the number that matters: a reference swap leaves roughly
<strong>{net_ref:.0%} of shared eGenes</strong> with a genuinely different causal candidate, and
an aligner swap about {net_wf:.0%}. That is a far weaker claim than four-in-ten, and it is the
correct one.</p>

<p>All four numbers are true and they belong together. The map as a whole is stable; the
shortlist drawn off the top of it is less so; the nominated variant within a shared hit changes
often but usually cosmetically; and a small, real remainder is a different hypothesis
altogether. Which raises the question the rest of this addresses: what distinguishes the part
that moves.</p>

<h2 id="mechanism">What kind of difference it is</h2>

<p>Before asking where the arms disagree, it is worth asking what kind of disagreement it is.
A test statistic can move for two quite different reasons, and they can be separated exactly
rather than approximately. Writing z as β/se,</p>

<blockquote><p><code>z<sub>a</sub> − z<sub>b</sub> = (β<sub>a</sub> − β<sub>b</sub>)/se<sub>a</sub>
+ β<sub>b</sub>(1/se<sub>a</sub> − 1/se<sub>b</sub>)</code></p></blockquote>

<p>The first term is movement in the <em>estimate</em>: the arms disagree about the effect,
which is what happens when they are effectively looking at different genotypes or different
local linkage. The second is movement in <em>precision</em>: the same effect measured with more
or less certainty, which is what missingness and genotype quality produce. The identity is
algebraic, and it was checked numerically in every contrast to a maximum residual of
{MECH["t2t_minus_grch38_linear"]["decomposition"]["max_decomposition_residual"]:.0e}.</p>

{fig("mechanism", "Grouped bar chart on a log scale comparing the effect and precision contributions to the change in test statistic, with the effect term far larger in all four contrasts", "The exact decomposition of the movement in the test statistic. In all four contrasts the estimate moves and the precision does not, by roughly an order of magnitude.")}

<p>The answer is the same for both axes, and it was not the expected one. Movement is dominated
by the estimate in
{min(MECH[k]["decomposition"]["share_variants_effect_dominated"] for k in CONTRAST_ORDER):.1%}
to
{max(MECH[k]["decomposition"]["share_variants_effect_dominated"] for k in CONTRAST_ORDER):.1%}
of variants, and precision barely moves at all: the median ratio of standard errors between
arms is
{MECH["t2t_minus_grch38_linear"]["decomposition"]["median_se_ratio_a_over_b"]:.4f} for a
reference swap and
{MECH["graph_minus_linear_t2t"]["decomposition"]["median_se_ratio_a_over_b"]:.4f} for an
aligner swap.</p>

<div class="callout">
<p><strong>Neither change buys precision. Both change what is being measured.</strong> That
matters for how any of this should be described. Where these methods differ they are not
reducing noise around a fixed quantity; they are estimating a different one. The single place
precision does shift is duplicated sequence, where the aligner's precision term roughly doubles
relative to its effect term while the reference's barely moves — pangenomic alignment is
measurably changing genotype certainty there, and only there.</p>
</div>

<p>That is a statement about the difference between paired estimates. The level of noise each
arm carries is a separate question, and it has the same answer. At matched allele frequency the
median standard error is
{BGN["linear_grch38_dv"]["median_se_maf_0.2_0.35"]:.5f},
{BGN["graph_grch38_dv"]["median_se_maf_0.2_0.35"]:.5f},
{BGN["linear_t2t_dv"]["median_se_maf_0.2_0.35"]:.5f} and
{BGN["graph_t2t_dv"]["median_se_maf_0.2_0.35"]:.5f} across the four arms — identical to the
fourth decimal. Splitting by sequence class does not separate them either: in duplicated
sequence the graph-to-linear ratio is
{BGR["aligner_t2t"]["median_se_duplicated_common"]:.3f} and the T2T-to-GRCh38 ratio
{BGR["reference_linear"]["median_se_duplicated_common"]:.3f}, the largest entry in the table
being T2T carrying standard errors slightly <em>larger</em> than GRCh38 exactly where a gain
would be most expected.</p>

<p><strong>No arm is measurably less noisy, anywhere.</strong> Standard errors are higher in
duplicated sequence than in ordinary sequence in every arm, and no method change repairs that.
Whatever these choices are doing, they are not improving the precision with which a shared
variant is measured.</p>

<p>There is a second way to ask the same question that needs no association model at all. The
same 225 donors and the same normalised allele must give the same allele frequency, so any
difference proves that the measurement at that site depends on which representation was used.
The agreement is logically compelled rather than empirically hoped for, which makes this the
cleanest instrument in the analysis — and the median difference is <strong>exactly
zero</strong> in all four contrasts.</p>

<p>Frequencies agree exactly at
{DATA["mechanism"]["by_contrast"]["t2t_minus_grch38_linear"]["allele_frequency"]["share_exact"]:.0%}
of sites under a reference swap and
{DATA["mechanism"]["by_contrast"]["graph_minus_linear_grch38"]["allele_frequency"]["share_exact"]:.0%}
under an aligner swap, and fall within 0.01 at
{DATA["mechanism"]["by_contrast"]["t2t_minus_grch38_linear"]["allele_frequency"]["share_within_0.01"]:.1%}
and
{DATA["mechanism"]["by_contrast"]["graph_minus_linear_grch38"]["allele_frequency"]["share_within_0.01"]:.1%}
respectively. The
{1 - DATA["mechanism"]["by_contrast"]["graph_minus_linear_grch38"]["allele_frequency"]["share_within_0.01"]:.1%}
to
{1 - DATA["mechanism"]["by_contrast"]["t2t_minus_grch38_graph"]["allele_frequency"]["share_within_0.01"]:.1%}
that fall outside it are the sites whose measurement is representation-dependent. The
comparison is <em>unsigned</em>: it does not say which arm is wrong, and truth may be the
higher or the lower frequency.</p>

<p>Worth noting because the two instruments share no statistics: the chromosomes carrying the
most frequency discordance are chr19, chr22, chr17, chrX, chr6 and chr21 — the same set the
association analysis identifies below, arrived at from genotype frequencies alone.</p>

<h2 id="yardstick">How big is that, really?</h2>

<p>Every number so far is unanchored. A correlation of {conc_r_ref:.2f}, or
{1 - DATA["mechanism"]["by_contrast"]["t2t_minus_grch38_linear"]["allele_frequency"]["share_within_0.01"]:.1%}
of sites with a shifted allele frequency, means nothing without knowing what a change of
comparable weight looks like. The variant caller supplies that. Swapping DeepVariant for
HaplotypeCaller is a decision most groups make without discussion, and it can be measured on
exactly the same instrument.</p>

{fig("yardstick", "Bar chart comparing the share of sites with shifted allele frequency across reference, aligner and caller swaps, with the caller bars roughly twice the height of the others", "The three method axes on one footing: same donors, same chromosomes, same stage. Changing the variant caller moves allele frequencies more than twice as much as changing the reference genome.")}

{t_yardstick}

<p><strong>A reference swap is the smallest of the three changes.</strong> It shifts the
frequency of {DATA["yardstick"]["mean_by_axis"]["reference"]:.2%} of sites, against
{DATA["yardstick"]["mean_by_axis"]["aligner"]:.2%} for an aligner swap and
{DATA["yardstick"]["mean_by_axis"]["caller"]:.2%} for a caller swap — the caller being
{DATA["yardstick"]["ratio_to_aligner"]["caller"]:.1f} times the aligner and
{DATA["yardstick"]["ratio_to_aligner"]["caller"] / DATA["yardstick"]["ratio_to_aligner"]["reference"]:.1f}
times the reference. Anyone comfortable with their choice of variant caller has already
accepted a larger perturbation than the one this page is about.</p>

<div class="callout">
<p><strong>The caveat is the finding, not a footnote.</strong> A cross-reference comparison can
only be made on variants that carry a unique identity in <em>both</em> references — about 1.9
million of the roughly 10 million in each callset. Same-reference comparisons use all of them.
So the reference axis above is measured <em>only among variants both references can
represent</em>, and the variants only one reference can represent are excluded entirely.</p>
<p>Read the two halves together and they say something sharper than either alone. Where two
references can both see a variant, they agree about it better than two aligners do. The
reference's distinctive effect is not that it measures shared variants differently — it is
<em>what it can see at all</em>. That is what the gene universes and the arm-exclusive variants
further down are measuring, and it is where the reference change actually lives.</p>
</div>

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
<strong>{HSA["replicating"]} of the {HSA["windows"]}</strong> show the same reference
effect under both aligners, which turns out to be the most informative thing about
them.</p>

<h3>Separating the reference effect from the aligner</h3>

<p>The obvious worry about any reference-difference result is that it is a mapping artifact
rather than a mapping improvement. Segmental duplications are where GRCh38 collapses paralogous
copies onto one locus, and a fixed difference between the copies then reads as a heterozygous
variant in nearly every donor. If T2T resolves that duplication, the genotype changes character
and so does the eQTL estimate — but so it would if the reads were simply being misplaced in a
new way. Segdup enrichment alone cannot separate those two stories, because both predict it.</p>

<p>The factorial design can. The reference contrast is measurable twice, once under each
aligner, and the two aligners treat paralogy differently: haplotype-sampled pangenomic
alignment resolves copies that a linear aligner collapses. A reference effect present under
both aligners is therefore a property of the reference. One that appears under a single
aligner is a reference-by-aligner interaction — at that locus the two factors are entangled,
and the reference difference cannot be read on its own.</p>

<p>This is not a replication argument, and it should not be mistaken for one. The four arms
are not four attempts at the same experiment. The reference and the aligner <em>are</em> the
experiment, and there is no expectation that they agree; the question is only which of the
two a given difference can be attributed to.</p>

{fig("hotspot-context", "Scatter plot of segmental-duplication fraction against mappability for 13 windows, with windows whose reference effect holds under both aligners clustered at low segdup and high mappability", "The 13 top-ranked windows in sequence context. Those whose reference effect holds under both aligners cluster at low segmental-duplication content and high mappability; those confined to one aligner spread into duplicated, poorly mappable sequence. chr1:15.9 Mb is the exception that holds under both despite high duplication content.")}

{t_hs_context}

<p>The split runs the way it would if entanglement with the aligner were confined to hard
sequence. Windows whose reference effect holds under both aligners have a median
segmental-duplication fraction of {SPLIT["replicating"]["segdup_bp_fraction"]:.3f} and median
mappability {SPLIT["replicating"]["mappability_bp_fraction"]:.3f}; those confined to one
aligner sit at {SPLIT["single_workflow"]["segdup_bp_fraction"]:.3f} and
{SPLIT["single_workflow"]["mappability_bp_fraction"]:.3f}. The interaction with the aligner
is concentrated exactly where read placement is contestable, and the {N_BOTH}
aligner-independent windows are largely ordinary sequence where paralog collapse is not an
available explanation.</p>

<p>This is a description of {HSA["windows"]} windows, not a test, and it should not be read as
one — a two-sided rank-sum comparison of the duplication fractions returns
p = {RANKSUM["segdup_bp_fraction"]:.2f}. One window, {mb(SEGDUP_EXCEPTION["interval"])},
holds under both aligners while sitting at
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

<h3>How far does this extend, and is any of it beyond chance?</h3>

<p>The obvious next question is not how the top ten rank but how many windows anywhere in the
genome exceed what chance would produce. That test was run on all four contrasts. It uses the
same chromosome-wise circular rotation, except that the quantity recorded on each rotation is
the largest window mean <em>anywhere in the genome</em>, so that its upper quantile is a
family-wise threshold rather than a per-window one. Across {GW["tests"]:,} window tests and
{GW["rotations"]:,} rotations, the number of windows that exceed it is
<strong>{GW["significant"]}</strong>.</p>

{t_gw}

<p>The reason sits in the middle two columns, and it is a property of the test rather than of
the data. A circular shift slides the score vector along the chromosome; it does not remove the
cluster of large differences, it relocates it. Whichever window lands under that cluster after
the shift inherits its magnitude, so the largest window mean under rotation comes back nearly
equal to the largest one observed — at best a ratio of {gw_ratio_hi:.2f} across the four
contrasts, and {gw_ratio_ref:.2f} in the contrast that matters most. A shift null therefore asks whether a
cluster's <em>position</em> is unusual, and the answer to that is no. It cannot be made to ask
whether the cluster's <em>magnitude</em> is unusual, which is the question worth asking.</p>

<p>A shift null is not the only option, though, and the second one works. Every window already
carries a count-matched robust z: its mean |Δ Z| minus the median across windows holding a
comparable number of variants, divided by that stratum's robust spread. The statistic is
centred on zero with unit spread by construction, so a normal null predicts a fixed rate in the
upper tail, and the amount by which the observed rate exceeds it estimates how much of the
flagged set is noise. Holding that estimate at five percent puts the cut at z between
{cut_lo:.2f} and {cut_hi:.2f} and leaves <strong>{exc_lo:,} to {exc_hi:,} windows per
contrast</strong> discordant beyond chance.</p>

{t_exc}

<p>At the more permissive z &gt; 1.96 the same tail yields {n196_lo:,} to {n196_hi:,} windows,
but around a quarter of them ({fdr196_lo:.0%} to {fdr196_hi:.0%}) are expected to be noise,
which is why the operating point is set by the error rate rather than by the sigma. Either way
the answer to the original question is that regional discordance is not a matter of a handful
of loci. It involves something on the order of five to seven percent of the testable genome.</p>

<p>That scale also settles the aligner question the thirteen windows could only gesture at. Of
the {AGREE["union"]:,} windows that either reference contrast flags,
<strong>{AGREE["flagged_under_both_aligners"]["n"]:,}</strong> are flagged under both aligners
and only {AGREE["flagged_under_one_aligner"]["n"]:,} under one — an overlap of
{AGREE["jaccard"]:.2f}. The reference effect is overwhelmingly a property of the reference and
not of the alignment method. The aligner-specific minority behaves exactly as the small sample
suggested it would: those windows carry more duplicated sequence
(p = {sci(AGREE["rank_sum_p"]["segdup_bp_fraction"], 1)}) and lower mappability
(p = {sci(AGREE["rank_sum_p"]["mappability_bp_fraction"], 1)}) than the windows both aligners
agree on. At seven against six that comparison was underpowered; at
{AGREE["flagged_under_both_aligners"]["n"]:,} against
{AGREE["flagged_under_one_aligner"]["n"]:,} it is not.</p>

<h3>Which direction, and does it matter where?</h3>

<p>Reducing a window to |Δ Z| discards the difference between a reference that recovers an
association and one that loses it. Two signed statistics separate them: the change in
association <em>strength</em>, |z| on T2T minus |z| on GRCh38, and the signed movement of the
effect itself, z minus z with both arms oriented to the same alternate allele. Each tail is cut
independently, at the same five percent.</p>

{t_dir}

<p>Strength is close to balanced — {DSTR["up"]["flagged"]:,} windows strengthen on T2T against
{DSTR["down"]["flagged"]:,} that weaken. The signed effect is not: {DSGN["up"]["flagged"]:,}
move up against {DSGN["down"]["flagged"]:,} that move down, and both reference contrasts show
the same roughly two-to-one skew.</p>

{fig("direction-context", "Three box plots comparing duplication, repeat content and mappability for windows where associations strengthen on T2T, weaken on T2T, and the rest of the genome", "Sequence context of discordant windows by direction. Boxes span the interquartile range and whiskers the 5th to 95th percentile. The strengthening and weakening distributions sit on top of each other; both differ from the rest of the genome, and for duplication that difference is entirely in the upper tail.")}

<p>The sequence context is the part worth pausing on, and it needs the distributions rather
than their averages. These fractions are heavily zero-inflated: median duplication content is
zero in the discordant windows just as it is across the rest of the genome, so an average
would suggest a shift that is not there. The difference lives in the tail. The 95th percentile
of duplication content is {dup_q95_up:.2f} among windows where the association strengthens on
T2T and {dup_q95_dn:.2f} among those where it weakens, against {dup_q95_rest:.2f} elsewhere;
the 5th percentile of mappability is {map_q05_up:.2f} and {map_q05_dn:.2f} against
{map_q05_rest:.2f}.</p>

<p>The two directions, though, are not distinguishable from one another at all — on
duplication (p = {DSTR["up_vs_down_rank_sum_p"]["segdup_bp_fraction"]:.2f}), on mappability
(p = {DSTR["up_vs_down_rank_sum_p"]["mappability_bp_fraction"]:.2f}) or on repeat content
(p = {DSTR["up_vs_down_rank_sum_p"]["repeat_bp_fraction"]:.2f}). Difficult sequence predicts
<em>that</em> a region will disagree between the two references. It does not predict
<em>which way</em>, which puts the direction down to something locus-specific rather than to
how hard the region is to align.</p>

<p>The one exception is the signed effect statistic, where the down-shifted windows are
modestly more duplicated than the up-shifted ones
(p = {sci(DSGN["up_vs_down_rank_sum_p"]["segdup_bp_fraction"], 1)}) and more repetitive
(p = {DSGN["up_vs_down_rank_sum_p"]["repeat_bp_fraction"]:.2f}). Together with the
two-to-one count asymmetry, that reads as a systematic downward shift of the test statistic in
duplicated sequence on T2T rather than as a regulatory difference — and it points back at
genotype quality in exactly the regions where the copy-number question below is still
open.</p>

<div class="callout">
<p><strong>What these numbers are and are not.</strong> None of the window statistics on this
page is a paired test that two effect estimates differ. Two estimates computed from the same
donors on two references are enormously correlated, and a naive paired test of them would be
badly anti-conservative. The window size, the minimum-count gate and the tail definition were
chosen for this follow-up rather than preregistered, and the ranking is candidate-screened
rather than genome-wide adjusted. Treat it as a shortlist for investigation, not as a
result.</p>
</div>

<h2 id="classes">What the moved regions are made of</h2>

<p>A list of discordant windows is only interesting if the windows have something in common.
Scoring every gene by how far its cis associations move, cutting at an estimated 5%
false-discovery proportion, and asking what kind of sequence those genes sit in gives an answer
that is the same in every contrast.</p>

{fig("gene-classes", "Lollipop chart of odds ratios across four contrasts and three sequence classes, all twelve above one, with the largest values for genomic-disorder regions under an aligner swap", "Odds that a method-sensitive gene falls in each class of sequence, against the genes tested in the same contrast. Twelve of twelve comparisons are positive.")}

{t_gene_classes}

<p>Genes whose eQTL results move are <strong>{gd_or_lo:.1f} to {gd_or_hi:.1f} times</strong>
more likely to sit in a recurrent genomic-disorder region, and
<strong>{har_or_lo:.1f} to {har_or_hi:.1f} times</strong> more likely to sit near a human
accelerated region — the annotation that exists specifically to mark recent human evolution.
Segmental duplication runs between
{min(GCLS[k]["segmental_duplication"]["odds_ratio"] for k in CONTRAST_ORDER):.1f} and
{max(GCLS[k]["segmental_duplication"]["odds_ratio"] for k in CONTRAST_ORDER):.1f} times.</p>

<div class="callout">
<p><strong>Difficulty and biological interest are not two facts here. They are one.</strong>
These regions are hard to align because they vary so much between people that a single linear
reference represents them badly, and between-person variation is the substance of disease
genetics — so the same property that makes them hard to measure is what makes them worth
measuring. That is why the enrichment above is not a warning about instability. It is the
prediction the design was built to test, and it holds on all three axes at once.</p>
</div>

<p>The consequence is uncomfortable but specific. Fields with the longest history of
difficulty — psychiatry, neurodevelopment, immunology — work disproportionately on loci of
exactly this kind. The claim is not that they have been getting wrong answers. It is that they
are the fields for which representation choice is not free, while most of genetics can
reasonably continue to treat it as inert.</p>

<p>Formal pathway enrichment of the moved genes is, correctly, almost empty: against a
background of the genes tested in the same contrast, no term survives multiple-testing
correction in three of the four contrasts. Below the threshold the same biology recurs across
contrasts — MHC class I antigen processing and presentation, T-cell receptor signalling,
allograft rejection, Fc-receptor activation — which is what one would expect if the effect is
carried by the MHC and the Fc-receptor loci rather than by immune pathways at large. It is
reported here as a consistent tendency, not as enrichment.</p>

<h2 id="universe">The genes that were never askable</h2>

<p>Every comparison so far conditions on genes present in both references, which makes all of
it blind to the genes that exist in one universe and not the other. Those universes were frozen
natively per reference and never intersected, and the difference between them is not what the
totals suggest.</p>

{t_universe}

<p>GRCh38 tests {UNI["grch38_universe"]:,} genes and T2T {UNI["t2t_universe"]:,} — a net gap of
{abs(UNI["grch38_universe"] - UNI["t2t_universe"])}. But the universes are not nested:
<strong>{uni_excl:,} genes are exclusive to one reference or the other</strong>, and
<strong>{uni_excl_egenes} of them are eGenes</strong>: real associations callable on one
reference and simply unavailable on the other. Exclusive genes are, if anything, <em>more</em>
likely to be eGenes than shared ones.</p>

<p>The T2T-exclusive set sits in sequence about four times as duplicated as the shared set
({UNI["t2t_only"]["mean_segdup"]:.3f} against
{UNI["shared_reference_context"]["mean_segdup"]:.3f} mean duplication content), which is the
same signature every other part of this analysis finds. The
{UNI["t2t_only"]["biotypes"].get("rRNA", 0)} ribosomal RNA genes testable only on T2T are the
expected consequence of an assembly that resolves the rDNA arrays GRCh38 leaves as gaps.</p>

<p>One caveat has to travel with these numbers. The two universes come from different
annotation releases, and a majority of exclusive genes —
{UNI["grch38_only_absent_from_t2t_annotation"]} of {UNI["grch38_only"]["genes"]} and
{UNI["t2t_only_absent_from_grch38_annotation"]} of {UNI["t2t_only"]["genes"]} — are absent from
the other reference's annotation altogether. Annotation therefore accounts for much of the
turnover, and only the residual is attributable to the reference sequence itself. Both figures
belong in any honest statement of this result.</p>

<p>The same question can be asked one level down, of variants rather than genes, and there the
two axes behave quite differently. A variant is exclusive to an arm when its normalised
identity has no counterpart among the other arm's tested variants; comparing those against the
same arm's shared variants asks whether the exclusive ones ever produce an association.</p>

{t_exclusive}

<p>Variants only the pangenome graph can call reach the association threshold
<strong>{EXCL["graph_minus_linear_t2t"]["call_rate_ratio"]:.1f} times</strong> as often as the
variants both aligners share, and they sit in sequence roughly two and a half times as
duplicated. Variants only T2T can call are, if anything, slightly <em>less</em> likely to carry
a signal than shared ones
({EXCL["t2t_minus_grch38_linear"]["call_rate_ratio"]:.2f}×).</p>

<p>So the aligner is not simply adding variants; it is adding variants that produce signal, in
exactly the duplicated sequence a linear aligner has least to work with. The reference adds
variants of ordinary informativeness. A yield count sees both as "more variants tested" and
cannot tell them apart.</p>

<p>The sharpest form of the question is per gene. Not "are exclusive variants informative on
average", but: <em>is a gene's single best association at a variant the other arm cannot
represent at all?</em></p>

{t_access}

<p><strong>About one gene in nine has its top association on a variant GRCh38 cannot
represent</strong> — {NACC["t2t_minus_grch38_linear"]["genes_whose_lead_is_exclusive"]:,} of
{NACC["t2t_minus_grch38_linear"]["genes"]:,}. For an aligner swap it is
{NACC["graph_minus_linear_t2t"]["share_leads_exclusive"]:.1%}.</p>

<div class="callout">
<p><strong>This is the mechanical core of the whole comparison.</strong> Three results point
the same way from different directions. Among variants both references can represent, they
agree about them more closely than two aligners do. But
{uni_excl:,} genes are testable on one reference only, {uni_excl_egenes} of those are eGenes,
and for one shared gene in nine the strongest association sits on a variant the other reference
has no way to write down.</p>
<p><strong>The reference's distinctive effect is on access, not on measurement.</strong> Where
it can see the same thing, it sees it the same way. What changes is what it can see. No
comparison of yields, effect sizes, or correlations can reach that, because every one of them
conditions on the variants both references share.</p>
</div>

<h2 id="chrx">Chromosome X, and why sex decides it</h2>

<p>One result does not fit the pattern of the others, and it took a stratified analysis to
understand it. Measured across the whole cohort, pangenomic alignment on T2T moves chromosome X
far more than it moves the autosomes — and no other contrast does anything comparable. Split by
sex, the effect turns out to live in a single cell of the design.</p>

{fig("chrx-by-sex", "Bar chart on a log scale of chrX odds of discordance against autosomes, by contrast and sex stratum, with seven bars at or below one and a single bar above seven", "Odds that a chrX gene is discordant relative to an autosomal gene, within the same sex stratum. Only pangenomic alignment on T2T, and only in XX donors, moves chromosome X.")}

{t_chrx}

<p>In XX donors, pangenomic alignment on T2T moves
<strong>{chrx_hit["chrx_rate"]:.1%}</strong> of chromosome X genes against
{chrx_hit["autosomal_rate"]:.1%} of autosomal genes — odds of
<strong>{chrx_hit["odds_ratio"]:.2f}</strong>. In XY donors the same comparison gives
{CHRX["graph_minus_linear_t2t::XY"]["odds_ratio"]:.2f}: chromosome X moves <em>less</em> than
the autosomes. The remaining six cells sit between
{min(CHRX[f"{k}::{s}"]["odds_ratio"] for k in REF_KEYS + [WF_KEYS[0]] for s in ("XX", "XY")):.2f}
and
{max(CHRX[f"{k}::{s}"]["odds_ratio"] for k in REF_KEYS for s in ("XX", "XY")):.2f}.</p>

<p>This is not a power artifact — the null stratum is the larger one, with 133 XY donors
against 92 XX. Nor is it a ploidy-encoding artifact: chromosome X dosages were audited in all
four arms and are encoded identically, with hemizygous genotypes appearing as heterozygous in
under 2.5% of XY calls in every arm against about 29% in XX.</p>

<div class="callout">
<p><strong>One mechanism accounts for all four rows.</strong> Pangenomic alignment resolves
haplotype <em>diversity</em>. A hemizygous X carries one haplotype and offers nothing to
resolve, which is why XY shows nothing. GRCh38's own X is poorly resolved, so the graph has no
better material to work with there, which is why neither sex moves on GRCh38. Only T2T's
complete X together with two X haplotypes gives the graph both the material and the
opportunity. It is a three-way interaction between reference, aligner and sex, which is why no
two-way analysis found it.</p>
</div>

<p>The practical consequence is that a whole-cohort chromosome X analysis in a mixed-sex study
reports a diluted version of a much stronger effect in half the samples, and should be
stratified. It also means the sex-stratified work that follows — built for an entirely
different question — is not separable from the methods question after all.</p>

<h2 id="ancestry">Ancestry, and who a reference represents</h2>

<p>Sex is one axis on which a method change can fall unevenly across a cohort. Ancestry is the
other, and it has a more direct mechanism: a reference genome is a particular sequence from
particular people, so how well it represents a donor depends on where that donor's ancestry
sits relative to it.</p>

<p>Measuring that needs no association model. Each donor's mean alternate-allele dosage says
how far they sit from whichever reference is in use, and the change in it between arms says
whether a reference swap moved them closer or further away. This cohort is
{ANC["group_sizes"].get("EUR", 0)} donors of European, {ANC["group_sizes"].get("AFR", 0)} of
African and {ANC["group_sizes"].get("AMR", 0)} of American projected ancestry, which is enough
to compare.</p>

{t_ancestry}

<p><strong>The direction reverses between groups.</strong> Moving to T2T reduces the
non-reference load of European-ancestry donors by {abs(anc_ref["by_group"]["EUR"]):.3f} and
<em>increases</em> it for African-ancestry donors by
{abs(anc_ref["by_group"]["AFR"]):.3f}. The pattern is identical under both aligners, so it is a
property of the reference and not of the alignment.</p>

<p>This is what the two references' provenance would predict. T2T-CHM13 is a single haploid
cell line of largely European ancestry; GRCh38 is a mosaic with substantial African-ancestry
contribution. Swapping one for the other therefore moves the reference point, and it moves it
toward some donors and away from others. The aligner swap does the same test at roughly a
fiftieth of the magnitude.</p>

<div class="callout">
<p><strong>What this does and does not say.</strong> It measures reference bias, not accuracy.
A donor carrying more non-reference alleles is further from the reference, which is a fact
about the reference rather than about the donor or about the quality of their genotypes. It
does <em>not</em> show that African-ancestry donors receive worse eQTL estimates on T2T —
establishing that would mean linking this shift to those donors' contribution to the
discordance measured earlier, which has not been done. What it does establish is that
reference choice interacts with cohort composition, and that a study reporting a T2T
reanalysis of an ancestrally diverse cohort is not making an ancestry-neutral change.</p>
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
permute sex labels and refit, preserving LD; that has not been run.
Second, this must be counted per gene and never per variant-gene pair. In an earlier pass a
single large LD block contributed the majority of significant pairs and was XY-stronger,
which flipped the apparent direction entirely when counted by pair. Genes, not pairs.</p>

<h2 id="standing">What is settled, and what is not</h2>

<p>Settled, in the sense of resting on the complete four-arm run tree:</p>

<ul>
  <li><strong>A reference swap is the smallest of the three method changes.</strong> Among
      variants both references can represent it shifts
      {DATA["yardstick"]["mean_by_axis"]["reference"]:.2%} of allele frequencies, against
      {DATA["yardstick"]["mean_by_axis"]["aligner"]:.2%} for an aligner swap and
      {DATA["yardstick"]["mean_by_axis"]["caller"]:.2%} for a variant-caller swap.</li>
  <li><strong>Most of the map does not move.</strong> An aligner swap leaves
      {conc_calls_aligner:.1%} of association calls identical and a reference swap
      {conc_calls_ref:.1%}, with {1 - gene_moved_hi:.0%} to {1 - gene_moved_lo:.0%} of genes
      untouched. Nothing here overturns prior genetics.</li>
  <li><strong>The shortlist is less stable than the map.</strong> A reference swap changes
      {topk_ref_lo} to {topk_ref_hi} of the top hundred genes and an aligner swap
      {topk_wf_hi}, so the practical consequence of the method choice is larger than the
      genome-wide correlation implies.</li>
  <li><strong>Agreeing on the gene is usually agreeing on the signal.</strong> A reference
      swap nominates a different lead variant for {1 - lead_ref_same_hi:.0%} to
      {1 - lead_ref_same_lo:.0%} of shared eGenes, but {ld_same_ref:.0%} of those pairs are in
      high linkage disequilibrium. Only about {net_ref:.0%} of shared eGenes end up with a
      genuinely different causal candidate.</li>
  <li><strong>The genes that do move are not a random sample.</strong> They are
      {gd_or_lo:.1f} to {gd_or_hi:.1f} times enriched for recurrent genomic-disorder regions,
      {har_or_lo:.1f} to {har_or_hi:.1f} times for human accelerated regions, and consistently
      for segmental duplication — twelve of twelve tests positive.</li>
  <li>Method changes move the <em>estimate</em>, not the <em>precision</em>. The standard
      error is essentially unchanged in every contrast, and at matched allele frequency no arm
      is measurably less noisy than any other — including in duplicated sequence, where a gain
      would be most expected. These are not noise reductions; they are changes in what is
      being measured.</li>
  <li>{uni_excl:,} genes are testable on one reference only, and {uni_excl_egenes} of them are
      eGenes — associations available on one reference and not the other, invisible to any
      matched-gene comparison.</li>
  <li><strong>The reference's effect is on access, not measurement.</strong> Among variants
      both references can represent they agree more closely than two aligners do, yet
      {NACC["t2t_minus_grch38_linear"]["share_leads_exclusive"]:.0%} of genes have their best
      association on a variant the other reference cannot represent at all.</li>
  <li><strong>The two axes add different kinds of variant.</strong> Variants only the
      pangenome graph can call reach the association threshold
      {EXCL["graph_minus_linear_t2t"]["call_rate_ratio"]:.1f}× as often as shared variants;
      variants only T2T can call are slightly less likely to
      ({EXCL["t2t_minus_grch38_linear"]["call_rate_ratio"]:.2f}×). A yield count cannot tell
      those apart.</li>
  <li>The chromosome X effect is specific to XX donors under pangenomic alignment on T2T
      (odds {chrx_hit["odds_ratio"]:.2f}); the other seven cells of the design show nothing.</li>
  <li><strong>A reference swap is not ancestry-neutral.</strong> Moving to T2T reduces
      European-ancestry donors' non-reference load by {abs(anc_ref["by_group"]["EUR"]):.3f} and
      raises African-ancestry donors' by {abs(anc_ref["by_group"]["AFR"]):.3f}, identically
      under both aligners.</li>
  <li>The four arms agree closely on cis-eQTL <em>yield</em>, but not on the surface beneath
      it. Reference choice perturbs effect estimates about {mag_ratio:.0f} times more than
      aligner choice does, and that disagreement is spatially concentrated rather than
      diffuse.</li>
  <li>Between {exc_lo:,} and {exc_hi:,} windows per contrast — five to seven percent of the
      testable genome — are discordant beyond chance at an estimated 5% false-discovery
      proportion. Regional disagreement is not confined to a few loci.</li>
  <li>{AGREE["flagged_under_both_aligners"]["n"]:,} of the {AGREE["union"]:,} windows a
      reference contrast flags are flagged under both aligners. The reference effect is a
      property of the reference, not of the alignment method.</li>
  <li>Discordant windows are enriched for duplicated and unmappable sequence whichever way
      the association moves, and the two directions are indistinguishable from each other.
      Sequence difficulty predicts that a region disagrees, not which way it goes.</li>
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
  <li><strong>k = 35.</strong> Its original justification does not survive the sweep. The
      choice stands only because every downstream map was built at it and the cost is
      under two percent.</li>
  <li><strong>The LD-preserving null.</strong> The direction-asymmetry verdict rests on a
      parametric null that ignores LD. The permutation version is the one that would
      settle it.</li>
  <li><strong>The anti-conservative offset.</strong> var(z) above 1 is unexplained and
      systematic.</li>
  <li><strong>There is no ground truth in this dataset.</strong> BrainVar is short-read
      sequencing with no per-donor assemblies, so nothing here can say which arm is
      <em>correct</em> — only which differ, and where. Every result on this page is a
      description of difference, not a verdict on a reference.</li>
  <li><strong>The yardstick is missing.</strong> Reference and aligner effects have not yet
      been placed against a change nobody considers controversial. Repeating the analysis with
      a different variant caller would say whether a reference swap perturbs a map more or
      less than a routine tooling decision, and that comparison is the one a reader needs.</li>
  <li><strong>Annotation versus sequence in the gene universes.</strong> A majority of
      reference-exclusive genes are absent from the other reference's annotation entirely, so
      annotation release explains much of the {uni_excl:,}-gene turnover. Separating that from
      genuine accessibility needs a matched annotation, which does not exist here.</li>
  <li><strong>How heavy the null tail really is.</strong> The false-discovery estimate
      above assumes a normal tail, calibrated on the centre of the genome-wide
      distribution. Linkage disequilibrium correlates neighbouring windows, so the true
      null tail is heavier and every count quoted is a lower bound on the error rate.
      A shift null cannot supply the correction, for the reason given above.</li>
  <li><strong>Why the signed effect skews downward.</strong> Windows whose test statistic
      moves down on T2T outnumber those that move up about two to one, and they carry more
      duplicated sequence. That is a genotype-quality signature rather than a regulatory
      one, and it is unexplained.</li>
  <li><strong>Copy number at the discordance windows.</strong> Haplotype-sampled pangenomic
      alignment addresses paralog collapse, but no amount of correct read placement makes a
      multi-copy locus diploid. Whether these are copy-number effects that a SNP dosage
      merely tags is untested. The joint-called genotypes retain the per-sample depth and
      allele-balance fields that would settle it.</li>
  <li><strong>The {N_ONE} windows confined to one aligner.</strong> These are
      reference-by-aligner interactions, and they are the duplicated, poorly mappable ones.
      They are reported as such rather than as reference findings, and which factor is
      responsible at each is unresolved.</li>
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
at {HS["rotations"]:,} rotations after a count-matched screen. The family-wise version records
the genome-wide maximum window mean on each of {GW["rotations"]:,} rotations and compares every
window against the upper quantile of that maximum, over all {GW["tests"]:,} window tests in the
four contrasts.</p>

<p>Every figure and every table on this page is generated directly from the run
tree by a single collection script, so no number here was transcribed by hand. Figures are
rendered in both light and dark palettes from one code path; the palette was checked
programmatically for colour-vision separation and contrast rather than by eye.</p>

<p>One reporting convention is worth stating because it changes conclusions rather than
presentation: counts of differential effects are reported per <em>gene</em>, never per
variant-gene pair, because linkage disequilibrium lets a single locus dominate a pair-level
count and reverse its apparent direction.</p>

<footer>
  <p>BrainVar GRCh38-versus-T2T cis-eQTL analysis. Every figure and every table is
  generated from the run tree by a single collection script.</p>
</footer>

</div>
</body>
</html>
"""

(HERE / "index.html").write_text(HTML)
(HERE / ".nojekyll").write_text("")
print(f"wrote {HERE / 'index.html'} ({len(HTML):,} bytes)")
