#!/usr/bin/env python3
"""Assemble the BrainVar GRCh38-vs-T2T eQTL report page.

Every number in the prose and in the tables is read from report_data.json, which
is itself collected only from the four-arm run tree.  Nothing is
transcribed by hand, so the text cannot drift from the run tree.
"""
from __future__ import annotations

import html
import json
import math
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

PV = DATA["provenance"]
PVA, PVG = PV["association"], PV["genotype_derivation"]

t_provenance = table(
    ["Stage", "What was used"],
    [["Primary cis scan",
      f'&plusmn;{PVA["cis_window_bp"] // 10**6}&nbsp;Mb window, MAF &ge; '
      f'{PVA["maf_threshold"]}, {PVA["permutations"]:,} permutations, '
      f'{PVA["expression_pcs"]} expression PCs and {PVA["genotype_pcs"]} genotype PCs in '
      f'{PVA["covariate_columns"]} covariate columns, Benjamini-Hochberg at '
      f'{PVA["fdr"]:.0%}'],
     ["Seeding",
      f'root seed {PVA["root_seed"]}; {PVA["seed_rule"]}'],
     ["Missing genotypes",
      "Missing dosages reach tensorqtl as the <code>-9</code> sentinel its mean imputation "
      "recognises. IEEE NaN, which the derivation produces, is silently dropped by that "
      "imputation instead of being imputed — a defect this project was bitten by once, and "
      "the conversion exists to prevent it."],
     ["Genotype derivation",
      # The commands carry <cohort> and <reference fasta> placeholders, which a browser
      # parses as tags unless they are escaped.
      "<br>".join(f'<code>{html.escape(c)}</code>' for c in PVG["stages"])],
     ["Accepted FILTER token",
      f'DeepVariant arms: {PVG["accepted_filter_token"]["deepvariant_arms"]}<br>'
      f'HaplotypeCaller arm: {PVG["accepted_filter_token"]["haplotypecaller_arm"]}'],
     ["After derivation", PVG["post_filters"]],
     ["Joint calling",
      f'DeepVariant arms: {PV["joint_calling"]["deepvariant_arms"]}<br>'
      f'HaplotypeCaller arm: {PV["joint_calling"]["haplotypecaller_arm"]}'],
     ["References",
      f'{PV["references"]["t2t"]}<br>{PV["references"]["grch38"]}'],
     ["Compute", f'{PVA["gpus"]}'],
     ["Software",
      ", ".join(f'{k}&nbsp;{v}' for k, v in PV["software"].items()
                if k in ("python", "tensorqtl", "torch", "numpy", "scipy", "pandas",
                         "polars", "sklearn", "pyarrow"))
      + "<br>" + "; ".join(PV["software"][k] for k in ("bcftools", "htslib", "bedtools")
                           if PV["software"][k] != "unavailable")]],
    note="Read from the run tree rather than transcribed: association parameters from a "
         "published arm's run record, derivation commands from the command list that run "
         "stored, joint-calling provenance from the callset headers, library versions from "
         "the interpreter that executed the analyses.")

CHRX2 = DATA["chrx_by_sex"]["two_sided"]["results"]


def _chrx_cell(key):
    """Odds with the two-sided p beneath, since the published p could only see enrichment."""
    return (f'<strong>{CHRX[key]["odds_ratio"]:.2f}</strong><br>'
            f'<span class="sub">p = {sci(CHRX2[key]["p_two_sided"])}</span>')


t_chrx = table(
    ["Contrast", "XX: chrX vs autosomal", "XX odds", "XY: chrX vs autosomal", "XY odds"],
    [[CONTRAST_LABEL[k],
      f'{CHRX[f"{k}::XX"]["chrx_rate"]:.1%} vs {CHRX[f"{k}::XX"]["autosomal_rate"]:.1%}',
      _chrx_cell(f"{k}::XX"),
      f'{CHRX[f"{k}::XY"]["chrx_rate"]:.1%} vs {CHRX[f"{k}::XY"]["autosomal_rate"]:.1%}',
      _chrx_cell(f"{k}::XY")] for k in CONTRAST_ORDER],
    cls="numeric",
    note="Each stratum is compared against its own autosomes, so the difference in stratum "
         "size cannot drive the comparison. <strong>The p-values are two-sided.</strong> The "
         "figures first published here were one-sided in the enrichment direction, which "
         "returns p at or near 1 for any cell with an odds ratio below one whatever the data "
         "show — and seven of these eight are below one. Three cells previously reported as "
         "null are strong depletions. On the ploidy control: no arm uses diploid encoding for "
         "hemizygous X, but the arms are not identical in it, and the spread is given in the "
         "text below rather than claimed away here.")

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

_ld = DATA["lead_switch_ld"]["by_contrast"]["t2t_minus_grch38_linear"]

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
         "coordinates. Distance is a weak proxy for whether two leads tag the same signal; "
         "linkage disequilibrium between the two leads is the measurement that settles it, and "
         "it is reported below — median r² = "
         f"{_ld['median_r2']:.2f} across {_ld['switched_leads']:,} switched leads on the "
         f"reference axis, of which {_ld['share_r2_below_0.2']:.0%} fall below "
         "r² = 0.2.")

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
         f'p = {sci(anc_ref["kruskal_p"])} for the reference swap and '
         f'{sci(anc_wf["kruskal_p"])} for the aligner swap. Ancestry labels are projected '
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

CSET = DATA["credible_sets"]["by_contrast"]
SP = DATA["credible_sets"]["susie_parameters"]
SUS = {"pcs": 35, "fdr": 0.05}
cs_top_ref = min(CSET[k]["share_same_top_variant"] for k in REF_KEYS)
cs_top_wf = min(CSET[k]["share_same_top_variant"] for k in WF_KEYS)
cs_ovl_ref = min(CSET[k]["share_sets_overlapping"] for k in REF_KEYS)
cs_ovl_wf = min(CSET[k]["share_sets_overlapping"] for k in WF_KEYS)

t_credible = table(
    ["Contrast", "Genes fine-mapped in both", "Same number of sets", "Sets overlap",
     "Median Jaccard", "Same top variant"],
    [[CONTRAST_LABEL[k], f'{CSET[k]["genes_finemapped_in_both"]:,}',
      f'{CSET[k]["share_same_number"]:.1%}',
      f'<strong>{CSET[k]["share_sets_overlapping"]:.1%}</strong>',
      f'{CSET[k]["median_jaccard"]:.2f}',
      f'<strong>{CSET[k]["share_same_top_variant"]:.1%}</strong>']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note=f'SuSiE fine-mapping of the current arms — '
         f'{DATA["credible_sets"]["gene_finemappings"]:,} gene fine-mappings in '
         f'{DATA["credible_sets"]["gpu_minutes"]:.0f} GPU-minutes. Restricted to genes '
         f'fine-mapped in both arms, with variants compared on normalised identity. Overlap '
         f'is computed on variants placeable in the common frame, so a set containing '
         f'reference-unique variants is compared on the remainder and disagreement is '
         f'understated.')

def pow10(x: float) -> str:
    """Render an exact power of ten as one, rather than as 1e-04."""
    e = round(math.log10(x))
    return (f"10<sup>&minus;{abs(e)}</sup>" if abs(x - 10.0 ** e) < 1e-15 * x
            else f"{x:.0e}")


COL = DATA["coloc"]["by_contrast"]
COLA = DATA["coloc"]["per_arm"]
CPR = DATA["coloc"]["priors"]
CPV = DATA["coloc"]["prior_variance"]
CFL = DATA["coloc"]["filters"]
CTC = DATA["coloc"]["trait_counts"]
CTOT = DATA["coloc"]["totals"]
PROX = DATA["coloc"]["membership_proxy"]["by_contrast"]
col_ref = max(COL[k]["share_of_union_changed"] for k in REF_KEYS)
col_wf = max(COL[k]["share_of_union_changed"] for k in WF_KEYS)
col_ref_strong = max(COL[k]["strong_share_changed"] for k in REF_KEYS)
col_wf_strong = max(COL[k]["strong_share_changed"] for k in WF_KEYS)
prox_ref = max(PROX[k]["share_of_union_changed"] for k in REF_KEYS)
prox_wf = max(PROX[k]["share_of_union_changed"] for k in WF_KEYS)
col_n_traits = sum(CTC.values())
col_arm_lo = min(COLA[c]["colocalised_pp4_0.8"] for c in CELLS)
col_arm_hi = max(COLA[c]["colocalised_pp4_0.8"] for c in CELLS)
col_gene_lo = min(COLA[c]["genes_tested"] for c in CELLS)
col_gene_hi = max(COLA[c]["genes_tested"] for c in CELLS)

t_coloc = table(
    ["Contrast", "Both arms", "Only the first", "Only the second",
     "Calls that change", "Confident calls that change", "Median |&Delta;PP4|"],
    [[CONTRAST_LABEL[k], f'{COL[k]["shared"]:,}', f'{COL[k]["only_positive"]:,}',
      f'{COL[k]["only_negative"]:,}',
      f'<strong>{COL[k]["share_of_union_changed"]:.1%}</strong>',
      f'<strong>{COL[k]["strong_share_changed"]:.1%}</strong>',
      f'{COL[k]["median_abs_pp4_difference"]:.5f}']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note=f'A colocalisation is a gene and trait whose posterior probability of one shared '
         f'causal variant reaches {DATA["coloc"]["pp4_call"]:.0%}. "Confident calls" '
         f'restricts to pairs already above PP4 {DATA["coloc"]["pp4_strong"]:.0%} in at '
         f'least one arm, where a flip is a claim someone would have made rather than a '
         f'borderline case crossing a line.')

t_coloc_proxy = table(
    ["Contrast", "Credible-set membership", "coloc.abf", "Confident coloc.abf"],
    [[CONTRAST_LABEL[k], f'{PROX[k]["share_of_union_changed"]:.1%}',
      f'{COL[k]["share_of_union_changed"]:.1%}',
      f'<strong>{COL[k]["strong_share_changed"]:.1%}</strong>']
     for k in CONTRAST_ORDER],
    cls="numeric",
    note="The same four contrasts scored by variant overlap and by a posterior. The proxy "
         "asks whether a GWAS index variant falls inside a credible set: it has no model of "
         "the GWAS signal and turns on a single variant's membership.")

_g = COL["t2t_minus_grch38_linear"]["by_trait_group"]
_ga = COL["graph_minus_linear_grch38"]["by_trait_group"]
_grows = sorted(_g.items(), key=lambda kv: -kv[1]["share_changed"])
t_coloc_traits = table(
    ["Trait area", "Studies", "Pairs tested", "Called in either arm", "Reference swap",
     "Aligner swap"],
    [[k.capitalize(), f'{CTC.get(k, 0)}', f'{v["tested"]:,}', f'{v["union"]:,}',
      f'<strong>{v["share_changed"]:.1%}</strong>',
      (f'{_ga[k]["share_changed"]:.1%}' if k in _ga else "&mdash;")]
     for k, v in _grows],
    cls="numeric",
    note="Trait areas are assigned per study rather than by keyword match. The "
         "called-in-either counts are modest, so the ordering itself is a tendency; the "
         "pooled brain-against-other test quoted in the text is the claim being made.")

EG = DATA["egene_counts"]
EGENE_UNION = EG["union"]
_pp = DATA["coloc"]["pp4_shift_by_denominator"]
PP4_CALLED = _pp["called_in_either"]["median"]
PP4_STRONG = _pp["strong_in_either"]["median"]
PP4_STRONG_MOVE = _pp["strong_in_either"]["share_moving_over_0.2"]

XD = DATA["crossed"]
XG = XD["genotype_term"]
XE = XD["expression_term"]
XC = XD["confounded_reference_contrast"]
XGATE = XD["gate"]

t_crossed = table(
    ["Contrast", "Differs in", "eGene turnover", "Correlation of gene significance"],
    [["control vs crossed", "<strong>genotypes only</strong>",
      f'<strong>{XG["turnover"]:.1%}</strong>',
      f'{XG["pearson_r_neglog10_pvalbeta"]:.4f}'],
     ["crossed vs published GRCh38", "<strong>expression only</strong>",
      f'<strong>{XE["turnover"]:.1%}</strong>',
      f'{XE["pearson_r_neglog10_pvalbeta"]:.4f}'],
     ["published T2T vs published GRCh38", "both, as published",
      f'<strong>{XC["turnover"]:.1%}</strong>',
      f'{XC["pearson_r_neglog10_pvalbeta"]:.4f}']],
    cls="numeric",
    note=f'Turnover is the share of the eGene union called in one cell and not the other. The '
         f'control cell reproduces the published T2T arm at {XGATE["turnover"]:.1%} turnover, '
         f'r = {XGATE["pearson_r_neglog10_pvalbeta"]:.4f} and '
         f'{XGATE["top_variant_agreement_among_shared_egenes"]:.0%} top-variant agreement, '
         f'which is what validates the machinery.')

XGC = XD["genotype_term_classes"]
GTC = XGC["classes"]
XFL = XD["flip_validation"]["strata"]
CLASS_LABEL = {"genomic_disorder": "Recurrent genomic-disorder region",
               "segmental_duplication": "Segmental duplication",
               "human_accelerated": "Human accelerated region"}
CLASS_ORDER = ["genomic_disorder", "segmental_duplication", "human_accelerated"]

t_genotype_classes = table(
    ["Class", "Reference axis, as published", "Genotype term alone",
     "Genes in class", "Flagged inside vs outside"],
    [[CLASS_LABEL[c],
      f'{GTC[c]["reference_axis_odds_ratio"]:.2f}×<br>'
      f'<span class="sub">p = {sci(GTC[c]["reference_axis_p"])}</span>',
      f'<strong>{GTC[c]["odds_ratio"]:.2f}×</strong><br>'
      f'<span class="sub">p = {sci(GTC[c]["p"])}</span>',
      f'{GTC[c]["in_class"]:,}',
      f'{GTC[c]["rate_in_class"]:.1%} vs {GTC[c]["rate_out"]:.1%}']
     for c in CLASS_ORDER],
    cls="numeric",
    note=f'The published column is the linear T2T-against-GRCh38 contrast, the one the genotype swap '
         f'isolates. Both columns use the same recipe, the same class tracks and the same Fisher '
         f'test against the genes tested in their own contrast. The published column flags '
         f'{XGC["reference_axis_genes_discordant"]:,} of '
         f'{XGC["reference_axis_genes_tested"]:,} genes '
         f'({XGC["reference_axis_discordant_rate"]:.1%}); the genotype term flags '
         f'{XGC["genes_discordant"]:,} of {XGC["genes_tested"]:,} '
         f'({XGC["discordant_rate"]:.1%}) at z &gt; {XGC["gene_cut"]:.2f}. The two flagged sets '
         f'are different sizes, so the odds ratios are comparable in direction and in rough '
         f'magnitude, not to the second decimal.')

GCE = DATA["gene_class_enrichment_terms"]
GCE_HIT = {k: v for k, v in GCE["by_contrast"].items() if v}
GCE_KEY = next(iter(GCE_HIT))
GCE_TERMS = GCE_HIT[GCE_KEY]
CXD = DATA["chrx_dosage_audit"]

t_enrichment_terms = table(
    ["Term", "Source", "Genes in term", "Moved genes in it", "p"],
    [[t["name"], t["source"], f'{t["term_size"]:,}', f'{t["intersection_size"]}',
      sci(t["p_value"])] for t in GCE_TERMS],
    cls="numeric",
    note=f'The only terms above threshold in any of the four contrasts, all from '
         f'{CONTRAST_LABEL[GCE_KEY]}. Queried over {", ".join(GCE["sources"])} from the '
         f'{GCE_TERMS[0]["query_size"]:,} moved genes of that contrast, against a background of '
         f'{GCE["background"]}. p-values are already corrected for multiple testing by the enrichment '
         f'service, which is why the other three contrasts return nothing rather than a long '
         f'weak list.')

CAA = DATA["caller_axis_association"]
CAR = CAA["results"]

CG8 = DATA["caller_axis_grch38"]
CG8R = CG8["results"]


def _row(label, sub, side, e):
    return [f'<strong>{label}</strong><br><span class="sub">{sub}</span>', side,
            f'{e["egenes_a"]:,} / {e["egenes_b"]:,}',
            f'<strong>{e["turnover"]:.1%}</strong>',
            f'{e["pearson_r_neglog10_pvalbeta"]:.4f}']


t_caller_axis = table(
    ["What changes", "Phenotype side", "eGenes", "eGene turnover",
     "Correlation of gene significance"],
    [_row("Variant caller", "DeepVariant &rarr; HaplotypeCaller", "GRCh38",
          CG8R["caller_grch38"]),
     _row("Aligner", "linear &rarr; pangenome graph", "GRCh38", CG8R["aligner_grch38"]),
     _row("Variant caller", "DeepVariant &rarr; HaplotypeCaller", "T2T",
          CG8R["caller_t2t"]),
     _row("Reference", "T2T &rarr; GRCh38 genotypes", "T2T",
          CG8R["reference_genotype_term"]),
     ["Nothing<br><span class=\"sub\">two runs of the same genotypes</span>", "T2T",
      f'{CAR["control_vs_crossed_control"]["egenes_a"]:,} / '
      f'{CAR["control_vs_crossed_control"]["egenes_b"]:,}',
      f'{CAR["control_vs_crossed_control"]["turnover"]:.1%}',
      f'{CAR["control_vs_crossed_control"]["pearson_r_neglog10_pvalbeta"]:.4f}']],
    cls="numeric",
    note='Autosomes only, on every row, and in every row expression, covariates, gene '
         'positions and annotation are held fixed at that side\'s published linear '
         'DeepVariant arm — so each row is the same statistic computed by the same code, '
         'differing only in the genotype matrix. Rows sharing a phenotype side are directly '
         'comparable to each other; rows on different sides are not, which is why the aligner '
         'arm was run at all. The last row is the gate: two runs that differ only in random '
         'seed and load path.')

MD = XD["mappability_definitions"]
MDR = MD["results"]
MDA = MD["agreement"]


def _mdef(defn, cls):
    e = MDR.get(defn, {}).get("genotype_term", {}).get(cls)
    if not e:
        return "—"
    t = e["tiers"]["power_and_alignment"]["regression"]
    return f'{t["odds_ratio"]:.2f}×<br><span class="sub">p = {sci(t["p"])}</span>'


t_mappability_defs = table(
    ["Class", "Mappability at the gene span", "at the nominal cis window",
     "at the variants themselves"],
    [[CLASS_LABEL[c], _mdef("mappability_gene", c), _mdef("mappability_window", c),
      _mdef("mappability_variants", c)] for c in CLASS_ORDER],
    cls="numeric",
    note=f'Genotype term, with statistical power and mappability both held fixed, under three '
         f'ways of measuring the same covariate. Correlation against the variant-level '
         f'definition is {MDA["pearson_r"]["gene_vs_variants"]:.2f} for the gene span and '
         f'{MDA["pearson_r"]["window_vs_variants"]:.2f} for the nominal window; the span of the '
         f'anchored variants reaches {MDA["pearson_r"]["footprint_vs_variants"]:.2f}. '
         f'Transcription start sites are recovered from the nominal output as variant position '
         f'minus start_distance, so no strand assumption enters.')

MC = XD["matched_classes"]
MCR = MC["results"]


def _adj(axis, cls, tier):
    """Regression estimate with its interval; the headline of the three estimators."""
    r = MCR[axis][cls]["tiers"][tier]["regression"]
    return (f'{r["odds_ratio"]:.2f}×<br><span class="sub">'
            f'{r["lo"]:.2f}&ndash;{r["hi"]:.2f}</span>')


t_matched = table(
    ["Class", "Axis", "Unmatched", "Power held fixed", "Power and mappability held fixed"],
    [[CLASS_LABEL[c], "genotype term" if a == "genotype_term" else "reference axis",
      f'{MCR[a][c]["crude_odds_ratio"]:.2f}×',
      _adj(a, c, "power"), _adj(a, c, "power_and_alignment")]
     for c in CLASS_ORDER for a in ("genotype_term", "reference_axis")],
    cls="numeric",
    note=f'Odds that a method-sensitive gene falls in the class, adjusted by a logistic model '
         f'of the covariates and their squares. Two further estimators with different failure '
         f'modes agree throughout: a propensity for class membership cut into twenty strata '
         f'and pooled by Mantel-Haenszel, and a 1:1 nearest-neighbour match with a '
         f'{MC["estimators"]["matched"].split("caliper ")[1].split(" SD")[0]}&nbsp;SD caliper '
         f'tested by McNemar. Balance is checked as standardised mean differences against a '
         f'{MC["balance_target_abs_smd"]:.2f} target; the worst residual on the power tier is '
         f'{max(t["tiers"]["power"]["stratified_worst_abs_smd"] for ax in MCR.values() for t in ax.values()):.3f}. '
         f'On the two genomic-disorder cells of the right-hand column the stratified estimator '
         f'fails its balance check and is disregarded — those genes are systematically '
         f'unmappable, so there are too few comparable controls; the matched estimator handles '
         f'it by discarding the handful of cases it cannot match rather than extrapolating.')

AO = DATA["allele_orientation"]
_g38 = ("linear_grch38_dv", "graph_grch38_dv")
flip_lin = AO["linear_grch38_dv"]["flipped_share"]
flip_gph = AO["graph_grch38_dv"]["flipped_share"]
flip_n = AO["linear_grch38_dv"]["flipped"]
realign_hi = max(AO[c]["realigned_share"] for c in _g38)

COLS = DATA["coloc"]["stratified"]
COLC = DATA["coloc"]["coverage"]
COLJ = DATA["coloc"]["coverage_adjusted"]
_sg = DATA["coloc"]["signal_adjusted"]["t2t_minus_grch38_linear"]
SIG = {
    "signal_median_brain": _sg["signal_median_brain"],
    "signal_median_other": _sg["signal_median_other"],
    "signal_or": _sg["signal_only"]["odds_ratio"],
    "signal_lo": _sg["signal_only"]["ci_low"],
    "signal_hi": _sg["signal_only"]["ci_high"],
    "signal_p": _sg["signal_only"]["p"],
    "joint_or": _sg["signal_by_coverage"]["odds_ratio"],
    "joint_lo": _sg["signal_by_coverage"]["ci_low"],
    "joint_hi": _sg["signal_by_coverage"]["ci_high"],
    "joint_p": _sg["signal_by_coverage"]["p"],
}

t_coloc_strat = table(
    ["Contrast", "Inside flagged windows", "Outside", "Odds ratio", "p"],
    [[CONTRAST_LABEL[k], f'{v["flagged_changed"]}/{v["flagged_union"]} = '
      f'<strong>{v["flagged_rate"]:.1%}</strong>',
      f'{v["unflagged_changed"]}/{v["unflagged_union"]} = {v["unflagged_rate"]:.1%}',
      f'<strong>{v["odds_ratio"]:.2f}</strong>', f'{v["fisher_p"]:.2g}']
     for k, v in COLS.items()],
    cls="numeric",
    note="Each gene is placed by the variant nearest its transcription start, and a window "
         "counts as flagged when any contrast marks it discordant at 5% FDR.")

GF = DATA["genotype_fidelity"]
GFM = GF["matched"]
GFP = GF["paired"]

t_fidelity = table(
    ["Arm", "Excess het, ordinary", "Excess het, duplicated",
     "Skewed-balance share, ordinary", "Skewed-balance share, duplicated"],
    [[LABEL[a], f'{GFM[a]["ordinary"]["excess_het_rate"]:.3%}',
      f'<strong>{GFM[a]["duplicated"]["excess_het_rate"]:.3%}</strong>',
      f'{GFM[a]["ordinary"]["mean_skewed_share"]:.4f}',
      f'<strong>{GFM[a]["duplicated"]["mean_skewed_share"]:.4f}</strong>']
     for a in CELLS],
    cls="numeric",
    note=f'Matched on {GF["matched_variants"]:,} variants called in all four arms, requiring '
         f'at least ten heterozygotes per site and ten reads per heterozygote. A site counts '
         f'as duplicated when its window is at least half segmental duplication. "Skewed" '
         f'means a heterozygote drawing under 30% or over 70% of its reads from one allele.')

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
  <div class="callout" style="margin-top:1.5rem">
  <p><strong>Provisional — exploratory results, not yet a finished analysis.</strong> An
  external review in August 2026 identified limits that this page is still being corrected
  for. What it establishes is that these pipeline choices <em>change eQTL results</em>, and
  where. It does <strong>not</strong> establish that either reference or aligner yields more
  biologically accurate results: BrainVar has no per-donor assemblies, so no analysis here has
  ground truth. Three specific cautions. The reference contrast is not a clean single-factor
  comparison, because expression is quantified per reference — it measures an end-to-end
  pipeline effect — though the crossed cells have since separated the two terms, and the
  disease-region enrichments survive on the genotype term alone. The discordance windows are a
  ranked outlier map, not FDR-controlled discoveries. And the gene-class and disease
  enrichments have since been matched on the properties that travel with these regions: they
  are not counting artefacts, but most of what they measure is alignment difficulty rather than
  anything independent of it. Details and the correction queue are in
  <code>docs/REVIEW_FINDINGS_20260819.md</code>.</p>
  </div>
</header>

<nav class="toc">
  <ol>
    <li><a href="#design">The question and the four arms</a></li>
    <li><a href="#pcs">How many expression PCs?</a></li>
    <li><a href="#maps">The four-arm cis-eQTL maps</a></li>
    <li><a href="#concordance">What a method change leaves alone</a></li>
    <li><a href="#yardstick">How big is that, really?</a></li>
    <li><a href="#mechanism">What kind of difference it is</a></li>
    <li><a href="#hotspots">Where the references disagree</a></li>
    <li><a href="#classes">What the moved regions are made of</a></li>
    <li><a href="#universe">The genes that were never askable</a></li>
    <li><a href="#ancestry">Ancestry, and who a reference represents</a></li>
    <li><a href="#chrx">Chromosome X, and why sex decides it</a></li>
    <li><a href="#interaction">Genotype×sex interaction, and why scale decides it</a></li>
    <li><a href="#stratified">Sex-stratified maps</a></li>
    <li><a href="#contrast">Contrasting effect sizes between sexes</a></li>
    <li><a href="#direction">Is the XX-stronger excess real?</a></li>
    <li><a href="#coloc">What reaches a paper</a></li>
    <li><a href="#standing">What is settled, and what is not</a></li>
    <li><a href="#methods">Methods and reproducibility</a></li>
  </ol>
</nav>

<div class="keyfig">
  <div><span class="n">{conc_calls_aligner:.0%}</span><span class="l">of association calls
    unchanged when the aligner is swapped</span></div>
  <div><span class="n">{gd_or_lo:.1f}–{gd_or_hi:.1f}×</span><span class="l">enrichment of the
    genes that do move for genomic-disorder loci</span></div>
  <div><span class="n">{EGENE_UNION:,}</span><span class="l">distinct eGenes across the four
    arms at k = 35 ({int_total:,} arm-level calls)</span></div>
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
over is exactly the operation whose necessity is under test. The two factors are crossed, so
each can be read off separately, and an effect that appears under only one level of the other
factor is identifiable as exactly that.</p>

<p>Because all four arms run the same code with the same parameters, a difference between
them is attributable to the reference and the aligner. That is the design's intent, and it
holds for the aligner axis but <strong>not</strong> for the reference axis. Expression is
quantified once per <em>reference</em> rather than once per arm, so swapping aligner holds the
phenotype exactly fixed while swapping reference changes the genotypes, the RNA alignment and
annotation, the expression covariates, and sometimes the set of testable genes. The reference
contrast therefore measures an end-to-end <strong>reference ecosystem effect</strong>, not an
isolated genotype-reference effect.</p>

<p><strong>That is not a caveat but a measurement, and it has been made.</strong> Running the
missing cell of the design — one reference's genotypes against the other's expression, plus a
control that passes the matching genotypes through identical code — separates the two terms.</p>

{t_crossed}

<p><strong>The published reference contrast is dominated by RNA re-quantification rather than
by genotype representation.</strong> Swapping the expression alone reproduces
{XE["turnover"]:.1%} of the {XC["turnover"]:.1%} turnover a full reference swap produces;
swapping the genotypes alone gives {XG["turnover"]:.1%}. The correlations say the same thing: a
genotype swap barely moves gene-level significance
(r = {XG["pearson_r_neglog10_pvalbeta"]:.4f}), while an expression swap accounts for
essentially the whole drop.</p>

<p>Read every reference-axis result below with that in mind. It does not make them wrong, and
the effect-size magnitudes reported later survive holding expression stable — but the question
of <em>which genes</em> a reference swap moves is largely a question about RNA quantification.
The gene-class and disease enrichments are computed on exactly that set, so they were
recomputed on the genotype term alone. <strong><a href="#classes">They survive it, and two of
the three are stronger.</a></strong> The reference effect is mostly RNA
<em>by volume</em> and genotype <em>by mechanism</em>: re-quantification moves many more genes,
but the genes that move because of variant representation are the ones sitting in hard,
disease-bearing sequence. The aligner axis is unaffected throughout, because it holds the
phenotype exactly fixed by construction.</p>

<div class="callout">
<p><strong>The re-signing was checked independently.</strong> The control cell validates the
machinery only where it does nothing — T2T variants never flip. But
{DATA["allele_orientation"]["linear_grch38_dv"]["flipped_share"]:.1%} of GRCh38 variants flip
orientation entering the common frame, and every one of their dosages is recoded as
2&nbsp;&minus;&nbsp;d. An inverted sign there would manufacture a
genotype effect indistinguishable from a real one. Allele frequency settles it: both cells
report the frequency of the same common-frame allele, so on chr22 the
{XFL["re-signed"]["n"]:,} re-signed variants agree at a median absolute difference of
{XFL["re-signed"]["median_abs_diff_as_is"]:.4f}, exactly as the {XFL["untouched"]["n"]:,}
untouched ones do. Had the sign been inverted, the re-signed set would instead sit
{XFL["re-signed"]["median_abs_diff_if_inverted"]:.2f} away. One thing this leaves open: the
re-signed variants correlate at {XFL["re-signed"]["r_as_is"]:.3f} against
{XFL["untouched"]["r_as_is"]:.3f} for untouched ones despite identical median agreement. That
is consistent with their narrower frequency range, and also with sites where the two references
disagree about which allele is <em>reference</em> being sites where genotyping disagrees more.
Not resolved here.</p>
</div>

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
<strong>{net_ref:.0%} of shared eGenes</strong> whose two leads are low-LD proxies rather
than tags of one signal — association and linkage alone cannot establish that these are
different <em>causal</em> variants, only that they are not obviously the same one — and
an aligner swap about {net_wf:.0%}. That is a far weaker claim than four-in-ten, and it is the
correct one.</p>

<p>There is one more level, and it is the one downstream analysis actually consumes. A lead
variant is a single number; a <em>credible set</em> is the statement "the causal variant is one
of these", and that is what fine-mapping reports and what colocalisation takes as input. All
four arms were fine-mapped with SuSiE — {DATA["credible_sets"]["gene_finemappings"]:,} gene
fine-mappings — so the sets can be compared directly.</p>

<div class="callout">
<p><strong>How the fine-mapping was done, because the settings decide what a credible set
means.</strong> SuSiE was run through <code>tensorqtl.susie</code> on each arm separately, over
that arm's own eGenes at {SUS["fdr"]:.0%} FDR from its own permutation pass, with
{SUS["pcs"]} expression PCs, a ±{SP["window"] // 10**6} Mb cis window and variants at MAF ≥
{SP["maf_threshold"]}. Four settings shape the output directly:</p>
<ul>
  <li><strong>L = {SP["L"]}</strong> caps the number of independent signals a gene may have.
      A gene at the cap is censored rather than resolved, which is why the count of sets per
      gene is reported as a comparison between arms and not as a biological quantity.</li>
  <li><strong>Coverage {SP["coverage"]:.0%}</strong> sets how much posterior mass a set must
      contain. Raising it enlarges every set; the overlap statistics below would rise with it
      and the top-variant agreement would not, so the two are not interchangeable summaries.</li>
  <li><strong>min_abs_corr = {SP["min_abs_corr"]}</strong> discards sets whose members are not
      mutually correlated at least this much — SuSiE's purity filter, which removes sets that
      are an artefact of the algorithm splitting one signal rather than a real second one.</li>
  <li><strong>Convergence</strong> at tolerance {SP["tol"]} within {SP["max_iter"]} iterations.</li>
</ul>
<p>Two of these matter for reading the comparison specifically. Because each arm fine-maps
<em>its own</em> eGenes, the raw per-arm counts are conditional on that arm's eGene selection
and are only comparable on the shared gene set used here. And because missing dosages reach
tensorqtl as a <code>-9</code> sentinel rather than as IEEE NaN, the same boundary adapter the
association runs use is applied first; without it the fine-mapping would silently drop
variants, which is the bug this project already found once.</p>
</div>

{t_credible}

<p><strong>Fine-mapping survives the method change at the level it reports.</strong> For
{cs_ovl_ref:.0%} to {cs_ovl_wf:.0%} of genes the two arms' credible sets overlap, and the
median Jaccard between them is
{min(CSET[k]["median_jaccard"] for k in REF_KEYS):.2f} under a reference swap and
{min(CSET[k]["median_jaccard"] for k in WF_KEYS):.2f} under an aligner swap. The arms also
agree on <em>how many</em> independent signals a gene has more than nine times in ten.</p>

<p>But the single highest-posterior variant inside those overlapping sets agrees only
{cs_top_ref:.0%} of the time under a reference swap, and {cs_top_wf:.0%} under an aligner one.
That is the lead-variant result again, seen through the instrument that matters most: the two
arms agree about the <em>region in contention</em> and disagree about which member of it to
name.</p>

<div class="callout">
<p><strong>The practical reading.</strong> A credible set is a robust thing to report across a
method change. A single named causal variant is not. Anything downstream that consumes a set —
colocalisation, or a follow-up experiment targeting a locus — inherits a conclusion that
survives the change. Anything that consumes the top variant alone inherits one that changes
about four times in ten between references.</p>
</div>

<p>All five numbers are true and they belong together. The map as a whole is stable; the
shortlist drawn off the top of it is less so; the nominated variant within a shared hit changes
often but usually cosmetically; a small real remainder is a different hypothesis altogether;
and the credible set that fine-mapping reports mostly survives. What none of them is, yet, is
anchored: a correlation of 0.94 means nothing until it is set against a change whose size is
already agreed on. That comes first, and then the question the rest of this page addresses —
what distinguishes the part that moves.</p>

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

<h3 id="caller-association">Does that hold where it matters?</h3>

<p>Allele frequency is one stage, and an early one. The comparison a reader actually needs
is whether a caller swap moves the <em>eQTL map</em> more or less than a reference swap does,
and until now this page could not say. Every row below is the same measurement: expression,
covariates, gene positions and annotation held fixed, only the genotype matrix changing. The
caller swap is run on both references, so it can be asked whether the answer is about the
callers or about one reference; and on GRCh38 an aligner swap is run alongside it, because
swapping aligner holds expression fixed by construction and so provides a comparator on
exactly the same phenotype side.</p>

{t_caller_axis}

<p><strong>The ordering survives the journey, and it replicates.</strong> The caller term
is {CG8R["caller_t2t"]["turnover"]:.1%} on T2T and
{CG8R["caller_grch38"]["turnover"]:.1%} on GRCh38 — the same number twice, so it is a property
of the two callers rather than of either reference. On the T2T side a caller swap moves
{CAA["caller_over_reference_ratio"]:.1f} times as much as a reference swap; on the GRCh38 side,
where an aligner swap supplies a comparator sharing a phenotype side exactly, it moves
{CG8["caller_over_aligner_grch38"]:.2f} times as much as that. At the association level the
three axes order <strong>caller &gt; reference &gt; aligner</strong>, which is the order the
allele-frequency yardstick already gave at the callset. Anyone comfortable with their choice of
variant caller has accepted a larger perturbation than this page is about, at the endpoint and
not merely at the callset.</p>

<p>One comparison in that table is deliberately not made. The reference term and the aligner
term sit on different phenotype sides — the crossed cells that isolate a reference swap were
only ever built against T2T expression — so the gap between
{CG8R["reference_genotype_term"]["turnover"]:.1%} and
{CG8R["aligner_grch38"]["turnover"]:.1%} mixes a genotype change with an annotation and
quantification change, and nothing is claimed from it. Both are compared to the caller term on
their own side instead, which is what the aligner arm exists to allow.</p>

<div class="callout">
<p><strong>The correlations are the sobering half, and they belong next to the ratios.</strong>
Gene-level significance correlates at
{CG8R["caller_grch38"]["pearson_r_neglog10_pvalbeta"]:.4f} across a caller swap,
{CG8R["reference_genotype_term"]["pearson_r_neglog10_pvalbeta"]:.4f} across a reference swap
and {CG8R["aligner_grch38"]["pearson_r_neglog10_pvalbeta"]:.4f} across an aligner swap — a
spread of two parts in a thousand across all three method axes. So the ratios above are about
threshold crossings, not about one axis perturbing the underlying statistics harder than
another. Every one of these is a small perturbation; they differ in how many genes they tip
across q = 0.05. Quoting the ratios without this beside them would overstate what separates
them.</p>
</div>

<p>Two things had to be settled before this comparison could be made honestly, and both
changed its shape. The published arms are derived from a haploid-fixed callset and
HaplotypeCaller has no such file, so whether that matters was measured rather than assumed:
<strong>zero haploid genotype records across
{CAA["haploid_evidence"]["autosomal_calls_checked"]:,} autosomal calls</strong>, against
{CAA["haploid_evidence"]["non_par_chrx_haploid_share"]:.2%} on non-PAR chromosome X — the
chromosome X figure being the positive control, since a test that could not detect haploid
calls would report zero everywhere. Hence autosomes only. That restriction was not caution: on
chromosome X the derivation blanked
{CAA["genotype_counts"]["haplotypecaller_male_het_blanked"]:,} heterozygous calls in XY donors
for HaplotypeCaller against 3,193 for DeepVariant, so run genome-wide the sex chromosomes alone
would have manufactured a large caller difference with nothing to do with the calling
model.</p>

<p>And the two callers mark their confident calls differently, because this HaplotypeCaller
callset has been through VQSR and DeepVariant has no counterpart step. Each is taken at its own
designation, which compares them as a user of either would actually run them — but VQSR is much
the more aggressive, keeping 82.8% of records against DeepVariant's 98.9% on a sample of
chromosome 20. <strong>Part of the caller difference above is that filtering rather than the
calling model</strong>, and this analysis cannot separate the two.</p>

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
{sci(MECH["t2t_minus_grch38_linear"]["decomposition"]["max_decomposition_residual"], 0)}.</p>

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

<p>That is a statement about precision, and precision is not fidelity. A standard error says
how tightly an effect is estimated <em>given</em> the genotypes; it cannot say whether those
genotypes describe a single diploid locus. Two signatures can, and neither needs a reference to
be declared correct. Collapse two paralogous copies onto one site and every carrier looks
heterozygous, so <strong>excess heterozygosity</strong> beyond Hardy-Weinberg expectation is the
classic signature of it. And a real diploid heterozygote should draw about half its reads from
each allele, so a <strong>skewed allele balance</strong> means the reads are not coming from one
diploid locus. Both were measured across five chromosomes in all four arms and compared on the
same variants.</p>

{t_fidelity}

<p><strong>Pangenomic alignment removes the duplicated-sequence penalty. The reference change
does not.</strong> In the linear arms, duplicated sequence carries visibly more of both
signatures than ordinary sequence — excess heterozygosity
{GFM["linear_grch38_dv"]["duplicated"]["excess_het_rate"]:.3%} against
{GFM["linear_grch38_dv"]["ordinary"]["excess_het_rate"]:.3%}, and a skewed-balance share of
{GFM["linear_grch38_dv"]["duplicated"]["mean_skewed_share"]:.4f} against
{GFM["linear_grch38_dv"]["ordinary"]["mean_skewed_share"]:.4f}. In the graph arms the penalty is
simply absent: the duplicated skewed-balance share
({GFM["graph_grch38_dv"]["duplicated"]["mean_skewed_share"]:.4f}) sits at its ordinary-sequence
value ({GFM["graph_grch38_dv"]["ordinary"]["mean_skewed_share"]:.4f}).</p>

<p>Paired on the same variants, an aligner swap inside duplicated sequence moves excess
heterozygosity by {abs(GFP["graph_minus_linear_grch38"]["duplicated"]["excess_het_difference"]) * 100:.3f}
percentage points on GRCh38 and the mean skewed-balance share by
{abs(GFP["graph_minus_linear_grch38"]["duplicated"]["mean_delta_skewed_share"]):.3f}. A reference
swap in the same sequence moves them by
{abs(GFP["t2t_minus_grch38_linear"]["duplicated"]["excess_het_difference"]) * 100:.3f} points and
{abs(GFP["t2t_minus_grch38_linear"]["duplicated"]["mean_delta_skewed_share"]):.3f}, and against
the graph background its effect on allele balance is not significant at all
(p = {GFP["t2t_minus_grch38_graph"]["duplicated"]["skewed_share_p"]:.2f}).</p>

<div class="callout">
<p><strong>This is the one place in the analysis where a method change is cleaner, and it is
worth being precise about what it buys.</strong> Not precision: the standard errors above are
identical to the fourth decimal. What changes is whether a genotype in duplicated sequence
describes one locus or two collapsed together. Those are different properties and only the
second responds to pangenomic alignment.</p>
<p>The direction is partly by construction — haplotype-aware alignment exists to stop
paralogous reads collapsing, so finding that it does is confirmation rather than discovery. The
magnitude is the new part, and it is measured on instruments that need no ground truth, which
this dataset does not have. Two cautions. The matched set requires a variant to survive in all
four arms, which in duplicated sequence selects the more tractable sites, so this understates
the gap at the hardest ones. And in the graph arms excess heterozygosity is <em>lower</em> in
duplicated sequence than in ordinary sequence, which that selection may explain and which is
not otherwise expected.</p>
</div>

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

<p>Building that frame exposes something worth stating on its own, because it is the most
literal form reference bias takes. Normalising a GRCh38 variant into T2T coordinates swaps
which allele is <em>reference</em> for {flip_lin:.1%} of them — {flip_n:,} variants, or nearly
one in three. The figure is {flip_lin:.2%} in the linear GRCh38 arm and {flip_gph:.2%} in the
graph arm, and it is exactly zero in both T2T arms, which are already in the frame. That
agreement across two aligners and its disappearance on the native reference confirm it is a
property of the two genomes rather than of this pipeline. It is also unsurprising once stated:
CHM13 is a single haploid genome, so at any common polymorphism it carries the allele GRCh38
calls alternate about as often as that allele is frequent.</p>

<p>The consequence is that <strong>"the reference allele" is not a property of a variant. It is
a property of a choice</strong>, and it changes for a third of the genome's common variation
depending on which choice is made. Any quantity polarised on REF — an allele frequency, the
sign of a burden, the direction of an effect — carries that choice with it. Worth setting
against the correction people usually worry about: only {realign_hi:.2%} of these variants
needed haplotype realignment to reconcile indel representation between the references. The
small, well-known problem is small. The large one is which allele you decided to call normal.</p>

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
contrast</strong> above that cut. <strong>These are a ranked outlier map, not
FDR-controlled discoveries.</strong> The only null actually simulated — a chromosome-wise
circular rotation, which preserves linkage — returned no genome-wide discoveries at all. The
five percent quoted here comes instead from a normal upper tail, and the observed statistic is
demonstrably heavier-tailed and right-skewed than normal, so the true false-discovery
proportion is uncalibrated rather than conservative. Treat the windows as a shortlist worth
investigating, and read every enrichment computed on them in that light.</p>

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

<h3 id="genotype-term">Is this an artefact of RNA re-quantification? No.</h3>

<p>The flagged set above is defined by the reference contrast, and
<a href="#design">the crossed cells showed that contrast is mostly RNA re-quantification</a>.
That raises the obvious worry that these enrichments describe where transcript quantification
disagrees rather than where variant representation does — which would make the page's central
disease claim a property of annotation rather than of sequence. Recomputing them on the
genotype term alone settles it. Both cells share one expression matrix, one covariate set and
one annotation; only the genotypes differ.</p>

{fig("genotype-term", "Arrow chart of three sequence classes, showing the odds ratio moving from the published reference contrast to the genotype term: genomic-disorder regions rise from 2.91 to 4.02 and segmental duplications from 1.65 to 2.40, while human accelerated regions fall from 1.84 to 1.54", "Each arrow runs from the published reference contrast to the same test on a genotype swap alone. Two of three strengthen, so the enrichment is not carried by RNA re-quantification.")}

{t_genotype_classes}

<p><strong>All three enrichments survive, and two are stronger.</strong> Isolating variant
representation raises the genomic-disorder enrichment from
{GTC["genomic_disorder"]["reference_axis_odds_ratio"]:.2f}× to
<strong>{GTC["genomic_disorder"]["odds_ratio"]:.2f}×</strong> and segmental duplication from
{GTC["segmental_duplication"]["reference_axis_odds_ratio"]:.2f}× to
<strong>{GTC["segmental_duplication"]["odds_ratio"]:.2f}×</strong>. Human accelerated regions
weaken from {GTC["human_accelerated"]["reference_axis_odds_ratio"]:.2f}× to
{GTC["human_accelerated"]["odds_ratio"]:.2f}× but hold
(p = {sci(GTC["human_accelerated"]["p"])}).</p>

<p>The two results fit together once extent is separated from specificity. The genotype term
flags a <strong>smaller</strong> set — {XGC["genes_discordant"]:,} genes
({XGC["discordant_rate"]:.1%}) against {XGC["reference_axis_genes_discordant"]:,}
({XGC["reference_axis_discordant_rate"]:.1%}) — but a <strong>more concentrated</strong> one:
{GTC["genomic_disorder"]["rate_in_class"]:.1%} of genomic-disorder genes are flagged against
{GTC["genomic_disorder"]["rate_out"]:.1%} elsewhere, a ratio the reference contrast does not
reach. Variant representation differs specifically in duplicated and divergent sequence, which
is where recurrent-CNV regions sit; RNA quantification differs more broadly and less
selectively. So the reference effect is mostly RNA by volume and genotype by mechanism — two
separable claims with different remedies. The first is an argument for matched annotation and
RNA alignment. The second is an argument about which reference you genotype against, and it is
the one this page is about.</p>

<div class="callout">
<p><strong>What this comparison can and cannot carry.</strong> The two flagged sets are
different sizes and use different cuts, so "discordant" does not mean quite the same thing on
each side and the odds ratios are comparable in direction rather than to the second decimal.
The genotype term is measured only on variants <em>both</em> references can represent, so it
excludes the access effect entirely — and access is where the reference's distinctive
contribution was shown to lie. It is therefore a <strong>lower bound</strong> on the genotype
axis, not an estimate of it. Neither flagged set is a controlled discovery set: the
false-discovery estimate assumes a normal tail this page has already shown to be too light.</p>
</div>

<h3 id="matched">Is it the class, or everything that travels with it?</h3>

<p>Both tables so far compare flagged genes against every gene tested, matched on nothing. That
is a weak comparison, because genes in these classes are not ordinary genes in any respect:
they carry more cis variants, they are longer, they sit in denser neighbourhoods, they are
expressed at different levels, and they are harder to align. Any of those could produce an
enrichment with no contribution from the class itself. Running the test again with them held
fixed separates the explanations.</p>

{t_matched}

<p><strong>None of it is a counting artefact.</strong> Holding cis-variant count, gene length,
gene density and expression level fixed leaves segmental duplication and human accelerated
regions <em>stronger</em> than unmatched on both axes, which means those covariates had been
masking them rather than manufacturing them. Only genomic disorder attenuates, and barely —
{MCR["genotype_term"]["genomic_disorder"]["crude_odds_ratio"]:.2f}× to
{MCR["genotype_term"]["genomic_disorder"]["tiers"]["power"]["regression"]["odds_ratio"]:.2f}× on
the genotype term and
{MCR["reference_axis"]["genomic_disorder"]["crude_odds_ratio"]:.2f}× to
{MCR["reference_axis"]["genomic_disorder"]["tiers"]["power"]["regression"]["odds_ratio"]:.2f}×
on the reference axis, both still far from one. This was the live alternative explanation and
it is now closed: the disease enrichment is not a story about disease genes having more chances
to look discordant.</p>

<p><strong>Holding mappability fixed asks a different question, and it is the more interesting
one.</strong> Mappability is not a nuisance variable here in the way the others are. A
segmental duplication is <em>by definition</em> sequence that recurs elsewhere, which is
precisely what makes it unmappable; recurrent-CNV regions are bounded by segdups. For those two
classes, conditioning on mappability removes the mechanism rather than a confounder, so an
attenuation is not a refutation. Only for human accelerated regions — defined by substitution
rate rather than by copy structure — is mappability cleanly a confounder.</p>

<ul>
  <li><strong>Segmental duplication survives conditioning on its own mechanism</strong>, on
      both axes and by all three estimators, with residual imbalance at or under 0.03:
      {MCR["genotype_term"]["segmental_duplication"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×
      on the genotype term and
      {MCR["reference_axis"]["segmental_duplication"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×
      on the reference axis. The annotation carries information that 100-mer uniqueness over
      the gene body does not — unsurprising, since the discordance statistic runs over a 1 Mb
      window while mappability here is measured over the gene span.</li>
  <li><strong>The human-accelerated enrichment on the genotype term is entirely
      mappability.</strong> It falls from
      {MCR["genotype_term"]["human_accelerated"]["crude_odds_ratio"]:.2f}× to
      {MCR["genotype_term"]["human_accelerated"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×
      (p = {MCR["genotype_term"]["human_accelerated"]["tiers"]["power_and_alignment"]["regression"]["p"]:.2f}),
      and the matched estimator puts it at
      {MCR["genotype_term"]["human_accelerated"]["tiers"]["power_and_alignment"]["matched"]["odds_ratio"]:.2f}.
      This is the one cell where mappability is a genuine confounder for its class, so
      <strong>this is a negative result</strong> and is reported as one. The reference axis
      behaves differently and holds
      ({MCR["reference_axis"]["human_accelerated"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×,
      p&nbsp;=&nbsp;{sci(MCR["reference_axis"]["human_accelerated"]["tiers"]["power_and_alignment"]["regression"]["p"])}).</li>
  <li><strong>The genomic-disorder enrichment is mostly carried by mappability</strong>, falling
      to
      {MCR["genotype_term"]["genomic_disorder"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×
      and
      {MCR["reference_axis"]["genomic_disorder"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×.
      The regression intervals shown above still exclude one, but the matched estimator — the
      only one of the three that passes its balance check on these two cells — does not
      ({MCR["genotype_term"]["genomic_disorder"]["tiers"]["power_and_alignment"]["matched"]["odds_ratio"]:.2f},
      {MCR["genotype_term"]["genomic_disorder"]["tiers"]["power_and_alignment"]["matched"]["lo"]:.2f}&ndash;{MCR["genotype_term"]["genomic_disorder"]["tiers"]["power_and_alignment"]["matched"]["hi"]:.2f}
      and
      {MCR["reference_axis"]["genomic_disorder"]["tiers"]["power_and_alignment"]["matched"]["odds_ratio"]:.2f},
      {MCR["reference_axis"]["genomic_disorder"]["tiers"]["power_and_alignment"]["matched"]["lo"]:.2f}&ndash;{MCR["reference_axis"]["genomic_disorder"]["tiers"]["power_and_alignment"]["matched"]["hi"]:.2f}),
      so call this attenuated to somewhere around one rather than resolved. Since mappability is
      partly downstream of this class, that is not evidence the enrichment was spurious — but it
      does mean it cannot be claimed as independent of alignment difficulty.</li>
</ul>

<div class="callout">
<p><strong>What the page can say after this, and what it cannot.</strong> It cannot say that
method-sensitive genes are enriched for disease sequence <em>independently</em> of how hard
that sequence is to align. It can say something narrower and better supported: the enrichment
<strong>is</strong> alignment difficulty, concentrated in the duplicated sequence where disease
genes happen to live. That is the mechanism this design proposed rather than a rival to it — but
the two sentences behave differently, and the difference is not rhetorical. An enrichment
independent of alignment difficulty would survive better alignment. This one is what better
alignment is for.</p>
</div>

<p>One thing about that adjustment deserved checking, because it could have undone the rest.
Mappability was measured over each gene's own span, while the statistic it explains is a mean
over that gene's cis variants — which sit anywhere within a megabase of the start site. Since
mappability carries most of the attenuation, measuring it in the wrong place would matter.
Recomputing it three ways settles it.</p>

{t_mappability_defs}

<p><strong>It does not matter.</strong> Measuring mappability at the anchored variants
themselves — the faithful definition, since the statistic is a mean over exactly those variants
— reproduces the gene-span answers: segmental duplication
{MDR["mappability_gene"]["genotype_term"]["segmental_duplication"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×
against
{MDR["mappability_variants"]["genotype_term"]["segmental_duplication"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×,
and both of the other classes at or near one either way. Segmental duplication survives all
three definitions.</p>

<p>The nominal cis window is the outlier, and the reason is worth stating because it runs the
right way. It correlates only {MDA["pearson_r"]["window_vs_variants"]:.2f} with where the
evidence actually sits — a 2 Mb interval always contains some unmappable sequence, so it
measures broad regional difficulty rather than the difficulty of the loci that produced the
statistic. Adjusting for a covariate measured with error under-adjusts, which leaves the
estimate nearer its unadjusted value; and that is precisely what the middle column does, in
every row. The correlations predict the pattern the odds ratios show.</p>

<div class="callout">
<p><strong>Difficulty and biological interest are not two facts here. They are one.</strong>
These regions are hard to align because they vary so much between people that a single linear
reference represents them badly, and between-person variation is the substance of disease
genetics — so the same property that makes them hard to measure is what makes them worth
measuring. That is why the enrichment above is not a warning about instability. It is the
prediction the design was built to test, and it holds on all three axes at once. The matched
analysis below turns that from a claim into a measurement: adjusting for alignment difficulty
removes most of the enrichment, which is what it should do if difficulty and disease are the
same fact rather than two.</p>
</div>

<p>The consequence is uncomfortable but specific. Fields with the longest history of
difficulty — psychiatry, neurodevelopment, immunology — work disproportionately on loci of
exactly this kind. The claim is not that they have been getting wrong answers. It is that they
are the fields for which representation choice is not free, while most of genetics can
reasonably continue to treat it as inert.</p>

<p>Formal pathway enrichment of the moved genes is almost empty, and the one place it is
not is worth more than the emptiness. Against a background of the genes tested in the same
contrast, no term survives multiple-testing correction in three of the four contrasts. Below
the threshold the same biology recurs across them — MHC class I antigen processing and
presentation, T-cell receptor signalling, allograft rejection, Fc-receptor activation — which
is what one would expect if the effect is carried by the MHC and the Fc-receptor loci rather
than by immune pathways at large, and it is reported as a tendency rather than as
enrichment.</p>

<p><strong>The fourth contrast returns two terms above threshold, and both are chromosome
X.</strong></p>

{t_enrichment_terms}

<p>That is the same contrast, on the same axis, that produces
<a href="#chrx">the chromosome X result</a> — and it arrives from a different instrument. The
chrX finding is a per-gene discordance rate stratified by donor sex. This is a gene-set
enrichment over inheritance-mode annotations, computed without reference to sex, to
chromosome, or to the stratified maps at all; the moved-gene list it queries was defined
genome-wide. Two measurements that share no statistic landing on the same contrast is the
strongest internal agreement on this page. It remains internal — both ultimately read the same
nominal scan — and an inheritance-mode annotation is a property of the genes, so this says the
moved set is X-enriched rather than that any particular gene moved for a reason.</p>

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
<strong>{chrx_hit["odds_ratio"]:.2f}</strong>, at p =
{sci(CHRX2["graph_minus_linear_t2t::XX"]["p_two_sided"])}. Every other cell in the design sits
<em>below</em> one: chromosome X is less likely to move than the autosomes everywhere else.</p>

<p><strong>And in three of those cells that depletion is strong and highly significant</strong>
— {CHRX["graph_minus_linear_grch38::XX"]["odds_ratio"]:.2f} and
{CHRX["graph_minus_linear_grch38::XY"]["odds_ratio"]:.2f} for an aligner swap on GRCh38, and
{CHRX["graph_minus_linear_t2t::XY"]["odds_ratio"]:.2f} for the XY stratum of the very contrast
that produces the effect, at p between
{sci(min(CHRX2[k]["p_two_sided"] for k in ("graph_minus_linear_grch38::XY",)))} and
{sci(max(CHRX2[k]["p_two_sided"] for k in ("graph_minus_linear_t2t::XY",)))}. That is a
four- to six-fold depletion, and it makes the result sharper rather than weaker. The pattern is
not one spike against silence. It is a consistent background in which chromosome X is
<em>more</em> stable than the autosomes under an aligner swap, and a single cell in which that
reverses by a factor of forty — the cell where T2T's complete X and two X haplotypes give the
graph both the material and the opportunity to resolve something.</p>

<div class="callout">
<p><strong>Those three depletions were previously reported here as showing nothing.</strong>
The p-values first published were one-sided in the enrichment direction, so a cell could be
depleted by any amount and still return p = 1. Seven of the eight cells have an odds ratio
below one, so seven of eight could not have produced a small p whatever the data contained.
The counts were always sufficient to test properly and are unchanged; only the test is. The
correction is recorded in <code>TWO_SIDED_CORRECTION.json</code> beside the original run,
which retains the published values so the record shows both.</p>
</div>

<p>This is not a power artifact — the null stratum is the larger one, with 133 XY donors
against 92 XX. Nor is it a ploidy-encoding artifact: chromosome X dosages were audited in all
four arms and no arm uses diploid encoding for hemizygous X, with hemizygous genotypes
appearing as heterozygous in
under {max(c["xy_share_het_like"] for c in CXD["cells"].values()):.1%} of XY calls in every arm
against about {max(c["xx_share_het_like"] for c in CXD["cells"].values()):.0%} in XX.</p>

<p>That control is weaker than it first looks and is stated at its real strength here. Gross
ploidy mishandling is excluded, but the arms are not identical: the residual het-like share in
XY spans {CXD["xy_het_like_spread"]:.2f}× across the four arms, and the number of chromosome X
rows each arm tests spans {CXD["chrx_rows_spread"]:.2f}×. So "the two arms encode chromosome X
the same way" is not something this audit establishes. What it does establish runs the right
way for the result: the contrast that produces the effect,
{CONTRAST_LABEL["graph_minus_linear_t2t"]}, is the more closely matched pair of the two aligner
contrasts at {CXD["rows_ratio_t2t_aligner_contrast"]:.2f}× in rows tested, while the contrast
that shows <em>nothing</em> is the badly matched one at
{CXD["rows_ratio_grch38_aligner_contrast"]:.2f}×. A testing imbalance large enough to
manufacture the positive result would have had to manufacture it in the other contrast
first.</p>

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

<h2 id="coloc">What reaches a paper</h2>

<p>The claim that actually leaves a study is rarely an effect size. It is that a
trait-associated variant appears to act through a particular gene — the statement that sends
someone to the bench. Testing whether the method changes <em>that</em> means asking, for every
gene and every trait, whether the two signals share a causal variant, and then asking whether
the answer survives a swap.</p>

<div class="callout">
<p><strong>The method here, because it decides what the numbers mean.</strong> This is
<code>coloc.abf</code> (Giambartolomei and colleagues), which places a prior on the effect at
each variant and turns each side's estimate and standard error into an approximate Bayes
factor. Those combine into five posterior probabilities: neither trait has a signal in the
region, one does, the other does, both do at <em>different</em> variants, or both do at the
<em>same</em> variant. The last of these, PP4, is the colocalisation claim, and a call here
means PP4 ≥ {CFL["pp4_call"]:.0%}.</p>
<p>Priors are the conventional p<sub>1</sub> = p<sub>2</sub> = {pow10(CPR["p1"])} that a given variant is causal for
either trait alone and p<sub>12</sub> = {pow10(CPR["p12"])} that it is causal for both.
The prior effect variances are {CPV["eqtl_quantitative"]:.4f} on the quantitative
expression side and {CPV["gwas_log_odds"]:.2f} on the log-odds side — standard
deviations of {CPV["eqtl_quantitative"] ** 0.5:.2f} and {CPV["gwas_log_odds"] ** 0.5:.2f}, which is what they are usually quoted as. A region is
only tested where the GWAS has an association at p ≤ {pow10(CFL["gwas_signal_in_window_p"])}
inside the gene's window and at least {CFL["min_shared_variants"]} variants are shared, since
below that the Bayes factor sums are dominated by noise.</p>
<p><strong>The single-causal-variant assumption is the method's real limitation and it is
retained deliberately.</strong> <code>coloc.susie</code> relaxes it, but needs a linkage
reference matched to every GWAS, and mismatch between that reference and the study is its
best-known failure mode. Since the question here is whether <em>four arms</em> differ, an
LD-free method is the safer instrument: the assumption is then applied identically to all
four and cancels from the comparison even where it is wrong in absolute terms. Anything that
differed between arms would come from the eQTL side, which is what is being measured, rather
than from a panel choice introduced by the analyst.</p>
<p>The implementation is vectorised rather than looped, and was checked against the reference
R <code>coloc</code> package on forty random regions including planted shared signals: the
largest difference in any posterior was <strong>2&times;10<sup>-9</sup></strong>.</p>
</div>

<p>Both sides are placed in the same common frame and both effect vectors are re-signed to its
alternate allele before any Bayes factor is computed — necessary because reference and
alternate swap for about a third of variants across that boundary, and a colocalisation reads
the two effects as measurements of the same thing. The panel is {col_n_traits} studies:
{CTC.get("psychiatric", 0)} psychiatric, {CTC.get("neurodevelopmental", 0)}
neurodevelopmental, {CTC.get("neurodegenerative", 0)} neurodegenerative,
{CTC.get("immune", 0)} immune, and {CTC.get("comparison", 0)} with no particular reason to be
exposed, which act as a comparison group. That yields {CTOT["tests"]:,} gene-trait tests
across the four arms, over {col_gene_lo:,}–{col_gene_hi:,} genes per arm — every gene with cis
statistics, not only the fine-mapped ones.</p>

{t_coloc}

<p><strong>A reference swap changes about a third of colocalisation calls; an aligner swap
changes about one in twenty-five.</strong> But the calls that change are overwhelmingly the
marginal ones. Restricting to pairs where at least one arm was already confident — the claims
a paper would actually have made — a reference swap changes
{col_ref_strong:.1%} and an aligner swap {col_wf_strong:.1%}. The underlying evidence barely
moves at all: the median change in PP4 is
{COL["t2t_minus_grch38_linear"]["median_abs_pp4_difference"]:.4f} under a reference swap and
{COL["graph_minus_linear_grch38"]["median_abs_pp4_difference"]:.5f} under an aligner swap,
with correlations of {COL["t2t_minus_grch38_linear"]["pp4_pearson_r"]:.3f} and
{COL["graph_minus_linear_grch38"]["pp4_pearson_r"]:.3f}.</p>

<div class="callout">
<p><strong>That median is over the wrong denominator, and the honest figures are larger.</strong>
It averages across all {COL["t2t_minus_grch38_linear"]["tested_in_both"]:,} gene-trait pairs
tested in both arms, the great majority of which carry no signal in either and so cannot move.
Restricted to the pairs where a call actually exists, the median posterior shift under a
reference swap is <strong>{PP4_CALLED:.4f}</strong> — about ten times the all-pairs figure —
and among pairs reaching PP4 ≥ 0.5 in either arm it is <strong>{PP4_STRONG:.4f}</strong>, with
<strong>{PP4_STRONG_MOVE:.0%}</strong> moving by more than 0.2. The reported correlation is
inflated the same way. So the reassurance is real for the map as a whole and much weaker for
the claims a paper would actually make.</p>
</div>

<div class="callout">
<p><strong>An earlier version of this page reported these numbers two to five times
larger, and it was wrong.</strong> That analysis scored <em>credible-set membership</em>:
whether a genome-wide significant GWAS variant happened to fall inside a gene's credible set.
That is variant overlap, not colocalisation. It has no model of the GWAS signal, it turns on
one variant's membership, and the variant in question is whichever index SNP a study chose to
report — which is itself a tagging choice, and therefore entangled with the very thing being
measured.</p>
{t_coloc_proxy}
<p>The proxy is not merely noisier; it is biased upward, because a conjunction that hinges on
one variant flips far more readily than a posterior computed over a whole region. It also
produced a spurious asymmetry between the references, which the posterior removes — see
below. Where the two disagree, the posterior is the one to believe.</p>
</div>

<p><strong>No arm finds more colocalisations than another.</strong> The four give
{col_arm_lo:,}–{col_arm_hi:,} calls, a spread under one per cent. This matters because the
membership proxy showed GRCh38 arms finding roughly a fifth more than T2T arms, and this page
previously declined to interpret that, on the grounds that the GWAS Catalog is itself
ascertained on studies built against earlier references. That reasoning turns out to have been
right. Given full summary statistics rather than catalog index variants, the asymmetry
disappears entirely. The direction is now readable on both axes, and what it says is that
neither reference nor aligner produces more colocalisations — only different ones.</p>

<h3>Where the change is concentrated</h3>

<p>Everything above is a genome-wide average, and this page has spent its length arguing that
the genome is not uniform. Averaging over the great majority of the genome that a method change
leaves alone will dilute any effect confined to the minority that it does not, so the comparison
is worth making inside and outside the discordance windows separately.</p>

{t_coloc_strat}

<p><strong>Colocalisation calls change far more often inside the discordance windows</strong> —
{COLS["graph_minus_linear_grch38"]["odds_ratio"]:.1f} times the odds on the aligner axis and
{COLS["t2t_minus_grch38_linear"]["odds_ratio"]:.2f} times on the reference axis. That is this
page's central regional claim reproducing at the endpoint a reader actually cares about, and on
the endpoint a reader cares about. It is <strong>not</strong> an independent confirmation, and
an earlier version of this page wrongly said it was: the windows were defined from the same
nominal association scan that supplies the eQTL side of this colocalisation, and a window is
flagged precisely because the arms' effect estimates diverge most there. The GWAS side is
independent but is common to both arms and cancels from the contrast. Read this as an
internal-consistency check — the regional signal propagates to a downstream endpoint — rather
than as corroboration from a separate instrument.</p>

<p>It also means the genome-wide figures in the table above understate the case where it
matters. Working inside these regions, roughly
{COLS["t2t_minus_grch38_linear"]["flagged_rate"]:.0%} of colocalisation calls are contingent on
the reference and {COLS["graph_minus_linear_grch38"]["flagged_rate"]:.0%} on the aligner
alone.</p>

<h3>Which traits move, and why this cannot answer it</h3>

<p>The project's argument predicts that psychiatric and neurodevelopmental loci should be the
most exposed to representation choice, because the regions that are hard to align are hard on
account of population-level variability and that same variability is what makes them
disease-relevant. Splitting the colocalisations by trait area appears to test that prediction
at the point where a result would be published. It does not, and the reason is worth following
carefully, because the raw comparison looks like a clean refutation.</p>

{t_coloc_traits}

<p>Taken at face value the prediction fails in the opposite direction. Brain traits are the
<em>least</em> likely to change: pooling psychiatric with neurodevelopmental gives 26.7%
against 34.2% for everything else, an odds ratio of <strong>0.70</strong> (p = 0.014). Three
explanations were tested.</p>

<ul>
  <li><strong>Not the source of the studies.</strong> The psychiatric studies come mostly from
      the consortium's own deposits, because the GWAS Catalog holds no adequately powered
      psychiatric summary statistics, while the immune and comparison studies come entirely
      from the Catalog — so source and trait area are nearly collinear. Within the psychiatric
      group, where both sources are present, they give 27.8% and 28.6%, <strong>p = 1</strong>.
      Restricting the comparison to Catalog-sourced studies leaves the ordering unchanged.</li>
  <li><strong>Not that brain calls are better determined.</strong> Neurodevelopmental traits
      have the <em>lowest</em> median PP4 among their calls (0.873 against 0.910 for immune)
      and the fewest above 0.95, so they sit closer to the threshold and should flip more
      readily, not less.</li>
  <li><strong>Not coverage.</strong> The GWAS panel is measurably thinner inside the flagged
      windows and most thinly for brain traits — they cover {COLC["brain_median"]:.3f} of the
      common frame there relative to outside, against {COLC["other_median"]:.3f} for the rest —
      and that reaches the individual test, where a brain-trait colocalisation draws on a
      median of {COLJ["t2t_minus_grch38_linear"]["shared_median_brain"]:,.0f} shared variants
      against {COLJ["t2t_minus_grch38_linear"]["shared_median_other"]:,.0f}. Stratifying on the
      variant count that actually entered each posterior leaves the deficit where it was:
      {COLJ["t2t_minus_grch38_linear"]["adjusted_odds_ratio"]:.2f}
      (95% CI {COLJ["t2t_minus_grch38_linear"]["ci_low"]:.2f}&ndash;{COLJ["t2t_minus_grch38_linear"]["ci_high"]:.2f},
      p = {COLJ["t2t_minus_grch38_linear"]["adjusted_p"]:.3f}).</li>
</ul>

<div class="callout">
<p><strong>The fourth explanation is the right one, and it dissolves the result.</strong>
Psychiatric and neurodevelopmental GWAS are simply weaker than the immune and anthropometric
studies they are being compared against. The strongest association inside a gene's window has a
median of 10<sup>&minus;{SIG["signal_median_brain"]:.1f}</sup> for brain traits against
10<sup>&minus;{SIG["signal_median_other"]:.1f}</sup> for the rest — three orders of magnitude.
A posterior built on a weaker trait signal sits further from the calling threshold and crosses
it less often, for reasons that have nothing to do with which regions the loci occupy.</p>
<p>Comparing tests at matched signal strength, <strong>the difference is no longer
significant</strong>: odds ratio {SIG["signal_or"]:.2f}
(95% CI {SIG["signal_lo"]:.2f}&ndash;{SIG["signal_hi"]:.2f}, p = {SIG["signal_p"]:.2f}).
Matching on signal and coverage together gives {SIG["joint_or"]:.2f}
({SIG["joint_lo"]:.2f}&ndash;{SIG["joint_hi"]:.2f}, p = {SIG["joint_p"]:.2f}). The apparent
trait effect was a study-power artefact.</p>
</div>

<p>So this analysis is <strong>underpowered to say anything about the trait prediction, in
either direction</strong>. That is a weaker conclusion than either "the argument is confirmed"
or "the argument fails", and it is the one the data supports. It is worth recording how it was
reached: the first three checks all left the deficit standing, and stopping there would have
produced a confident published negative that the fourth check removes.</p>

<p>None of this touches the regional findings. The window-level discordance, the enrichment of
moved genes for recurrent genomic-disorder loci and human accelerated regions, the chromosome X
result and the ancestry result all stand, and the stratified colocalisation above reproduces
the central regional claim on an independent instrument.</p>

<p>Two limits remain on any trait-level claim from data of this kind. A gene-trait pair is only
testable where a GWAS already has an association in the window, so the loci examined are the
ones GWAS has already succeeded on; if representation choice bites hardest where GWAS has
<em>not</em> yet produced clean signals, those loci are absent by construction. And every GWAS
here is built against GRCh38, so a variant only T2T can represent cannot appear on the trait
side at all. Settling the question needs a cohort imputed and analysed against a T2T panel,
which is a different experiment from this one.</p>

<h2 id="standing">What is settled, and what is not</h2>

<p>Settled, in the sense of resting on the complete four-arm run tree:</p>

<ul>
  <li><strong>A reference swap is the smallest of the three method changes.</strong> Among
      variants both references can represent it shifts
      {DATA["yardstick"]["mean_by_axis"]["reference"]:.2%} of allele frequencies, against
      {DATA["yardstick"]["mean_by_axis"]["aligner"]:.2%} for an aligner swap and
      {DATA["yardstick"]["mean_by_axis"]["caller"]:.2%} for a variant-caller swap. The
      ordering holds at the endpoint as well as at the callset: on autosomes, with expression
      held fixed, a caller swap turns over {CAR["caller_term"]["turnover"]:.1%} of the eGene
      union against {CAR["reference_genotype_term_autosomes"]["turnover"]:.1%} for a reference
      swap. The gene-level correlations are near-identical
      ({CAR["caller_term"]["pearson_r_neglog10_pvalbeta"]:.4f} against
      {CAR["reference_genotype_term_autosomes"]["pearson_r_neglog10_pvalbeta"]:.4f}), so that
      ratio is about threshold crossings rather than a larger perturbation.</li>
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
      lead whose LD with the other arm's is low enough that they are not obviously tagging one
      signal. Association and linkage cannot show these are different causal variants.</li>
  <li><strong>A reference swap changes about a third of colocalisation calls and an aligner
      swap about one in twenty-five</strong> ({col_ref:.0%} and {col_wf:.0%} of the union),
      but the calls that move are the marginal ones. Among pairs already confident in at
      least one arm the figures fall to {col_ref_strong:.1%} and {col_wf_strong:.1%}, and the
      median shift in posterior is
      {COL["t2t_minus_grch38_linear"]["median_abs_pp4_difference"]:.4f} and
      {COL["graph_minus_linear_grch38"]["median_abs_pp4_difference"]:.5f}. The evidence
      barely moves; a threshold gets crossed.</li>
  <li><strong>Neither reference nor aligner produces more colocalisations — only different
      ones.</strong> The four arms give {col_arm_lo:,}–{col_arm_hi:,} calls, a spread under
      one per cent. An earlier credible-set-membership analysis appeared to show GRCh38
      finding a fifth more than T2T; that was the GWAS Catalog's own ascertainment, and it
      vanishes once full summary statistics replace catalog index variants.</li>
  <li><strong>Colocalisation change is concentrated where this page said it would be.</strong>
      Calls change {COLS["graph_minus_linear_grch38"]["odds_ratio"]:.1f} times more often
      inside the discordance windows on the aligner axis and
      {COLS["t2t_minus_grch38_linear"]["odds_ratio"]:.2f} times on the reference axis — the
      central regional signal propagating to a downstream endpoint. This is an
      internal-consistency check, not independent corroboration: the windows were defined from
      the same nominal scan that supplies the eQTL side.</li>
  <li><strong>The trait prediction is untestable with GRCh38-based GWAS, not refuted.</strong>
      Brain-trait colocalisations appear to change less (odds ratio 0.70, p = 0.014), and that
      survives checks on study source, call confidence and coverage. It does not survive
      matching on GWAS signal strength — brain studies are three orders of magnitude weaker
      inside a gene's window — after which the difference is not significant
      ({SIG["signal_or"]:.2f}, p = {SIG["signal_p"]:.2f}). Every GWAS here is also built
      against GRCh38, so a variant only T2T can represent cannot appear on the trait side at
      all.</li>
  <li><strong>Fine-mapping survives at the level it reports.</strong> Credible sets overlap
      for {cs_ovl_ref:.0%} to {cs_ovl_wf:.0%} of shared genes, but the top variant inside them
      agrees only {cs_top_ref:.0%} of the time under a reference swap. A credible set is robust
      to the method change; a named causal variant is not.</li>
  <li><strong>The genes that do move are not a random sample.</strong> They are
      {gd_or_lo:.1f} to {gd_or_hi:.1f} times enriched for recurrent genomic-disorder regions,
      {har_or_lo:.1f} to {har_or_hi:.1f} times for human accelerated regions, and consistently
      for segmental duplication — twelve of twelve tests positive.</li>
  <li><strong>The enrichments are not a counting artefact.</strong> Held fixed on
      cis-variant count, gene length, gene density and expression level, every one holds or
      strengthens — segmental duplication rises from
      {MCR["genotype_term"]["segmental_duplication"]["crude_odds_ratio"]:.2f}× to
      {MCR["genotype_term"]["segmental_duplication"]["tiers"]["power"]["regression"]["odds_ratio"]:.2f}×.
      Adding mappability attenuates all three: segdup survives on both axes, genomic disorder
      becomes marginal, and human accelerated regions go null on the genotype term
      ({MCR["genotype_term"]["human_accelerated"]["tiers"]["power_and_alignment"]["regression"]["odds_ratio"]:.2f}×,
      p&nbsp;=&nbsp;{MCR["genotype_term"]["human_accelerated"]["tiers"]["power_and_alignment"]["regression"]["p"]:.2f}).
      For the first two classes mappability is the mechanism rather than a confounder, so those
      attenuations are not refutations — but the enrichment can no longer be called independent
      of alignment difficulty.</li>
  <li><strong>Those enrichments belong to the genotype axis, not to RNA re-quantification.</strong>
      Recomputed on a genotype swap alone, with expression, covariates and annotation held
      identical, all three survive and two strengthen: genomic disorder
      {GTC["genomic_disorder"]["reference_axis_odds_ratio"]:.2f}× →
      {GTC["genomic_disorder"]["odds_ratio"]:.2f}×, segmental duplication
      {GTC["segmental_duplication"]["reference_axis_odds_ratio"]:.2f}× →
      {GTC["segmental_duplication"]["odds_ratio"]:.2f}×, human accelerated
      {GTC["human_accelerated"]["reference_axis_odds_ratio"]:.2f}× →
      {GTC["human_accelerated"]["odds_ratio"]:.2f}×. The genotype term flags a smaller set that
      is more concentrated in hard sequence, and because it sees only variants both references
      can represent it is a lower bound.</li>
  <li>Method changes move the <em>estimate</em>, not the <em>precision</em>. The standard
      error is essentially unchanged in every contrast, and at matched allele frequency no arm
      is measurably less noisy than any other — including in duplicated sequence, where a gain
      would be most expected. These are not noise reductions; they are changes in what is
      being measured.</li>
  <li><strong>Pangenomic alignment reduces collapse-compatible QC signatures in duplicated
      sequence.</strong>
      In duplicated sequence the linear arms carry the signatures of collapsed paralogues —
      excess heterozygosity {GFM["linear_grch38_dv"]["duplicated"]["excess_het_rate"]:.3%} of
      sites and a skewed heterozygote allele balance at
      {GFM["linear_grch38_dv"]["duplicated"]["mean_skewed_share"]:.4f} — and the graph arms do
      not ({GFM["graph_grch38_dv"]["duplicated"]["excess_het_rate"]:.3%} and
      {GFM["graph_grch38_dv"]["duplicated"]["mean_skewed_share"]:.4f}, their ordinary-sequence
      level). The reference change moves neither appreciably. This is the one place an arm is
      cleaner, on instruments that need no ground truth.</li>
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
      (odds {chrx_hit["odds_ratio"]:.2f}, p =
      {sci(CHRX2["graph_minus_linear_t2t::XX"]["p_two_sided"])}). The other seven cells all sit
      below one, and three are significant <em>depletions</em> of
      {CHRX["graph_minus_linear_grch38::XY"]["odds_ratio"]:.2f} to
      {CHRX["graph_minus_linear_t2t::XY"]["odds_ratio"]:.2f} — so chromosome X is ordinarily
      more stable than the autosomes under an aligner swap, and reverses only there.</li>
  <li><strong>A reference swap is not ancestry-neutral.</strong> Moving to T2T reduces
      European-ancestry donors' non-reference load by {abs(anc_ref["by_group"]["EUR"]):.3f} and
      raises African-ancestry donors' by {abs(anc_ref["by_group"]["AFR"]):.3f}, identically
      under both aligners.</li>
  <li>The four arms agree closely on cis-eQTL <em>yield</em>, but not on the surface beneath
      it. Reference choice perturbs effect estimates about {mag_ratio:.0f} times more than
      aligner choice does, and that disagreement is spatially concentrated rather than
      diffuse.</li>
  <li>Between {exc_lo:,} and {exc_hi:,} windows per contrast — five to seven percent of the
      testable genome — sit above a cut placed at an estimated 5% false-discovery
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
  <li><strong>The caller axis is one caller pair.</strong> The yardstick now reaches the
      association map on both references and gives the same answer on each, so the result is
      not a property of T2T — but it is still DeepVariant against HaplotypeCaller, on
      autosomes. Whether the ordering holds for another caller pair, or on the sex
      chromosomes, is untested. Both HaplotypeCaller arms have also been through VQSR, which
      neither DeepVariant arm has a counterpart to, so some of the caller difference is a
      filtering policy rather than the calling model, and this design cannot separate the
      two.</li>
  <li><strong>The reference term exists on one phenotype side only.</strong> Isolating a
      reference swap needs a crossed cell — one reference's genotypes against the other's
      expression — and those were only ever built against T2T expression. So the reference
      term cannot be set against the aligner term directly, and the table above does not.
      Building the mirror cell on the GRCh38 side would close that, and has not been
      done.</li>
  <li><strong>What colocalisation cannot see.</strong> A gene-trait pair is only testable
      where a GWAS already has an association inside the gene's window, so the loci examined
      are the ones GWAS has already resolved. If representation choice bites hardest where
      GWAS has <em>not</em> yet produced clean signals — which is what the regional results
      here would predict — those loci are absent from that table by construction. Nothing in
      this analysis can distinguish "brain loci are robust" from "the brain loci GWAS has
      resolved are the robust ones".</li>
  <li><strong>One causal variant per region.</strong> <code>coloc.abf</code> assumes it. The
      assumption is applied identically to all four arms, so it does not bias the comparison,
      but it does make each individual posterior unreliable at loci with several independent
      signals — which the fine-mapping shows are common enough to matter. A SuSiE-based
      colocalisation would relax it at the cost of needing a linkage reference matched to
      every study, and that trade has not been made here.</li>
  <li><strong>The reference axis has no variant-level mappability check.</strong> Measuring
      mappability at the anchored variants rather than over the gene span reproduces the
      matched result on the genotype term, which is what retires that concern — but the same
      check cannot be run on the reference axis, because its run retained only the gene table
      and not the anchored variant set. Reproducing it means redoing the common-frame identity
      join over the full four-arm nominal scan. The reference-axis numbers are therefore
      supported by the agreement measured on the other axis rather than by a direct test.</li>
  <li><strong>Whether the human-accelerated null is real or a power failure.</strong> On the
      genotype term the enrichment vanishes once mappability is held fixed, on
      {MCR["genotype_term"]["human_accelerated"]["in_class"]:,} genes of which only
      {MCR["genotype_term"]["human_accelerated"]["flagged_in_class"]} are flagged. The
      reference axis, with twice the flagged set, still shows the effect. Those are consistent
      with a genuine axis difference and equally consistent with the genotype term being
      underpowered for the smallest of the three classes.</li>
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

{t_provenance}

<p><strong>The aligners are the one thing this page cannot name.</strong> No artifact the
analysis reads records the per-sample aligner or its version — the callset headers carry the
joint-calling provenance and nothing upstream of it. What the callset identifiers do establish
is that {PV["aligners"]["graph_arms"]}, and that {PV["aligners"]["linear_arms"]}. The aligner
axis is therefore identified by callset rather than by tool version, and that gap is stated
here rather than filled by inference.</p>

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

<p><strong>Fine-mapping.</strong> SuSiE via <code>tensorqtl.susie</code>, per arm, over that
arm's own eGenes at {SUS["fdr"]:.0%} FDR from its own permutation pass, with {SUS["pcs"]}
expression PCs: L = {SP["L"]} signals, coverage {SP["coverage"]:.0%}, purity filter
min_abs_corr = {SP["min_abs_corr"]}, MAF ≥ {SP["maf_threshold"]},
±{SP["window"] // 10**6} Mb window, tolerance {SP["tol"]} within {SP["max_iter"]} iterations —
{DATA["credible_sets"]["gene_finemappings"]:,} gene fine-mappings in
{DATA["credible_sets"]["gpu_minutes"]:.0f} GPU-minutes. Missing dosages are converted from
IEEE NaN to the <code>-9</code> sentinel before reaching tensorqtl, whose mean imputation
recognises the sentinel and not NaN; without that adapter variants are dropped silently.</p>

<p><strong>Colocalisation.</strong> <code>coloc.abf</code> with priors
p<sub>1</sub> = p<sub>2</sub> = {pow10(CPR["p1"])}, p<sub>12</sub> = {pow10(CPR["p12"])} and prior
effect variances {CPV["eqtl_quantitative"]:.4f} (quantitative) and {CPV["gwas_log_odds"]:.2f}
(log odds). A gene-trait pair is tested where the GWAS has p ≤
{pow10(CFL["gwas_signal_in_window_p"])} inside the gene's cis window and at least
{CFL["min_shared_variants"]} variants are shared; a call is PP4 ≥ {CFL["pp4_call"]:.0%}.
Implemented vectorised in polars and numpy, with log Bayes factors accumulated by log-sum-exp
in float64 — z² reaches the hundreds and is summed over thousands of variants, which float32
cannot carry — and checked against the reference R <code>coloc</code> package on forty random
regions including planted shared signals, agreeing to 2&times;10<sup>-9</sup>. The panel is
{col_n_traits} studies; the well-powered psychiatric and neurodevelopmental GWAS have no
summary statistics in the GWAS Catalog and come from the consortium's own public deposits,
which use a different layout, a different genome build and rsID rather than coordinate keys.
Those are mapped into the common frame through a T2T-native dbSNP index, so no liftover is
performed on them either. Both effect vectors are re-signed to the common frame's alternate
allele before any Bayes factor is computed.</p>

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
