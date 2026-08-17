# Current state

**Updated:** 2026-08-17. Update this file at the end of every work session, in the same
session as the work. A run root without a line here is invisible to the next contributor.

## Next action

**SuSiE fine-mapping is running** on both GPUs — GRCh38 arms on device 0, T2T arms on device 1,
logs at `logs/susie_gpu{0,1}_20260817.log`, run root
`runs/susie_finemapping_v4_k35_20260817`. Shards resume individually, so an interruption is
safe: relaunch the same command and completed shards are skipped.

When it finishes: count credible sets per gene per arm, compare credible-set membership across
arms conditional on a shared gene set, and fold the result into the page. That also unblocks
GWAS colocalization, which needs a trait source in addition.

### After that

- **Add the variant-caller axis at the association level.** The allele-frequency stage gives
  the yardstick cheaply; scoring `hc_minus_dv_*` through the same window and gene pipelines as
  the other axes needs a genotype derivation and association run, which is gated and expensive.
- **Count credible sets per gene.** Locate and validate the SuSiE outputs first — the
  candidates are e36-era and may not correspond to the current v4/k35 arms.
- **Test ancestry-dependence of discordance.**
- **Measure signal in newly accessible sequence**, and **test whether arm-exclusive variants
  carry signal**, stratified by region difficulty.
- **Test whether GWAS colocalization changes.** Last, and dependent on the two above.

## Complete

| Workstream | What it established | Run root |
|---|---|---|
| Concordance baseline | Concordance. Aligner swap r = 0.997, 97.4% identical calls at \|z\|>4, median \|Δz\| = 0.002. Reference swap r = 0.94, 73–75% identical calls, median \|Δz\| = 0.17. 91–95% of genes unmoved. | `gene_discordance_disease_20260816` |
| Discordant regions | Window discordance at 5% FDR: cut z 2.70–2.81, 1,350–1,858 windows per contrast (5–7% of testable genome). | `window_discordance_exceedance_20260816` |
| Discordant regions | Direction: strength balanced (580 up / 595 down); signed effect skewed ~2:1 down. Up and down windows are indistinguishable in sequence context (p = 0.54–0.94); both differ from the rest of the genome. | `directional_window_discordance_20260816` |
| Gene classes | Method-sensitive genes enriched for ClinGen recurrent-CNV regions (OR 2.9–9.1), human accelerated regions (OR 1.7–1.9), segmental duplications (OR 1.7–2.5) — all four contrasts, 12/12 positive. MHC OR 2.5–3.8. | `discordant_gene_classes_20260816` |
| Gene classes | Aligner effect is only ~50% portable across references (gene Jaccard 0.495); reference effect is 92% portable across aligners (0.919). | `gene_discordance_disease_20260816` |
| Calibration | **Negative, and worth keeping:** genome-wide family-wise rotation test returns 0 of 106,342 window tests. A circular shift relocates the signal rather than removing it, so the null max ≈ observed max (ratio 1.02). A shift null tests *position*, not *magnitude*, and cannot certify these windows. | `genomewide_window_discordance_fwer_20260816` |
| chrX by sex | The chrX aligner effect is **XX-specific**: graph−linear·T2T moves 34.6% of chrX genes in XX against 6.8% autosomal (OR 7.27, p = 3.8×10⁻⁸⁵), while the other seven contrast-by-stratum cells show nothing (OR 0.16–0.91). Not power — XY has 133 donors to XX's 92. chrX dosage encoding is identical across all four arms, so differential ploidy handling is excluded. | `chrx_discordance_by_sex_20260816` |
| Effect vs precision | Both axes are **effect-dominated**, not precision-dominated — 93.9% (reference) and 94.9% (aligner) of variants; median se ratio 1.0017 and 1.0000. Neither method change buys precision; both change what is measured. In duplicated sequence the aligner's precision term doubles relative to its effect term while the reference's barely moves. Refutes the hypothesis stated in the plan. | `effect_precision_and_af_20260816` |
| Allele-frequency concordance | Median \|ΔAF\| exactly 0 in all four contrasts. Exactly equal at 70% (reference) and 82–83% (aligner); within 0.01 at 96.2% and 98.4%. The **1.6–3.9% beyond 0.01 is the caution map**, unsigned. Worst chromosomes (chr19, chr22, chr17, chrX, chr6, chr21) match the eQTL-derived hotspots, from an instrument sharing no statistics with them. | `effect_precision_and_af_20260816` |
| Gene universe | The 159-gene net difference between universes hides a **1,175-gene turnover**: 667 genes testable only on GRCh38, 508 only on T2T, and **278 of them are eGenes**. Exclusive genes are more often eGenes than shared ones (26.2% / 21.7% vs 20.8%). T2T-exclusive genes are ~4x as duplicated as shared genes. Majority are absent from the other annotation entirely, so annotation release explains much of the turnover. | `gene_universe_asymmetry_20260816` |
| Top-k concordance | Genome-wide stability does not carry to the head of the ranking. A reference swap changes **19–22 of the top 100** genes (rank correlation 0.83); an aligner swap changes 5 (0.96). Five of the reference-swap changes are genes not testable in the other arm at all. | `topk_rank_concordance_20260816` |
| Lead-variant switching | Among genes that are eGenes in **both** arms, a reference swap changes the lead variant for **41%** and an aligner swap for 22–24%. Median move 11–14 kb; a fifth to a quarter cross the TSS. Whether the two leads tag the same signal is unresolved — LD between them is not yet computed. | `lead_variant_switching_20260816` |
| Lead-switch LD | Most lead switches are relabelling, not new hypotheses. Median r² between the two leads is 0.85 for a reference swap and 0.93 for an aligner swap; 56% and ~70% exceed r² 0.8. Only ~12% (reference) and ~5% (aligner) are effectively independent, so roughly **5% of shared eGenes** get a genuinely different causal candidate under a reference swap. Corrects the emphasis of the distance stage. | `lead_switch_ld_20260816` |
| Ancestry dependence | **The reference swap is not ancestry-neutral.** Moving to T2T changes median alternate-allele load by −0.047 in EUR donors and **+0.017 in AFR** donors (Kruskal p = 1.3×10⁻³⁸), identical under both aligners — European donors come to look more like the reference and African donors less. The aligner axis shows the same test at ~1/50 the magnitude. | `ancestry_dependence_20260816` |
| Arm-exclusive variants | The two axes add different kinds of variant. Variants only the graph can call reach \|z\| ≥ 4 **2.1× more often** than shared variants and sit in 2.4× more duplicated sequence; variants only T2T can call are slightly *less* likely to carry signal (0.92×). Invisible to any yield count. chr1 and chr19. | `arm_exclusive_variants_20260816` |
| Three-axis yardstick | **A reference swap is a smaller perturbation than an aligner swap, and both are under half a caller swap.** Sites differing by more than 0.01 in allele frequency: reference 2.34–2.61%, aligner 2.56–2.91%, caller 5.84–6.50% (2.26× the aligner). Measured on identical donors, contigs and stage. Cross-reference contrasts are restricted to variants both references can represent, which is the point rather than a limitation. | `three_axis_af_yardstick_20260816` |
| Newly accessible signal | **About one gene in nine (11.2%) has its best association at a variant GRCh38 cannot represent**; 4.8% for an aligner swap. With the yardstick and the gene universes this establishes that the reference's distinctive effect is on **access, not measurement** — where it can see the same variant it sees it the same way. | `newly_accessible_signal_20260816` |
| Background noise per arm | **No arm is measurably less noisy.** Median standard error at matched allele frequency is identical to the fourth decimal across all four arms; the largest ratio is T2T carrying 0.8% *larger* errors than GRCh38 in duplicated sequence. Confirms and extends the decomposition: these methods do not buy precision, anywhere. | `background_noise_by_arm_20260816` |

## Known risks

- **chrX must be reported as XX-specific.** The whole-cohort 22.3% figure is a
  dilution of a 34.6% XX-specific effect; stating it as a property of chrX rather
  than of chrX-in-XX would be wrong.
- **Enrichments not yet controlled for expression level.** Duplicated-region genes are often
  lowly expressed; some of the *Characterise discordant regions and gene classes* enrichment could be a power artifact. Open sub-task on *Characterise discordant regions and gene classes*,
  resolved by *Describe background-noise properties per arm*.
- **All FDR estimates are lower bounds.** Normal-null tail assumption; LD correlates
  neighbouring windows and genes.
- **No ground truth exists in this dataset.** No result can be phrased as accuracy. See
  `ANALYSIS_PLAN.md`.
- **A claim in the plan was wrong and has been corrected.** The e36 fine-mapping was
  described as fit on genotypes missing many variants; in fact its adapter imputed
  NaN directly and was never exposed to that defect. It remains unusable because it
  used 36 expression PCs and selected its eGenes from a pre-rebuild map.
- **The ancestry shift is described, not linked to outcomes.** Alternate load measures reference bias; whether it degrades eQTL estimates for AFR donors specifically is untested and is the natural follow-up.
- **`collect_report_data.py` lives outside version control**, in the analysis workspace. It is
  the single point of failure for reproducing the page. Moving it into this repo is a small,
  worthwhile task nobody has done.

## Page status

**Published.** `main` and `hotspot-section` are both at `abe4ad2`; the live site
(`josephlalli.github.io/brainvar-t2t-eqtl-report`) serves it and has been verified to carry the
new title and none of the removed material. Publishing is a fast-forward of `main` to the
working branch and a push — **do not merge or force-push, and ask before publishing.**

Already live: the missing-genotype narrative removed and the page retitled *Mapping Brain
eQTLs on Two Reference Genomes*; the "replication" framing replaced with the factorial
reference-by-aligner reading; the genome-wide FDR and directional sections with their figures.

**Not yet on the page**, though the analysis is complete and in the run tree: the concordance
baseline, the gene-class enrichments, the MHC result, the chrX result, and the aligner-versus-
reference portability asymmetry. These are held back until the framing rewrite below. The chrX
result is now resolved and must be written as XX-specific.

**Framing the page must adopt** (agreed with J.L., 2026-08-16):

1. **Most genetics is unaffected, and say so first.** 91–95% of genes and 97% of association
   calls are untouched. The widely held belief that reference and aligner do not much matter
   is correct for most of the genome, and confirming it is the project's first job.
2. **Then the remainder, not as a footnote.** The minority that moves is enriched 3–9× for
   recurrent genomic-disorder loci and ~1.8× for human accelerated regions, concentrated in
   the MHC, and on T2T reshapes a quarter of chromosome X including a quarter of X-linked
   recessive disease genes.
3. **Then the reason, which is the argument.** These regions are hard to align because they
   vary between people, and that same variability is why they are disease-relevant. Fields
   with the longest history of difficulty — psychiatry, neurodevelopment, immunology — are
   disproportionately exposed to a variable the field treats as inert. Not "their answers are
   wrong", but "they are the fields where representation choice is not free".
4. **The deliverable is a caution map, not a verdict.** Describe differences. Do not rank the
   arms. See `CONVENTIONS.md` → Prose.
