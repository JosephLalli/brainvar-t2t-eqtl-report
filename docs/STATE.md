# Current state

> **STOP — read `REVIEW_FINDINGS_20260819.md` first.** Three adversarial reviews were run on
> 2026-08-19. Two claims currently live on the public page are demonstrably wrong (the
> colocalisation "independent instrument" claim, and "chrX dosages encoded identically"), and a
> confound was found that threatens the reference-axis flagged set and the gene-class
> enrichments computed on it: **expression is quantified per reference, not per arm**, so a
> reference swap changes the phenotype as well as the genotypes. The headline 9× magnitude
> ratio survives, but the crossed association has since **measured** the confound directly: a reference swap turns over 21.7% of eGenes, of which swapping expression alone accounts for 20.1% and swapping genotypes alone only 6.6%. **The reference effect is mostly RNA by volume but genotype by mechanism**: the turnover is largely re-quantification, but the disease-region enrichments computed on the flagged set *survive* when recomputed on the genotype term alone, and two of three are stronger (genomic disorder 2.91× → **4.02×**, segmental duplication 1.65× → **2.40×**). Do not add new analyses before reading it.

**Updated:** 2026-08-20. Update this file at the end of every work session, in the same
session as the work. A run root without a line here is invisible to the next contributor.

## Next action

**Every workstream on the plan is complete or explicitly deferred, and colocalisation is now
Bayesian rather than a membership proxy.** What remains:

- the page. It carries fifteen analyses across eighteen sections and about 14,000 words, and a
  **proposed reordering sits under Page status below, awaiting a decision.** That is a
  presentation change, not another result;
- the variant-caller axis reaches the association level on both references and replicates
  (`caller_axis_association_20260821`, `caller_axis_grch38_20260821`); what remains is that it is
  one caller pair, and that VQSR has no DeepVariant counterpart.

**The open scientific question is now a specific one.** At the colocalisation endpoint the
project's central prediction fails: psychiatric and neurodevelopmental colocalisations change
*less* under a reference swap than other traits (OR 0.70, p = 0.014), and that survives both a
source check and a confidence check. The result is bounded by a selection effect the analysis
cannot see past — colocalisation is only testable where a GWAS has already resolved a clean
signal, so the loci where representation choice might bite hardest are excluded by
construction. Testing that would mean asking whether the discordance windows are depleted of
GWAS signal in the first place, which is a tractable question against the existing run tree
and is the single most informative thing left to do.

### After that

- **Test whether the discordance windows carry GWAS signal at all.** If they are depleted, the
  colocalisation negative is a selection effect rather than a refutation, and that is
  measurable with what is already on disk.
- **Add the variant-caller axis at the association level.** The allele-frequency stage gives
  the yardstick cheaply; scoring `hc_minus_dv_*` through the same window and gene pipelines
  needs a genotype derivation and association run, which is gated and expensive.
- **Finish the background-noise workstream**: Hardy-Weinberg and allele-balance stages need
  passes over the callsets.
- **A SuSiE-based colocalisation**, which would relax the single-causal-variant assumption at
  the cost of needing an LD reference matched to each GWAS. Worth it only if the multi-signal
  loci turn out to matter; the fine-mapping says 11-12% of genes carry more than one credible
  set, so they are not rare.

## Complete

| Workstream | What it established | Run root |
|---|---|---|
| Concordance baseline | Concordance. Aligner swap r = 0.997, 97.4% identical calls at \|z\|>4, median \|Δz\| = 0.002. Reference swap r = 0.94, 73–75% identical calls, median \|Δz\| = 0.17. 91–95% of genes unmoved. | `gene_discordance_disease_20260816` |
| Discordant regions | Window discordance at 5% FDR: cut z 2.70–2.81, 1,350–1,858 windows per contrast (5–7% of testable genome). | `window_discordance_exceedance_20260816` |
| Discordant regions | Direction: strength balanced (580 up / 595 down); signed effect skewed ~2:1 down. Up and down windows are indistinguishable in sequence context (p = 0.54–0.94); both differ from the rest of the genome. | `directional_window_discordance_20260816` |
| Gene classes | Method-sensitive genes enriched for ClinGen recurrent-CNV regions (OR 2.9–9.1), human accelerated regions (OR 1.7–1.9), segmental duplications (OR 1.7–2.5) — all four contrasts, 12/12 positive. MHC OR 2.5–3.8. | `discordant_gene_classes_20260816` |
| Gene classes | Aligner effect is only ~50% portable across references (gene Jaccard 0.495); reference effect is 92% portable across aligners (0.919). | `gene_discordance_disease_20260816` |
| Calibration | **Negative, and worth keeping:** genome-wide family-wise rotation test returns 0 of 106,342 window tests. A circular shift relocates the signal rather than removing it, so the null max ≈ observed max (ratio 1.02). A shift null tests *position*, not *magnitude*, and cannot certify these windows. | `genomewide_window_discordance_fwer_20260816` |
| chrX by sex | The chrX aligner effect is **XX-specific**: graph−linear·T2T moves 34.6% of chrX genes in XX against 6.8% autosomal (OR 7.27, p = 3.8×10⁻⁸⁵), while the other seven contrast-by-stratum cells show nothing (OR 0.16–0.91). Not power — XY has 133 donors to XX's 92. chrX dosage encoding was audited in all four arms: hemizygous calls appear heterozygous in 1.7–2.5% of XY genotypes against ~28.9% in XX, so no arm uses diploid encoding — but the arms differ by 1.47x in that residual and by 1.49x in chrX rows tested, so gross ploidy mishandling is excluded and 'identical' is **not** supported (see REVIEW_FINDINGS_20260819.md). The null cells were also tested one-sided, so the 4–6x depletions reported as 'nothing' are untested. | `chrx_discordance_by_sex_20260816` |
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
| SuSiE fine-mapping | 12,189 gene fine-mappings across four arms in 40 GPU-minutes; 2,540–2,602 genes with credible sets per arm, 1.11–1.12 sets per gene, median set size 7. | `susie_finemapping_v4_k35_20260817` |
| Credible-set comparison | **Fine-mapping survives at the level it reports.** Sets overlap for 96–98% of shared genes (median Jaccard 0.82 reference, 0.97 aligner), but the top variant within them differs **39%** of the time under a reference swap and 22% under an aligner swap. A credible set is robust to method choice; a named causal variant is not. | `credible_set_comparison_20260817` |
| Reference-allele orientation | **"The reference allele" is a property of a choice, not of a variant.** Normalising a GRCh38 variant into the common T2T frame swaps REF and ALT for **30.08%** of them (30.07% in the graph arm, exactly 0% in both T2T arms). Only 0.38% needed indel realignment — the correction the field worries about is small, the one it does not discuss is large. Consequence: GWAS-to-credible-set matching must use **unordered** allele pairs. | `four_arm_variant_identity_comparison_20260811` |
| GWAS variant placement | **The old placement was the bug, not T2T.** The 2022 GWAS VCF placed only 39.9% of catalog rsIDs; LiftoverIndel from the catalog's own GRCh38 coordinates places 97.7%, and an independent native-T2T dbSNP155 placement agrees with it on the exact base for **99.91%** of 331,177 shared rsIDs. Dropout bias toward divergent regions was predicted but is weak (OR 1.09). | `gwas_variant_placement_crosscheck_20260817` |
| GWAS colocalisation (Bayesian) | **Real `coloc.abf`, genome-wide, 40 GWAS in 4 arms.** A reference swap changes 31.8% of calls and an aligner swap 4.0%; among confident calls, 12.2% and 1.6%. Median |ΔPP4| 0.0022 / 0.00005. **All four arms give 926–933 colocalisations** — the direction confound seen under the membership proxy was catalog ascertainment and vanishes with full summary statistics. **The thesis prediction fails here:** brain traits change *least* (OR 0.70, p = 0.014), robust to source (p = 1) and not explained by confidence. Bounded by a selection effect: only loci a GWAS has already resolved are testable. | `gwas_coloc_bayesian_20260817` |
| GWAS colocalisation (membership proxy) | **Superseded.** Credible-set membership overstated instability 2–5× (56.6% vs 31.8% reference; 21.0% vs 4.0% aligner) and is biased upward, not merely noisier. Retained only to quantify how far a variant-overlap proxy departs from the posterior it approximates. | `gwas_coloc_v2_20260817` |

| Selection-effect tests | **The regional claim survives; the trait claim is untestable here.** Colocalisation calls change far more inside FDR-flagged discordance windows (reference 37.1% vs 29.9%, OR 1.38 p=0.028; aligner 8.0% vs 2.5%, **OR 3.44 p=0.00028**) — the genome-wide test was diluting this. The apparent brain-trait deficit (OR 0.70, p=0.014) survives checks on source (p=1), call confidence, and coverage (adjusted OR 0.68, p=0.010) but **not on GWAS signal strength**: brain studies are ~3 orders of magnitude weaker in-window (median 10⁻⁶·² vs 10⁻⁹·⁴), and matched on signal the difference is n.s. (OR 0.80, p=0.15; joint with coverage 0.78, p=0.12). **The disease-gene-set test is retracted** as circular. | `gwas_coloc_bayesian_20260817`, `neuro_gene_sets_20260817` |

| Background noise: HWE and allele balance | **Pangenomic alignment removes the duplicated-sequence penalty; the reference change does not.** On 1,654,230 matched variants, linear arms show elevated collapse signatures in duplicated sequence (excess het 0.283%/0.262%, skewed-het share 0.065/0.056) that the graph arms do not (0.087%/0.082%, 0.039/0.039 — at their ordinary-sequence level). Aligner swap moves excess het −0.195pp in duplicated sequence; reference swap −0.021pp, and its allele-balance effect is n.s. on the graph axis. **First result where an arm is cleaner** — and it separates precision (unchanged) from genotype fidelity (improved). | `background_noise_hwe_ab_20260817` |

| Crossed association | **The reference contrast is dominated by RNA re-quantification.** Running the missing cells of the design (GRCh38 genotypes x T2T expression, and a T2T control through identical code) decomposes a 21.7% eGene turnover into 20.1% from expression alone and only 6.6% from genotypes alone. The control reproduces the published linear_t2t_dv arm exactly (0.9% turnover, r = 0.9999, 100% top-variant agreement), so the machinery is validated. Gene-class enrichments were recomputed on the genotype term — next row. | `crossed_reference_association_20260821` |

| Genotype-term enrichment | **The disease-region enrichments belong to variant representation, not RNA quantification.** Recomputed on the genotype term alone (control vs crossed; expression, covariates and annotation identical): genomic disorder 2.91× → **4.02×** (p = 4e-24), segmental duplication 1.65× → **2.40×** (p = 3e-24), human accelerated 1.84× → 1.54× (p = 0.0072). The genotype term flags a smaller set (579 vs 1,256 genes) that is far more concentrated in hard sequence. Dosage re-signing independently validated against allele frequency (median |Δaf| = 0.0000 for re-signed variants; 0.39 had the sign been inverted). **Bounded:** different z-cuts per side, and the genotype term excludes reference-exclusive variants, so it is a lower bound. | `genotype_term_gene_classes_20260821` |

| Matched gene classes | **The enrichments are not a counting artefact; they are alignment difficulty.** Holding cis-variant count, gene length, gene density and expression level fixed changes nothing -- every estimate holds or strengthens (segdup 2.40x -> 2.94x). Adding mappability attenuates all three: segdup survives on both axes (2.06x / 1.39x), genomic disorder becomes marginal, and **human accelerated regions go null on the genotype term** (1.15x, p = 0.41). For segdup and genomic disorder mappability is the mechanism, not a confounder, so those attenuations are not refutations. | `matched_gene_class_enrichment_20260821` |

| Page truth pass | **Four claims the page made that the run tree does not support are corrected.** The self-contradiction on the yardstick (missing / already measured) and the stale lead-variant LD caveat are replaced by the measured results. The chrX ploidy control now carries its 1.47x and 1.49x between-arm spreads. And the two enrichment terms that survive correction — X-linked recessive inheritance (p = 3e-06) and X-linked inheritance (p = 3e-06), both in graph - linear T2T — are reported rather than omitted; they corroborate the chrX result from an instrument sharing no statistic with it. | `gene_discordance_disease_20260816`, `chrx_discordance_by_sex_20260816`, `lead_switch_ld_20260816` |

| Mappability definition | **The weakest joint in the matched result turns out not to matter.** Mappability recomputed at the gene span, over the nominal cis window, and at the anchored variants themselves. Gene span and variant level agree (r = 0.76) and give the same adjusted odds ratios: segdup 2.06x vs 2.08x, genomic disorder 1.37x vs 1.27x. **The nominal cis window is the poor proxy** (r = 0.44) and under-adjusts, leaving every estimate nearer its crude value -- which is what measurement error in an adjustment covariate predicts. | `window_mappability_20260821` |

| Caller axis at the association level | **A caller swap moves the eQTL map about 1.5x as much as a reference swap does** -- 9.7% eGene turnover against 6.7%, with expression, covariates and annotation held fixed on both, autosomes only. The gate holds: two runs of the same genotypes agree to 1.2%. But the gene-level correlations are near-identical (0.9917 against 0.9927), so the ratio is about threshold crossings rather than a larger perturbation. **The yardstick the page said was missing now exists at the association level, and it does not favour the page's framing.** Bounded by VQSR, which HaplotypeCaller has and DeepVariant has no counterpart to. | `caller_axis_association_20260821` |

| Caller axis on GRCh38 | **The caller effect is a property of the callers, not of T2T.** The caller term replicates across references (9.7% on T2T, 9.4% on GRCh38). On one phenotype side, with expression held fixed by construction, a caller swap moves 1.64x as much as an aligner swap (9.4% vs 5.7%) -- the matched comparison the T2T arm could not make. Association-level order: caller > reference (6.7%) > aligner, matching the allele-frequency yardstick. Correlations span two parts in a thousand across all three, so these differ in threshold crossings rather than in how hard they push. | `caller_axis_grch38_20260821` |

| chrX two-sided, finding 5, Methods | **Three chromosome X cells reported as null are strong depletions.** The published p-values were one-sided in the enrichment direction and seven of eight cells have an odds ratio below one, so they could not return a small p. Retested two-sided: 0.18, 0.16 and 0.24 at p ≤ 7e-09. **Sharpens the result**: chrX is ordinarily more stable than the autosomes under an aligner swap and reverses by 40x in one cell, rather than one spike against silence. Also repaired the manifest that denied its own dosage audit (review finding 5, the last verified one), and gave Methods a provenance table with every parameter, version and joint-calling detail read from an artifact. | `chrx_discordance_by_sex_20260816` |

## Known risks

- **chrX must be reported as XX-specific.** The whole-cohort 22.3% figure is a
  dilution of a 34.6% XX-specific effect; stating it as a property of chrX rather
  than of chrX-in-XX would be wrong.
- **Enrichments are now controlled** for expression level, cis-variant count, gene
  length, gene density and mappability -- see `matched_gene_class_enrichment_20260821`.
  That gap -- mappability over the gene span against a statistic computed over the cis
  window -- was tested in `window_mappability_20260821` and does not matter: measuring
  mappability at the anchored variants themselves reproduces the gene-span answers. The
  reference axis is supported by that agreement rather than by a direct check, because its
  run retained no variant-level output.
- **Superseded note.** Duplicated-region genes are often
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

Publishing is a fast-forward of `main` to the working branch and a push — **do not merge or
force-push, and ask before publishing.**

The branch is at `8f9730c` with uncommitted work in `build_site.py`, `index.html` and both
docs. The page is ~110 KB across **eighteen `h2` sections and about 13,000 words**, and every
analysis on the plan is now written into it — including the concordance baseline, gene-class
enrichments, MHC, chrX-as-XX-specific, ancestry, the genome-wide FDR and directional sections,
credible sets, and colocalisation. *(An earlier version of this note listed several of those as
"not yet on the page"; that was stale and is corrected here.)*

### Proposed structural revision (not yet applied — needs a decision)

A section-by-section read says the page has outgrown its running order rather than its content.
Three specific problems:

1. **`hotspots` is 3,239 words, 3 figures and 6 tables** — a quarter of the page in one
   section. It carries the genome-wide FDR sweep, the directional split and the regional
   detail, which are three separable arguments.
2. **The sex material is fragmented into five consecutive short sections** — `chrx`,
   `interaction`, `stratified`, `contrast`, `direction`, 2,061 words between them — with
   `ancestry` sitting in the middle of the run and breaking it.
3. **The argument peaks in the middle.** `coloc` ("What reaches a paper") is the endpoint the
   whole page builds toward, and nine sections follow it.

The order that matches the agreed framing below would be: setup → *most genetics is
unaffected* (concordance, yardstick) → *but not uniformly* (mechanism, hotspots) → *and the
exceptions are the interesting regions* (classes) → *some genes are not even askable*
(universe) → *and this lands unevenly on people* (ancestry, then the sex cluster together) →
*and here is what it does to a published claim* (coloc) → standing → methods. That moves
`coloc` to the end and `universe` earlier, and groups ancestry with sex.

This is a presentation judgement rather than a factual correction, so it is recorded here
rather than applied.

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
