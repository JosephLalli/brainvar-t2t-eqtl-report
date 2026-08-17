# Current state

**Updated:** 2026-08-16. Update this file at the end of every work session, in the same
session as the work. A run root without a line here is invisible to the next contributor.

## Next action

**Separate the chrX effect by sex.** Do this before anything else. graph-vs-linear on T2T moves 22.3% of
chrX genes (OR 4.47, p = 4×10⁻³⁰) and no other contrast does anything on chrX. It is
unexplained and suspect. XY donors are hemizygous, the T2T build uses a masked Y, and X
alignment in males is a different problem than in females — so this may be a ploidy-handling
artifact rather than a result. Stratified XX and XY maps already exist at
`runs/sex_stratified_maps_v4_nanfix_20260810` and
`runs/sex_stratified_nominal_v4_nanfix_20260810`.

If it is an artifact, a headline currently on the page has to come down.

Then **decompose Δz into Δβ and Δse** together with **map representation-dependence by
allele-frequency concordance** — both cheap, both mechanistic, and together they establish what
kind of difference each axis produces. Then **characterise genes testable in one reference
only**, then **add the variant-caller axis**.

## Complete

| Workstream | What it established | Run root |
|---|---|---|
| Concordance baseline | Concordance. Aligner swap r = 0.997, 97.4% identical calls at \|z\|>4, median \|Δz\| = 0.002. Reference swap r = 0.94, 73–75% identical calls, median \|Δz\| = 0.17. 91–95% of genes unmoved. | `gene_discordance_disease_20260816` |
| Discordant regions | Window discordance at 5% FDR: cut z 2.70–2.81, 1,350–1,858 windows per contrast (5–7% of testable genome). | `window_discordance_exceedance_20260816` |
| Discordant regions | Direction: strength balanced (580 up / 595 down); signed effect skewed ~2:1 down. Up and down windows are indistinguishable in sequence context (p = 0.54–0.94); both differ from the rest of the genome. | `directional_window_discordance_20260816` |
| Gene classes | Method-sensitive genes enriched for ClinGen recurrent-CNV regions (OR 2.9–9.1), human accelerated regions (OR 1.7–1.9), segmental duplications (OR 1.7–2.5) — all four contrasts, 12/12 positive. MHC OR 2.5–3.8. | `discordant_gene_classes_20260816` |
| Gene classes | Aligner effect is only ~50% portable across references (gene Jaccard 0.495); reference effect is 92% portable across aligners (0.919). | `gene_discordance_disease_20260816` |
| Calibration | **Negative, and worth keeping:** genome-wide family-wise rotation test returns 0 of 106,342 window tests. A circular shift relocates the signal rather than removing it, so the null max ≈ observed max (ratio 1.02). A shift null tests *position*, not *magnitude*, and cannot certify these windows. | `genomewide_window_discordance_fwer_20260816` |

## Known risks

- **chrX result may not survive *Separate the chrX effect by sex*.** See next action.
- **Enrichments not yet controlled for expression level.** Duplicated-region genes are often
  lowly expressed; some of the *Characterise discordant regions and gene classes* enrichment could be a power artifact. Open sub-task on *Characterise discordant regions and gene classes*,
  resolved by *Describe background-noise properties per arm*.
- **All FDR estimates are lower bounds.** Normal-null tail assumption; LD correlates
  neighbouring windows and genes.
- **No ground truth exists in this dataset.** No result can be phrased as accuracy. See
  `ANALYSIS_PLAN.md`.
- **SuSiE reuse is unverified.** Candidate outputs are e36-era, not the current v4/k35 arms.
  Confirm arm correspondence before building *Count credible sets per gene* on them.
- **`collect_report_data.py` lives outside version control**, in the analysis workspace. It is
  the single point of failure for reproducing the page. Moving it into this repo is a small,
  worthwhile task nobody has done.

## Page status

Branch `hotspot-section` in this repo, **not merged and not pushed**. Live site
(`josephlalli.github.io/brainvar-t2t-eqtl-report`) is at commit `4890ce5` and reflects none of
the following:

- NaN-bug narrative removed; page retitled *Mapping Brain eQTLs on Two Reference Genomes*
- "Replication" framing replaced with the factorial reference-by-aligner reading
- Genome-wide FDR section, directional section, and their figures added

Publishing is a fast-forward of `main` to this branch and a push. Do not merge or push without
asking.

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
