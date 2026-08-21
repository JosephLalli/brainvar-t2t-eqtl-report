# Review findings, 2026-08-19 — read before editing the page

Three adversarial reviews of the results page (statistics; claims-versus-evidence; methods and
reproducibility) were run against the page, the generator, the docs and the run tree. This
records what was **verified against the run tree**, what is **reported but unverified**, and
what survives. Roughly 45 findings were returned; the ones that change conclusions are here.

**The page is currently live at the state described below and carries at least two claims that
are demonstrably wrong.** Fix those before adding anything new.

---

## Verified — these were checked directly against the run tree

### 1. Expression is quantified per *reference*, not per *arm* — the axes are not commensurable

`runs/native_expression_all_contigs_20260808/` contains **two** matrices, `grch38/` and `t2t/`,
for four arms. An aligner swap therefore holds the phenotype exactly fixed; a reference swap
changes genotypes **and** re-quantifies RNA. Nothing on the page says so.

Quantified in `runs/expression_confound_decomposition_20260819` (14,080 genes with both
contrasts and a cross-reference expression correlation):

| Stratum | genes | ref median \|Δz\| | aligner median \|Δz\| | ratio | ref flag | aligner flag |
|---|---|---|---|---|---|---|
| all genes | 14,080 | 0.2152 | 0.0225 | **9.56** | 8.62% | 5.19% |
| r < 0.99 | 1,489 | 0.4252 | 0.0253 | 16.84 | 50.03% | 16.86% |
| r > 0.99 | 12,591 | 0.2061 | 0.0223 | **9.23** | 3.72% | 3.81% |
| r > 0.999 | 3,805 | 0.1768 | 0.0200 | **8.84** | 0.39% | 2.16% |
| r > 0.9999 | 75 | 0.1618 | 0.0185 | **8.76** | 0.00% | 1.33% |

**The magnitude ratio survives** (9.56 → 8.76, ~8% attenuation): the reference really does
perturb effect estimates about nine times as much as the aligner, and that is not an RNA
artefact. **The flagged set does not.** The reference flag rate collapses 22-fold as expression
stabilises while the aligner's falls 2.4-fold, and at r > 0.999 the reference flags *fewer*
genes than the aligner. The two statistics disagree because the flag is a count-matched robust
z — an outlier score computed within a contrast against its own baseline — so the reference sits
uniformly higher (real) while the genes standing out relative to that baseline are
disproportionately the re-quantified ones (confounded).

**What this threatens:** the identity of the reference-axis method-sensitive set, and therefore
everything computed on it — the ClinGen genomic-disorder, HAR and segdup enrichments *on the
reference axis*, and the reference-axis hotspot map. Multi-mapping RNA quantification fails
hardest in duplicated sequence, which is exactly where those enrichments live, so the
alternative explanation is live. 61.9% of reference-flagged genes have expression r < 0.99
against 5.9% of concordant genes.

**What it does not threaten:** anything on the aligner axis, where expression is held exactly
fixed — the genotype-fidelity result, the chrX XX-specific finding, the aligner-axis
colocalisation. Nor the allele-frequency or variant-identity work, which involve no expression.

**The clean test has now been run** — `runs/crossed_reference_association_20260821`, recorded in ANALYSIS_PLAN.md. It confirms the confound and quantifies it: of a 21.7% eGene turnover under a reference swap, expression alone accounts for 20.1% and genotypes alone for 6.6%. The control cell reproduces the published arm at 0.9% turnover and r = 0.9999, validating the machinery.

**The threatened enrichments have now been retested directly, and they hold.**
`runs/genotype_term_gene_classes_20260821` recomputes the gene-class enrichments on the
genotype term alone — one expression matrix, one covariate set, one annotation, only the
genotypes swapped — using the same recipe, the same class tracks and the same Fisher test.
All three survive and two strengthen:

| Class | Reference axis, as published | Genotype term alone |
|---|---|---|
| recurrent genomic-disorder region | 2.91x | **4.02x** (p = 4e-24) |
| segmental duplication | 1.65x | **2.40x** (p = 3e-24) |
| human accelerated region | 1.84x | 1.54x (p = 0.0072) |

The alternative explanation this finding raised — that multi-mapping RNA quantification fails
hardest in duplicated sequence, which is where the enrichments live — is therefore **not** what
carries them. The genotype term flags a smaller set (579 genes against
1,256) that is markedly more concentrated in those regions
(12.8% of genomic-disorder genes against 3.5% elsewhere).
The reference effect is mostly RNA **by volume** and genotype **by mechanism**.

Three limits stay attached to that. The two flagged sets differ in size and cut, so the odds
ratios are comparable in direction rather than to the second decimal. The genotype term sees
only variants both references can represent, so it excludes the access effect and is a **lower
bound**. And neither side is matched on mappability, gene density or expression level, so what
is shown is that disease genes are present in this sequence, not why. That last point is now
carried on the page and listed under "not settled".

The re-signing this required was validated independently:
`runs/crossed_reference_nominal_20260821/FLIP_VALIDATION.json`. Re-signed and untouched
variants agree on allele frequency to the same median absolute difference (0.0000); had the
sign been inverted the re-signed set would sit 0.39 away.

### 2. The colocalisation "independent instrument" claim is false

The page says the stratified colocalisation reproduces the regional claim "on an instrument
that shares no statistics with the one that defined the windows". Both sides read the **same
nominal scan**: `runs/gwas_coloc_bayesian_20260817/MANIFEST.json` gives `eqtl_source:
runs/four_arm_nominal_v4_k35_nanfix_20260811`, and
`runs/noninteraction_int_hotspots_100kb_20260812/MANIFEST.json` lists `nominal_run_record`
pointing into that same run. A window is flagged *because* the arms' β/SE vectors differ most
there, so "calls change 3.4× more often inside flagged windows" is close to algebraic. The GWAS
side is independent but is common to both arms and cancels from the contrast.
**Fix:** present it as an internal-consistency check, or re-test using a window definition
disjoint from the eQTL statistics (sequence context, mappability, or held-out chromosomes).

### 3. The median |ΔPP4| reassurance statistic is conditioned on the wrong denominator

The page reports median |ΔPP4| = 0.0021 over all 250,410 tested gene-trait pairs, most of which
carry no signal in either arm. Restricted to pairs where a call exists it is **0.0204** (~10×),
and among PP4 ≥ 0.5 pairs, 24% move by more than 0.2. The reported Pearson r = 0.915 is
inflated the same way.

### 4. "chrX dosages are encoded identically" is not what the audit shows

`runs/chrx_discordance_by_sex_20260816/dosage_audit.json`: XY het-like share is 0.0252 /
0.0176 / 0.0172 / 0.0194 (**1.47× spread**) and `chrX_rows_seen` is 207,004 / 307,491 /
304,294 / 261,596 (**1.49× spread**). The two arms of the GRCh38 aligner contrast differ by 48%
in how much chrX they test. That sentence is the control guarding the OR 7.27 result.

### 5. A manifest contradicts its own directory

The same run's `MANIFEST.json` records `"dosage_audit": {"skipped": "no sex labels resolved"}`
while the completed audit (92 XX, 133 XY) sits beside it.

### 6. A surviving enrichment result goes unreported

The page says no term survives correction "in three of the four contrasts" and then discusses
sub-threshold immune terms at length. `runs/gene_discordance_disease_20260816` shows the fourth
contrast (`graph − linear · T2T`) returns **HP:0001419 X-linked recessive inheritance,
p = 2.8e-6** and HP:0001417 X-linked inheritance. The above-threshold result is omitted.

### 7. The page contradicts itself on the yardstick

"The yardstick is missing" and "a reference swap is the smallest of the three method changes"
both appear. The first is stale text that survived the yardstick workstream. The lead-variant
LD caveat ("has not been computed") is stale in the same way.

---

## Reported by reviewers, not independently verified

Worth checking before relying on any of them, but each is specific enough to act on:

- **The window FDR is uncalibrated.** The only simulated null (rotation) returned zero; its
  replacement is a normal upper tail applied to a demonstrably right-skewed statistic. Reported
  FDP runs 5% → 27% across plausible null scales.
- **chrX null cells use a one-sided Fisher test** (`alternative="greater"`), so a 4–6×
  *depletion* under an aligner swap is reported as "nothing" at p = 1.
- **Mantel–Haenszel pooling across heterogeneous strata** (odds ratios 0.55, 0.49, 1.38, 1.01,
  0.78) with no homogeneity test, and severe stratum imbalance.
- **A single coloc prior variance** W = 0.2² applied to all 40 studies including ~11
  quantitative traits used in native units; brain-trait effect sizes are additionally
  *reconstructed* by two different formulas, and source is near-collinear with trait area.
- **"was a study-power artefact"** overstates a widened interval whose CI still contains the
  original estimate; the page says the correct thing two paragraphs later.
- **The ancestry table shows one of two aligner contrasts** and the other is null (p = 0.19).
- **The stratified colocalisation was run for only 2 of 4 contrasts.**
- **Methods omit the primary scan's cis window, MAF ≥ 0.05, permutation count and seed rule**;
  no software versions anywhere; aligners never named.
- **49 of 77 run roots have no manifest**; 22 of 27 existing manifests carry no SHA-256 despite
  `CONVENTIONS.md` requiring it.
- **The gene-class enrichment confound is stated in one manifest and not the other.**
  `neuro_gene_sets_20260817` records that flagged windows differ in mappability, segdup content
  and gene density and "this test does not match on those" — but
  `discordant_gene_classes_20260816`, which backs the page's central claim, carries no such
  caveat. Same confound, applied only to the test that produced an unwelcome answer.
  *Partly addressed:* `genotype_term_gene_classes_20260821` carries the caveat explicitly,
  and the page now states it in the classes section and under "not settled". The original
  `discordant_gene_classes_20260816` manifest is unchanged, and no matched test has been
  run on either axis.
  *Now closed:* `matched_gene_class_enrichment_20260821` runs the matched test on both
  axes, matching on mappability, gene density, expression level, cis-variant count and
  gene length. The power covariates explain none of the enrichment; mappability explains
  most of it, and all of the human-accelerated enrichment on the genotype term.
- **eGene counts disagree** between the PC-sweep and maps tables at the same k; the "12,189"
  headline is a sum over four overlapping arms, not a distinct-gene count.
- **Several quoted ranges do not match their tables** (2.2× vs 2.0–3.0; 56–69% vs 56.0–70.8%).

---

## What reviewers confirmed is sound

Orientation handling across the common frame is correct and consistent throughout (genotypes
recoded on flip, AF complemented, unordered allele keys, both effect vectors re-signed). SuSiE
and coloc parameters on the page match the scripts exactly. The inverse-variance-weighted
selection rule for the between-sex contrast is a textbook handling of selection-induced bias.
Boundary propagation from manifests to page is better than typical.

---

## Suggested order of work

1. **The crossed association** (§1). It changes conclusions rather than wording, so do it before
   rewriting text that may need to change again.
2. **The verified factual errors** (§2–§7). Quick, and two of them are live on a public page.
3. **The statistical items**, in the order listed above.
4. Provenance and Methods completeness last — tedious, but nothing else depends on it.
