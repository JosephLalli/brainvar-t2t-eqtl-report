# Brief: can a GRCh38-based GWAS see into hypervariable regions?

**Purpose of this document.** It is a prompt. Give it to a language model and ask for an
execution plan; the plan it returns should be detailed enough for a *second* model, with no
prior context, to carry out. Everything below is background the planner needs. The explicit
request is at the end.

---

## 1. The scientific claim being tested

A four-arm cis-eQTL study (2 references × 2 aligners) on the BrainVar developmental brain
cohort has established that reference and aligner choice leaves most of the genome alone but
moves a specific minority, and that the minority is not random. Method-sensitive genes are
enriched 2.9–9.1× for recurrent genomic-disorder regions, 1.7–1.9× for human accelerated
regions, and 1.7–2.5× for segmental duplications, consistently across all four contrasts.

The governing hypothesis is:

> Regions that are hard to align are hard *because of population-level variability*. That same
> variability is what makes them disease-relevant. So fields with the longest history of
> difficulty — psychiatry, neurodevelopment, immunology — are disproportionately exposed to a
> variable most of the field treats as inert. For most of genetics reference choice does not
> matter; for a small, interesting, historically difficult subset it does.

## 2. The result that motivates this experiment

Bayesian colocalisation (`coloc.abf`) was run genome-wide across all four arms against 40 GWAS
(12 psychiatric, 5 neurodevelopmental, 2 neurodegenerative, 10 immune, 11 comparison
phenotypes), giving 1,031,475 gene-trait tests. The prediction was that psychiatric and
neurodevelopmental colocalisations would be the most sensitive to a reference swap.

**They were the least sensitive.** Pooling psychiatric with neurodevelopmental, 26.7% of
colocalisation calls changed under a reference swap against 34.2% for everything else — odds
ratio 0.70, p = 0.014. Two obvious rescues were tested and both failed: the pattern is not
explained by the data source (within the psychiatric group, GWAS-Catalog-sourced and
consortium-sourced studies give 27.8% and 28.6%, p = 1), and it is not explained by brain
colocalisations being more confident (neurodevelopmental traits have the *lowest* median
posterior among their calls, 0.873 against 0.910 for immune, so they sit closer to the
calling threshold and should flip more readily).

**But the test may be structurally incapable of detecting the effect**, for two reasons.

*Ascertainment.* A gene-trait pair is only testable where the GWAS already has an association
inside the gene's cis window. The loci examined are therefore the ones GWAS has already
succeeded on. If representation choice bites hardest where GWAS has *not* yet produced clean
signals, those loci are excluded by construction.

*Coverage, measured.* GWAS coverage of the common variant frame is depleted inside the
project's FDR-flagged discordance windows, for every trait, and **most for brain traits**:

| Trait group | Coverage inside flagged windows ÷ outside |
|---|---|
| immune | 0.870 |
| comparison | 0.862 |
| neurodegenerative | 0.798 |
| psychiatric | 0.777 |
| neurodevelopmental | 0.765 |

Critically, the standard-error ratio at matched allele frequency is ≈1.000. **Conditional on a
variant being measured there is no power loss; the studies simply do not have variants there.**
That is the signature of an array-plus-imputation coverage artefact rather than of quiet
biology. Two confounds remain unresolved: flagged windows are segmental-duplication enriched,
so some depletion is expected for any technology; and the psychiatric studies are mostly older
consortium releases, so the trait difference could be an imputation-panel vintage effect.

*The structural ceiling.* Every GWAS in the panel is GRCh38-based. Any variant only T2T-CHM13
can represent is absent from the trait side by construction. The colocalisation analysis is
therefore blind to precisely the mechanism the hypothesis is about, which makes a null result
close to unfalsifiable.

## 3. The experiment to plan

**Determine whether GRCh38-based GWAS are blind in the regions where the hypothesis predicts
neuropsychiatric disease variants concentrate, or whether those regions are adequately
interrogated and the hypothesis is wrong.**

The definitive version is: re-impute a cohort against a T2T-CHM13 haplotype reference panel,
re-run the association analysis, and compare discovery against the GRCh38-based equivalent in
flagged versus unflagged regions.

A relevant prior result exists and should be built on:

> Lalli, Bortvin, McCoy & Werling (2025). *A T2T-CHM13 recombination map and globally diverse
> haplotype reference panel improves phasing and imputation.* bioRxiv 10.1101/2025.02.24.639687.
> Reports 38% fewer assembly-discordant genotypes and 16% fewer switch errors on T2T-CHM13,
> **with the largest gains on chromosome X and in regions flanking disease-causing CNVs** —
> i.e. in the same class of regions this project flags.

**The planner should treat "re-run a GWAS" as one of several designs and recommend among
them**, because the full version requires individual-level genotype and phenotype data for a
psychiatric cohort, which is controlled-access and may be a hard blocker. In particular
consider:

- **A. Full re-imputation and re-association.** Requires individual-level data (PGC via
  dbGaP/controlled access, or an institutional cohort). Highest evidential value, highest
  access risk.
- **B. Imputation-accuracy experiment without phenotypes.** Take a WGS truth set, mask to a
  common array manifest, impute against a GRCh38 panel and against the T2T panel, and compare
  recovery — call rate, r² against truth, and MAF-stratified accuracy — inside flagged versus
  matched unflagged windows. **This tests the blindness claim directly and needs no phenotype
  data at all.** Likely the highest value per unit of access risk; the planner should say
  whether it agrees.
- **C. Panel-site accounting.** Compare which sites exist in HRC / TOPMed / 1000G-GRCh38
  panels against the T2T panel, inside flagged windows. Cheapest; establishes an upper bound
  on what any GRCh38 imputation could recover.
- **D. Positive control from existing data.** Identify traits with both an array-imputed and a
  WGS-based GWAS and compare signal density in flagged windows.

## 4. Resources that exist

The planner should assume these are available but **must specify how to verify each** rather
than assuming paths are correct.

- BrainVar: 225 donors, developmental brain, short-read WGS plus RNA-seq, already processed
  through four alignment/reference arms with joint-called genotypes retaining depth and
  allele-balance fields.
- A four-arm variant identity table mapping every variant between native GRCh38/T2T
  coordinates and a LiftoverIndel-normalised common T2T frame, including a per-variant
  orientation sign (REF/ALT swap for ~30% of variants across that boundary).
- FDR-flagged discordance windows genome-wide (100-kb, 5–7% of the testable genome) with
  mappability, segmental-duplication, repeat and HAR fractions attached.
- 40 GWAS summary-statistic sets already harmonised into the common frame.
- `LiftoverIndel` (indel-aware liftover with haplotype realignment), GRCh38↔CHM13v2 chains,
  a reference-differences VCF, and both assemblies' FASTA.
- Compute: a shared 256-core, 1 TB machine. **Constraints: keep under 96 cores and 350 GB
  resident, and leave at least 3 TB of disk free.** No sudo. GPUs available (2× NVIDIA L4).

## 5. What to produce

Return a plan that a second model could execute without further context. It must contain:

1. **A recommendation** among designs A–D (or a better one), with the reasoning made explicit
   and the decisive weakness of each named.
2. **Source files and how to obtain them** — reference panels, array manifests, truth sets,
   cohort data — with URLs or accession identifiers, approximate sizes, and for each one
   whether it is open, registered-access or controlled-access. Flag every step that requires a
   human to sign something, since that is the rate-limiting path.
3. **Methods, concretely.** Named tools and versions for phasing and imputation (e.g. Eagle,
   Beagle, SHAPEIT, minimac, GLIMPSE), the exact accuracy metrics, and the statistical
   comparison — including how flagged and unflagged windows should be matched, given that
   flagged windows differ systematically in mappability, segmental-duplication content and
   variant density, so an unmatched comparison would be confounded.
4. **The discriminating predictions.** For each analysis, state what result would support
   "GRCh38 GWAS are blind here" and what would support "adequately powered, hypothesis
   refuted". A design that cannot produce both outcomes should be discarded.
5. **Correctness gates.** What must be checked before any result is believed — at minimum,
   how to confirm the two imputation runs are comparable (same samples, same masked sites,
   same MAF spectrum) so that a difference cannot be an artefact of the comparison itself.
6. **Cost estimate** in CPU-hours, wall-clock and disk, against the constraints in §4.
7. **A staged plan with abandonment criteria** — what to run first, and what result would
   justify stopping rather than continuing.

## 6. Standing constraints on how the work is framed

- The project describes differences; it does not crown a winner. No analysis should be
  designed to show that an arm is better, and no result should be written that way.
- There is **no ground truth in BrainVar** — short-read sequencing, no per-donor assemblies —
  so nothing can establish which arm is *correct*, only which differ and where. Any design
  that needs ground truth must bring its own (which is why an external WGS truth set matters
  in design B).
- Negative results are kept and reported, not discarded. The colocalisation result in §2 is a
  negative for the project's own hypothesis and is published as such.
- "Genotype error versus real biology" is a false dichotomy here: identifying where measurement
  is representation-dependent is *how* real biology becomes visible.
