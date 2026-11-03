# Phenotype-gene similarity (HPO)

First component of the ranking model: for each candidate gene, how well its
known phenotype (HPO) matches the patient's.

## Method

Resnik + best-match average (BMA) — the same principle used by Exomiser/Phenolyzer:

1. **IC(t) = -log(frequency(t))**, where frequency(t) = fraction of annotated
   genes (in `genes_to_phenotype.txt`) linked to `t` or a descendant of `t`.
   Rare, specific terms (e.g. "Rhabdomyosarcoma") have high IC; generic terms
   ("Abnormality of X") have IC ~0.
2. **sim(t1, t2)** = IC of the most informative common ancestor (MICA) between
   `t1` and `t2`.
3. For each patient term, the best match among the gene's terms is taken and
   averaged → **`HPO_SIM_P2G`** (the clinically relevant direction: how well
   the gene explains what the patient has). The symmetric version,
   **`HPO_SIM_BMA`**, is also reported.

Implemented in `bin/hpo_gene_similarity.py`, Python standard library only.

## Data (`refs/`, downloaded 2026-09-29)

- `hp.obo` — HPO ontology, release `2026-09-01` ([hpo.jax.org](https://hpo.jax.org))
- `genes_to_phenotype.txt` — gene↔HPO annotations (aggregates every disease linked to each gene)

## Usage

```bash
python3 bin/hpo_gene_similarity.py \
  --candidates ../annotation/results/candidates/WGS_EX2312012.candidates.tsv \
  --hpo-terms "HP:0002859,HP:0000121,HP:0004322,HP:0001508,HP:0003202,HP:0001622,HP:0001518,HP:0200067" \
  --out ranked_by_hpo.tsv
```

The HPO terms come from `Challenge_Clinical_Phenotype_1.docx` (rhabdomyosarcoma,
nephrocalcinosis, short stature, failure to thrive, muscle atrophy, prematurity,
SGA, the parents' recurrent miscarriages — the last one included because the
document asks for it to be treated as part of the phenotype, not just history).

Output: `candidates.tsv` + 4 columns (`HPO_SIM_P2G`, `HPO_SIM_BMA`,
`N_GENE_HPO_TERMS`, `HPO_BEST_MATCH`), sorted by `HPO_SIM_P2G` descending.
Genes with no HPO annotation (not linked to any known disease) are left as `NA`.

## Validation result (WGS_EX2312012)

`BUB1B` ranks **#1** out of 640 candidates, via an exact match on
`HP:0002859` (Rhabdomyosarcoma, IC=5.53 — a term rarely shared across genes).
This is consistent with mosaic variegated aneuploidy syndrome 1 (MVA1,
OMIM #257300), already flagged by ClinVar/InterVar for one of its two
candidate variants (see conversation). Confirms the similarity design is
capturing the expected signal.

---

# Composite score (`bin/composite_score.py`)

Combines five evidence axes into a 0-100 score. Each axis is reported as its
own column plus a readable `REASON` string, so the ranking stays
**auditable**: an analyst can see why a variant went up or down.

```
SCORE = 30*S_PHENO + 25*S_EVIDENCE + 25*S_IMPACT + 10*S_FREQ + 10*S_INHERIT − penalties
```

| Axis | Weight | Where it comes from |
|---|---|---|
| `S_PHENO` | 30 | `HPO_SIM_P2G` / 2.0 (fixed scale, not min-max, so runs stay comparable) |
| `S_EVIDENCE` | 25 | ClinVar + `InterVar_automated` + count of pathogenic ACMG criteria |
| `S_IMPACT` | 25 | **by consequence class**: truncating 0.95; splicing via dbscSNV; missense = average of REVEL/AlphaMissense/MetaRNN/CADD |
| `S_FREQ` | 10 | rarity: `(−log10(AF) − 2) / 3`, clamped; absent from gnomAD = 1.0 |
| `S_INHERIT` | 10 | hom/hemizygous = 1.0; het with ≥2 candidates in the gene = 0.8; isolated het = 0.4 (+bonus if the gene is constrained) |

### Design decisions that matter

- **In-silico scores only apply to missense.** REVEL/AlphaMissense don't exist
  for stopgain/frameshift. Blindly averaging them would sink truncating
  variants, which are the highest-impact ones. That's why `S_IMPACT` is
  computed per consequence class, and missing predictors on an LoF don't
  penalize it.
- **LOEUF/pLI never penalize a recessive candidate.** `BUB1B` has pLI ≈ 0
  precisely because heterozygous loss is tolerated: using pLI as a general
  filter would rule out this kind of diagnosis. Constraint only gives a
  **bonus** to the dominant hypothesis.
- **No parents means no phase and no de novo.** Two candidates in the same
  gene are flagged as "possible compound het (PHASE UNCONFIRMED)"; an
  isolated het gets partial credit.
- **A gene with no HPO annotation isn't zeroed out**, it gets a floor
  (`--pheno-floor 0.15`): it could be a novel finding, even if it's less
  likely to explain a recognizable presentation.
- **Quality flags** (penalize, don't eliminate — a real finding stays visible):
  - `FLAG_NOISY_GENE` (127/640): hyperpolymorphic genes or ones prone to
    artifactual LoF calls with short reads — HLA cluster, mucins,
    `FLG`/`HRNR`, olfactory receptors, `ZNF`, giant genes. Without this,
    16 of the top 60 were HLA.
  - `FLAG_AB` (79/640): allele balance inconsistent with the called genotype
    (het with AB < 0.25 or > 0.75; AB ≈ 0 or 1 on a het is contradictory).

### Exception: chromosomes with excess genome-wide skewed hets (`chrom_ab_scan.py`)

`FLAG_AB` assumes a skewed AB is a *local* artifact (mapping, PCR). But this
phenotype is a chromosomal instability syndrome — a skew *shared across
thousands of variants on the same chromosome* would be the expected signature
of a large mosaic aneuploidy/CNV (the same principle as the B-allele
frequency of a cytogenomic array), not a pinpoint artifact.

`bin/chrom_ab_scan.py` scans the **whole** normalized VCF (not the 640
candidates — a whole chromosome only has 0-2 of those, no statistical power)
and computes, per chromosome, what fraction of high-depth hets fall outside
0.35-0.65. A chromosome is flagged if that rate exceeds the other autosomes'
by a robust z-score (median+MAD) ≥3.

**Result in this case:** chr20, 21 and 22 stand out clearly (z=6.2-6.6) from
the rest (z between −0.9 and 2.2). When a candidate variant falls on one of
those chromosomes *and* also triggers `FLAG_AB`, the AB penalty is reduced to
0.3× (`--chrom-skew-factor`) instead of being dropped entirely — because a
real aneuploidy can't be reliably distinguished from a technical mapping bias
that **also** happens to be specific to those chromosomes (GC-rich,
acrocentric — see limitation below). It only affected 2 `CYP2D7` variants
(penalty 10→3, +7 points each); the top 20, including `BUB1B`, didn't change.

```bash
python3 bin/chrom_ab_scan.py \
  --vcf ../annotation/results/normalized/WGS_EX2312012.norm.vcf.gz \
  -o results/WGS_EX2312012.chrom_ab_summary.tsv
```

**Limitation that remains unresolved:** this signal alone does NOT
distinguish a real mosaic aneuploidy from a chromosome-specific technical
mapping bias — chr20/21/22 already showed, from the original genome-wide
check, more depth than expected (52-58x vs ~44x elsewhere), a pattern
typically associated with GC-content/mappability on acrocentric chromosomes
with short reads, not necessarily biological. Confirming it for real
requires: comparing against the same pipeline run on control samples (is it
universal to the pipeline, or specific to this patient?), a dedicated CNV
caller (binned depth + BAF, not just the mean), or karyotype/FISH — the real
clinical standard for MVA1, since the mosaicism is stochastic cell-by-cell
and may not leave a clean signature in bulk WGS.

## Usage

```bash
python3 bin/composite_score.py \
  --candidates results/WGS_EX2312012.ranked_by_hpo.tsv \
  --constraint refs/gene_constraint.tsv \
  --out results/WGS_EX2312012.ranked_composite.tsv --top 25
```

Weights and thresholds are parameters (`--w-pheno`, `--ab-lo`, `--no-quality-flags`, …).

`refs/gene_constraint.tsv` — LOEUF/pLI per gene, derived from
`gnomad.v4.1.constraint_metrics.tsv` (one row per gene, preferring MANE Select; 18,203 genes).

## Result (WGS_EX2312012)

Deliverables: `results/WGS_EX2312012.ranked_composite.tsv` (640 rows) and
`results/WGS_EX2312012.top25.tsv` (top-N with key columns).

| # | Score | Gene | Variant | Evidence |
|---|---|---|---|---|
| 1 | **91.5** | `BUB1B` | `c.2210T>G p.L737X` stopgain | ClinVar P/LP (MVA1), InterVar Pathogenic, PVS1+PM2+PP3, exact HPO match on rhabdomyosarcoma |
| 2 | 78.2 | `LZTR1` | `c.2244C>G p.Y748X` stopgain | ClinVar Pathogenic, virtually absent from gnomAD — **secondary finding to consider** (schwannomatosis/Noonan, dominant) |
| 3 | 72.4 | `BUB1B` | `c.3006T>G p.N1002K` missense | AlphaMissense 0.92, CADD 24.5, AF 6.8e-07 |
| 4 | 65.3 | `GNRHR` | `c.416G>A p.R139H` missense | ClinVar Pathogenic, but weak phenotype match (0.70) |

The two `BUB1B` variants land at #1 and #3, with a clear gap from the rest
(cohort median: 22.0). Both have clean allele balance (0.54 / 0.46) and GQ 99.

## Limitations (important)

- **Weights are hand-set, not learned.** With a single case and no labeled
  data there's nothing to fit a model against; they're a reasonable *prior*,
  not a calibration.
- **Validation with `BUB1B` is not independent.** `BUB1B` was first found
  through a manual panel review, and the score's design (in particular the
  decision that LOEUF/pLI shouldn't penalize recessive genes) was made
  knowing the candidate gene was recessive. Real validation needs a
  simulated-case benchmark: insert known pathogenic ClinVar variants into
  background genomes (e.g. 1000 Genomes) together with the disease's HPO
  terms, and measure where the score ranks them (top-1, top-10, median rank).
  That's the real next step, and it's also what would allow training a
  *learning-to-rank* model (LightGBM) instead of fixed weights.
- **Coverage**: exonic and ±10 bp splicing only. Deep intronic and regulatory
  variants, CNVs/SVs, repeat expansions, and the mitochondrial genome are all
  left out (excluded in `NORMALIZE`; needs its own pipeline because of
  heteroplasmy).
- **Phase is unresolved**: the two `BUB1B` variants are ~10.9 kb apart, too
  far to link with short reads. Confirming *trans* requires a parental
  sample (targeted Sanger would do) or long-read sequencing.
- `S_PHENO` depends on the gene already being linked to disease in HPO; novel
  genes are at a disadvantage by construction.
