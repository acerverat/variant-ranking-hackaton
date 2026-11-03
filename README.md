# Variant Ranking — Hackathon

Prioritization of variants from a whole-genome VCF for clinical diagnosis,
combining functional annotation (ANNOVAR), phenotype-gene similarity (HPO),
and an auditable composite score (clinical evidence + impact + frequency +
inheritance + quality control).

## Structure

```
annotation/   Nextflow pipeline: VCF → NORMALIZE → ANNOVAR → SELECT_CANDIDATES
ranking/      HPO-gene similarity + composite score → final ranked table
```

See each folder's `README.md` for step-by-step detail, parameters, and
design decisions (why each thing is or isn't penalized).

**This repository is code only.** Outputs from running the pipeline
(candidate variants, ranked tables, Nextflow reports) are not versioned here
(`results/` is in `.gitignore`) — they're generated locally when you run it.

## Pipeline, at a glance

```
VCF (WGS, GRCh38, 5M variants)
  │  NORMALIZE: PASS, splits multiallelics, aligns indels, basic QC
  │  ANNOVAR: 9 databases (refGene, ensGene, avsnp151, ClinVar, gnomAD v4.1
  │           genome+exome, dbNSFP47a, dbscSNV, InterVar/ACMG)
  │  SELECT_CANDIDATES: exonic/splicing, non-synonymous, rare (≤1% gnomAD)
  ▼
640 candidate variants (see results/candidates/ once you run it)
  │  hpo_gene_similarity.py: Resnik+BMA similarity against the HPO phenotype
  │  composite_score.py: 0-100 score across 5 axes + quality flags
  ▼
final ranked table (see results/*.top25.tsv once you run it)
```

## Result (case WGS_EX2312012)

| # | Score | Gene | Variant | Evidence |
|---|---|---|---|---|
| 1 | 91.5 | `BUB1B` | stopgain p.L737X | ClinVar P/LP, InterVar Pathogenic, exact HPO match (rhabdomyosarcoma) |
| 2 | 78.2 | `LZTR1` | stopgain p.Y748X | ClinVar Pathogenic — secondary finding |
| 3 | 72.4 | `BUB1B` | missense p.N1002K | AlphaMissense 0.92, same gene as #1 |
| 4 | 65.3 | `GNRHR` | missense p.R139H | ClinVar Pathogenic, weak phenotype match |

**Leading hypothesis:** mosaic variegated aneuploidy syndrome 1 (MVA1,
`BUB1B`), recessive — consistent with rhabdomyosarcoma, growth restriction,
and the parents' recurrent miscarriages described in the phenotype. Phase
between the two `BUB1B` variants is unconfirmed (no parental sample).

## Known limitations

- Score weights are hand-set, not learned — a simulated-case benchmark
  (ClinVar + 1000 Genomes) is still pending for independent validation.
- Coverage limited to exonic/splicing (±10 bp); no CNVs/SVs, deep intronic
  variants, or mitochondrial.
- No trio: no recessive variant has confirmed phase.

Full detail in `ranking/README.md`.

## Credits

The annotation pipeline is adapted from
[Pipelines_Inmegen/Annotation](https://github.com/INMEGEN/Pipelines_Inmegen/tree/Principal/Annotation)
(INMEGEN, CC BY-NC-SA 4.0).
