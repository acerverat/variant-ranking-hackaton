# Variant Ranking — Hackatón

Priorización de variantes de un VCF de genoma completo para diagnóstico clínico,
combinando anotación funcional (ANNOVAR), similitud fenotipo-gen (HPO) y un score
compuesto auditable (evidencia clínica + impacto + frecuencia + herencia + control
de calidad).

## Estructura

```
annotation/   Pipeline de Nextflow: VCF → NORMALIZE → ANNOVAR → SELECT_CANDIDATES
ranking/      Similitud HPO-gen + score compuesto → tabla final rankeada
```

Ver el `README.md` de cada carpeta para el detalle de cada paso, los parámetros,
y las decisiones de diseño (por qué se penaliza o no cada cosa).

**Este repositorio es solo código.** Las salidas de correr el pipeline (variantes
candidatas, tablas rankeadas, reportes de Nextflow) no se versionan aquí
(`results/` está en `.gitignore`) — se generan localmente al correrlo.

## Pipeline, de un vistazo

```
VCF (WGS, GRCh38, 5M variantes)
  │  NORMALIZE: PASS, separa multialélicos, alinea indels, QC básico
  │  ANNOVAR: 9 bases (refGene, ensGene, avsnp151, ClinVar, gnomAD v4.1
  │           genoma+exoma, dbNSFP47a, dbscSNV, InterVar/ACMG)
  │  SELECT_CANDIDATES: exónico/splicing, no sinónimo, raro (≤1% gnomAD)
  ▼
640 variantes candidatas (ver results/candidates/ al correrlo)
  │  hpo_gene_similarity.py: similitud Resnik+BMA contra el fenotipo HPO
  │  composite_score.py: score 0-100 con 5 ejes + banderas de calidad
  ▼
tabla final rankeada (ver results/*.top25.tsv al correrlo)
```

## Resultado (caso WGS_EX2312012)

| # | Score | Gen | Variante | Evidencia |
|---|---|---|---|---|
| 1 | 91.5 | `BUB1B` | stopgain p.L737X | ClinVar P/LP, InterVar Pathogenic, HPO exacto (rabdomiosarcoma) |
| 2 | 78.2 | `LZTR1` | stopgain p.Y748X | ClinVar Pathogenic — hallazgo secundario |
| 3 | 72.4 | `BUB1B` | missense p.N1002K | AlphaMissense 0.92, mismo gen que #1 |
| 4 | 65.3 | `GNRHR` | missense p.R139H | ClinVar Pathogenic, fenotipo débil |

**Hipótesis principal:** síndrome de aneuploidía variegada en mosaico tipo 1
(MVA1, `BUB1B`), recesivo — coherente con rabdomiosarcoma, restricción de
crecimiento y abortos recurrentes de los padres descritos en el fenotipo. Fase
entre las dos variantes de `BUB1B` no confirmada (sin muestra parental).

## Limitaciones conocidas

- Pesos del score fijados a mano, no aprendidos — pendiente un set de casos
  simulados (ClinVar + 1000 Genomes) para validación independiente.
- Solo cobertura exónica/splicing (±10 pb); sin CNVs/SVs, intrónicas profundas
  ni mitocondrial.
- Sin trío: ninguna variante recesiva está confirmada en fase.

Detalle completo en `ranking/README.md`.

## Créditos

El pipeline de anotación está adaptado de
[Pipelines_Inmegen/Annotation](https://github.com/INMEGEN/Pipelines_Inmegen/tree/Principal/Annotation)
(INMEGEN, CC BY-NC-SA 4.0).
