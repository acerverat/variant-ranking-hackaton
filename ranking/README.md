# Similitud fenotipo-gen (HPO)

Primer componente del modelo de ranking: para cada gen candidato, qué tan bien
coincide su fenotipo conocido (HPO) con el del paciente.

## Método

Resnik + best-match average (BMA) — el mismo principio que usan Exomiser/Phenolyzer:

1. **IC(t) = -log(frecuencia(t))**, con frecuencia(t) = fracción de genes anotados
   (en `genes_to_phenotype.txt`) a `t` o a un descendiente de `t`. Términos raros y
   específicos (p. ej. "Rhabdomyosarcoma") tienen IC alto; términos genéricos
   ("Abnormality of X") tienen IC ~0.
2. **sim(t1, t2)** = IC del ancestro común más informativo (MICA) entre `t1` y `t2`.
3. Para cada término del paciente se toma la mejor coincidencia entre los términos
   del gen, y se promedia → **`HPO_SIM_P2G`** (dirección clínicamente relevante:
   qué tan bien explica el gen lo que tiene el paciente). También se reporta
   **`HPO_SIM_BMA`**, la versión simétrica.

Implementación en `bin/hpo_gene_similarity.py`, solo biblioteca estándar de Python.

## Datos (`refs/`, descargados 2026-09-29)

- `hp.obo` — ontología HPO, release `2026-09-01` ([hpo.jax.org](https://hpo.jax.org))
- `genes_to_phenotype.txt` — anotaciones gen↔HPO (agrega todas las enfermedades ligadas a cada gen)

## Uso

```bash
python3 bin/hpo_gene_similarity.py \
  --candidates ../annotation/results/candidates/WGS_EX2312012.candidates.tsv \
  --hpo-terms "HP:0002859,HP:0000121,HP:0004322,HP:0001508,HP:0003202,HP:0001622,HP:0001518,HP:0200067" \
  --out ranked_by_hpo.tsv
```

Los términos HPO son los del `Challenge_Clinical_Phenotype_1.docx` (rabdomiosarcoma,
nefrocalcinosis, talla baja, falla de medro, atrofia muscular, prematuridad, PEG,
abortos recurrentes de los padres — este último incluido porque el documento pide
tratarlo como parte del fenotipo, no solo como antecedente).

Salida: `candidates.tsv` + 4 columnas (`HPO_SIM_P2G`, `HPO_SIM_BMA`,
`N_GENE_HPO_TERMS`, `HPO_BEST_MATCH`), ordenado por `HPO_SIM_P2G` descendente.
Genes sin anotación HPO (no ligados a ninguna enfermedad conocida) quedan con `NA`.

## Resultado de validación (WGS_EX2312012)

`BUB1B` queda en el **puesto #1** de 640 candidatas, por coincidencia exacta en
`HP:0002859` (Rabdomiosarcoma, IC=5.53 — término poco frecuente entre genes). Esto
es consistente con el síndrome de aneuploidía variegada en mosaico tipo 1 (MVA1,
OMIM #257300) ya señalado por ClinVar/InterVar para una de sus dos variantes
candidatas (ver conversación). Confirma que el diseño de la similitud está
capturando la señal esperada.

---

# Score compuesto (`bin/composite_score.py`)

Combina cinco ejes de evidencia en un score 0-100. Cada eje se reporta como columna
aparte más una cadena `REASON` legible, para que el ranking sea **auditable**: un
analista puede ver por qué subió o bajó cada variante.

```
SCORE = 30*S_PHENO + 25*S_EVIDENCE + 25*S_IMPACT + 10*S_FREQ + 10*S_INHERIT − penalizaciones
```

| Eje | Peso | De dónde sale |
|---|---|---|
| `S_PHENO` | 30 | `HPO_SIM_P2G` / 2.0 (escala fija, no min-max, para comparar entre corridas) |
| `S_EVIDENCE` | 25 | ClinVar + `InterVar_automated` + nº de criterios ACMG patogénicos |
| `S_IMPACT` | 25 | **por clase de consecuencia**: truncante 0.95; splicing vía dbscSNV; missense = promedio de REVEL/AlphaMissense/MetaRNN/CADD |
| `S_FREQ` | 10 | rareza: `(−log10(AF) − 2) / 3`, acotado; ausente en gnomAD = 1.0 |
| `S_INHERIT` | 10 | hom/hemicigoto = 1.0; het con ≥2 candidatas en el gen = 0.8; het aislado = 0.4 (+bono si el gen está restringido) |

### Decisiones de diseño que importan

- **Los scores in-silico solo aplican a missense.** REVEL/AlphaMissense no existen para
  stopgain/frameshift. Promediarlos a ciegas hundiría a las variantes truncantes, que son
  las de mayor impacto. Por eso `S_IMPACT` se calcula por clase de consecuencia y la
  ausencia de predictores en una LoF no penaliza.
- **LOEUF/pLI nunca penalizan a un candidato recesivo.** `BUB1B` tiene pLI ≈ 0 justamente
  porque la pérdida heterocigota se tolera: usar pLI como filtro general descartaría este
  tipo de diagnóstico. La restricción solo da un **bono** a la hipótesis dominante.
- **Sin padres no hay fase ni de novo.** Dos candidatas en el mismo gen se marcan como
  "posible het compuesto (FASE NO CONFIRMADA)"; un het aislado recibe crédito parcial.
- **Un gen sin anotación HPO no se anula a 0**, recibe un piso (`--pheno-floor 0.15`): podría
  ser un hallazgo nuevo, aunque es menos probable que explique un cuadro reconocible.
- **Banderas de calidad** (penalizan, no eliminan — un hallazgo real sigue visible):
  - `FLAG_NOISY_GENE` (127/640): genes hiperpolimórficos o propensos a LoF artefactual con
    lectura corta — cluster HLA, mucinas, `FLG`/`HRNR`, receptores olfatorios, `ZNF`, genes
    gigantes. Sin esto, 16 de las 60 primeras eran HLA.
  - `FLAG_AB` (79/640): balance alélico incompatible con el genotipo llamado (het con
    AB < 0.25 o > 0.75; AB ≈ 0 o 1 en un het es contradictorio).

### Excepción: cromosomas con exceso de hets desviados genoma-completo (`chrom_ab_scan.py`)

`FLAG_AB` asume que un AB sesgado es un artefacto *local* (mapeo, PCR). Pero este fenotipo
es un síndrome de inestabilidad cromosómica — un sesgo *compartido por miles de variantes en
el mismo cromosoma* sería la firma esperada de una aneuploidía/CNV grande en mosaico (el mismo
principio que la B-allele frequency de un arreglo citogenómico), no un artefacto puntual.

`bin/chrom_ab_scan.py` escanea el VCF normalizado **completo** (no las 640 candidatas — un
cromosoma entero tiene ahí 0-2 variantes, sin poder estadístico) y calcula, por cromosoma, qué
fracción de los hets de alta profundidad cae fuera de 0.35-0.65. Un cromosoma se marca si esa
tasa excede, con un z-score robusto (mediana+MAD) ≥3, la de los demás autosomas.

**Resultado en este caso:** chr20, 21 y 22 se separan con claridad (z=6.2-6.6) del resto
(z entre −0.9 y 2.2). Cuando una variante candidata cae en uno de esos cromosomas *y además*
dispara `FLAG_AB`, la penalización de AB se reduce a 0.3× (`--chrom-skew-factor`) en vez de
anularse — porque no se puede distinguir con certeza una aneuploidía real de un sesgo técnico
de mapeo que **también** es específico de esos cromosomas (ricos en GC, acrocéntricos — ver
limitación abajo). Solo afectó a 2 variantes de `CYP2D7` (penalización 10→3, +7 puntos cada
una); el top 20, incluido `BUB1B`, no cambió.

```bash
python3 bin/chrom_ab_scan.py \
  --vcf ../annotation/results/normalized/WGS_EX2312012.norm.vcf.gz \
  -o results/WGS_EX2312012.chrom_ab_summary.tsv
```

**Limitación que no queda resuelta:** esta señal, por sí sola, NO distingue una aneuploidía
mosaica real de un sesgo técnico de mapeo específico de cromosoma — chr20/21/22 ya mostraban,
desde el chequeo genoma-completo original, más profundidad de lo esperado (52-58x vs ~44x en
el resto), un patrón típicamente asociado a GC-content/mappability en cromosomas acrocéntricos
con lectura corta, no necesariamente biológico. Confirmarlo de verdad requiere: comparar contra
el mismo pipeline en muestras control (¿es universal del pipeline o específico de este
paciente?), un llamador de CNV dedicado (profundidad en bins + BAF, no solo la media), o
cariotipo/FISH — el estándar clínico real para MVA1, porque el mosaicismo es estocástico
célula por célula y no necesariamente deja una firma limpia en WGS en bulto.

## Uso

```bash
python3 bin/composite_score.py \
  --candidates results/WGS_EX2312012.ranked_by_hpo.tsv \
  --constraint refs/gene_constraint.tsv \
  --out results/WGS_EX2312012.ranked_composite.tsv --top 25
```

Pesos y umbrales son parámetros (`--w-pheno`, `--ab-lo`, `--no-quality-flags`, …).

`refs/gene_constraint.tsv` — LOEUF/pLI por gen, derivado de
`gnomad.v4.1.constraint_metrics.tsv` (una fila por gen, prefiriendo MANE Select; 18,203 genes).

## Resultado (WGS_EX2312012)

Entregables: `results/WGS_EX2312012.ranked_composite.tsv` (640 filas) y
`results/WGS_EX2312012.top25.tsv` (top-N con columnas clave).

| # | Score | Gen | Variante | Evidencia |
|---|---|---|---|---|
| 1 | **91.5** | `BUB1B` | `c.2210T>G p.L737X` stopgain | ClinVar P/LP (MVA1), InterVar Pathogenic, PVS1+PM2+PP3, HPO match exacto en rabdomiosarcoma |
| 2 | 78.2 | `LZTR1` | `c.2244C>G p.Y748X` stopgain | ClinVar Pathogenic, prácticamente ausente en gnomAD — **hallazgo secundario a considerar** (schwannomatosis/Noonan, dominante) |
| 3 | 72.4 | `BUB1B` | `c.3006T>G p.N1002K` missense | AlphaMissense 0.92, CADD 24.5, AF 6.8e-07 |
| 4 | 65.3 | `GNRHR` | `c.416G>A p.R139H` missense | ClinVar Pathogenic, pero coincidencia fenotípica débil (0.70) |

Las dos variantes de `BUB1B` quedan en #1 y #3, con brecha clara respecto al resto
(mediana del conjunto: 22.0). Ambas con balance alélico limpio (0.54 / 0.46) y GQ 99.

## Limitaciones (importantes)

- **Los pesos están fijados a mano, no aprendidos.** Con un solo caso y sin datos
  etiquetados no se puede ajustar un modelo; son un *prior* razonable, no un ajuste.
- **La validación con `BUB1B` no es independiente.** `BUB1B` se encontró primero con una
  revisión manual de panel, y el diseño del score (en particular la decisión de que
  LOEUF/pLI no penalice recesivos) se tomó sabiendo que el gen candidato era recesivo. Para
  validar de verdad hace falta un conjunto de casos simulados: insertar variantes
  patogénicas conocidas de ClinVar en genomas de fondo (p. ej. 1000 Genomes) junto con los
  HPO de su enfermedad, y medir en qué posición las deja el score (top-1, top-10, rango medio).
  Ese es el siguiente paso real, y es también lo que permitiría entrenar un
  *learning-to-rank* (LightGBM) en vez de pesos fijos.
- **Cobertura**: solo exónico y ±10 pb de splicing. Quedan fuera variantes intrónicas
  profundas y regulatorias, CNVs/SVs, expansiones de repetidos y el mitocondrial (excluido
  en `NORMALIZE`; requiere un pipeline propio por heteroplasmia).
- **La fase no está resuelta**: las dos variantes de `BUB1B` están a ~10.9 kb, demasiado
  lejos para ligarlas con lecturas cortas. Confirmar *trans* requiere muestra de los padres
  (Sanger dirigido bastaría) o lectura larga.
- `S_PHENO` depende de que el gen ya esté ligado a enfermedad en HPO; genes nuevos quedan
  en desventaja por construcción.
