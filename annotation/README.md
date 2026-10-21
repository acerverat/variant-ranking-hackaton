# Anotación de variantes para priorización clínica

`VCF (WGS, GRCh38) → NORMALIZE → ANNOVAR → SELECT_CANDIDATES`

Basado en [Pipelines_Inmegen/Annotation](https://github.com/INMEGEN/Pipelines_Inmegen/tree/Principal/Annotation)
(INMEGEN, licencia CC BY-NC-SA 4.0: dar crédito, uso no comercial, compartir derivados con la misma licencia).

| Paso | Qué hace | Salida |
|---|---|---|
| `NORMALIZE` | `PASS`, separa multialélicos, alinea indels (si hay `ref_fasta`), quita `*`, GQ/DP mínimos, BED opcional | `normalized/*.norm.vcf.gz` |
| `ANNOVAR` | `convert2annovar.pl` + `table_annovar.pl` (protocolos configurables, `-hgvs`, splicing ±10 pb) | `annovar/*.hg38_multianno.txt` |
| `SELECT_CANDIDATES` | exónicas/splicing no sinónimas, AF gnomAD ≤ 1 % (P/LP de ClinVar ≤ 5 %), genotipo, conteo por gen | `candidates/*.candidates.tsv` |

## Servidor (gpib)

Rutas ya fijadas en `nextflow.config`:

- `sample_info.tsv` → `/STORAGE/gpib/proyectos/HuggingFace/data/WGS_EX2312012_HGWCNDSX7.vcf.gz`
- `annovar_dir` → `/scratch/home/acerverat/programs/annovar` (`table_annovar.pl` presente, instalación 2026-09)
- `ref_fasta` → `/STORAGE/references/Pipelines_INMEGEN/star_v2.7.11b/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  (Ensembl GRCh38, sin `chr`; indexada con `samtools faidx`)
- `annovar_protocol` → 9 bases, todas verificadas con una corrida real de prueba:
  `refGeneWithVer, ensGene, avsnp151, clinvar_20240917, gnomad41_genome, gnomad41_exome, dbnsfp47a, dbscsnv11, intervar_20180118`
  - `avsnp151`, `clinvar_20240917`, `gnomad41_genome` están enlazadas (symlink) desde la instalación anterior en
    `/STORAGE/references/Pipelines_INMEGEN/Annovar/humandb` — no se volvieron a descargar.
  - `gnomad41_exome`, `dbnsfp47a`, `dbscsnv11`, `intervar_20180118` se descargaron el 2026-09-28
    (log: `${annovar_dir}/logs/download_missing_dbs.log`).
  - `dbnsfp47a` reemplaza a `dbnsfp33a` (REVEL, AlphaMissense, MetaRNN, etc.), por eso ya no está en el protocolo.
  - `intervar_20180118` agrega `InterVar_automated` y los 28 criterios ACMG (`PVS1`, `PS1-4`, `PM1-6`, `PP1-5`,
    `BA1`, `BS1-4`, `BP1-7`) como columnas — pasan directo a `candidates.tsv` sin cambios en el código.
- **`ANNOVAR` no usa `--vcfinput`.** Con las 9 bases sobre el genoma completo, la regeneración automática de
  VCF que hace `table_annovar.pl --vcfinput` revienta en ciertos sitios intergénicos/ncRNA con
  `prefield not defined ... with field=N and prefield=M` (bug conocido de ANNOVAR). El `.txt` se genera bien
  igual, y es lo único que se usa después, así que el módulo convierte a avinput con `convert2annovar.pl
  -includeinfo` primero y anota con `--otherinfo` en vez de `--vcfinput` — mismo resultado, sin regenerar VCF.
- **Servidor compartido — tope de la mitad de los recursos.** `gpib` tiene 76 núcleos / 1 TB RAM en total,
  pero es compartido con otros grupos. `process.executor` local está acotado a 38 núcleos / 500 GB para todo
  el pipeline en conjunto (`ANNOVAR` usa 16 de esos núcleos). Antes de lanzar, revisar `uptime`/`htop`: si el
  servidor ya está cargado por otros usuarios, bajar el tope aún más en `nextflow.config`.

## Uso

```bash
cd /STORAGE/gpib/proyectos/HuggingFace/annotation
nextflow run main.nf -stub                    # prueba de conexiones sin datos
nextflow run main.nf                          # corrida real (usa nextflow.config tal cual)
nextflow run main.nf -profile slurm           # si se manda a la cola en vez de local
```

Para sobrescribir un parámetro puntual sin editar `nextflow.config`: `--annovar_protocol ...`, `--regions panel.bed`, etc.

Antes de lanzar, el pipeline verifica que exista `humandb/hg38_<protocolo>.txt` para cada protocolo.
