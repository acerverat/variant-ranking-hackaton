# Variant annotation for clinical prioritization

`VCF (WGS, GRCh38) → NORMALIZE → ANNOVAR → SELECT_CANDIDATES`

Based on [Pipelines_Inmegen/Annotation](https://github.com/INMEGEN/Pipelines_Inmegen/tree/Principal/Annotation)
(INMEGEN, CC BY-NC-SA 4.0 license: give credit, non-commercial use, share derivatives under the same license).

| Step | What it does | Output |
|---|---|---|
| `NORMALIZE` | `PASS`, splits multiallelics, aligns indels (if `ref_fasta` is given), drops `*`, min GQ/DP, optional BED | `normalized/*.norm.vcf.gz` |
| `ANNOVAR` | `convert2annovar.pl` + `table_annovar.pl` (configurable protocols, `-hgvs`, splicing ±10 bp) | `annovar/*.hg38_multianno.txt` |
| `SELECT_CANDIDATES` | non-synonymous exonic/splicing, gnomAD AF ≤ 1% (ClinVar P/LP ≤ 5%), genotype, per-gene count | `candidates/*.candidates.tsv` |

## Server (gpib)

Paths already set in `nextflow.config`:

- `sample_info.tsv` → `/STORAGE/gpib/proyectos/HuggingFace/data/WGS_EX2312012_HGWCNDSX7.vcf.gz`
- `annovar_dir` → `/scratch/home/acerverat/programs/annovar` (`table_annovar.pl` present, 2026-09 install)
- `ref_fasta` → `/STORAGE/references/Pipelines_INMEGEN/star_v2.7.11b/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  (Ensembl GRCh38, no `chr` prefix; indexed with `samtools faidx`)
- `annovar_protocol` → 9 databases, all verified with a real test run:
  `refGeneWithVer, ensGene, avsnp151, clinvar_20240917, gnomad41_genome, gnomad41_exome, dbnsfp47a, dbscsnv11, intervar_20180118`
  - `avsnp151`, `clinvar_20240917`, `gnomad41_genome` are symlinked from the previous install at
    `/STORAGE/references/Pipelines_INMEGEN/Annovar/humandb` — not re-downloaded.
  - `gnomad41_exome`, `dbnsfp47a`, `dbscsnv11`, `intervar_20180118` were downloaded on 2026-09-28
    (log: `${annovar_dir}/logs/download_missing_dbs.log`).
  - `dbnsfp47a` replaces `dbnsfp33a` (REVEL, AlphaMissense, MetaRNN, etc.), which is why it's no longer in the protocol.
  - `intervar_20180118` adds `InterVar_automated` and the 28 ACMG criteria (`PVS1`, `PS1-4`, `PM1-6`, `PP1-5`,
    `BA1`, `BS1-4`, `BP1-7`) as columns — they pass straight through to `candidates.tsv` with no code changes.
- **`ANNOVAR` does not use `--vcfinput`.** With all 9 databases over the whole genome, the automatic VCF
  regeneration that `table_annovar.pl --vcfinput` does crashes on certain intergenic/ncRNA sites with
  `prefield not defined ... with field=N and prefield=M` (a known ANNOVAR bug). The `.txt` file is still
  generated correctly, and it's the only thing used downstream, so the module converts to avinput with
  `convert2annovar.pl -includeinfo` first and annotates with `--otherinfo` instead of `--vcfinput` — same
  result, without regenerating the VCF.
- **Shared server — capped at half the resources.** `gpib` has 76 cores / 1 TB RAM total, but it's shared
  with other groups. The local `process.executor` is capped at 38 cores / 500 GB for the whole pipeline
  (`ANNOVAR` uses 16 of those cores). Before launching, check `uptime`/`htop`: if the server is already
  loaded by other users, lower the cap further in `nextflow.config`.

## Usage

```bash
cd /STORAGE/gpib/proyectos/HuggingFace/annotation
nextflow run main.nf -stub                    # wiring test, no real data
nextflow run main.nf                          # real run (uses nextflow.config as is)
nextflow run main.nf -profile slurm           # if submitting to a queue instead of local
```

To override a single parameter without editing `nextflow.config`: `--annovar_protocol ...`, `--regions panel.bed`, etc.

Before launching, the pipeline checks that `humandb/hg38_<protocol>.txt` exists for every protocol.
