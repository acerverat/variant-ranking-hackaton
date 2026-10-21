#!/usr/bin/env python3
"""Detección de exceso de desbalance alélico por cromosoma (genoma completo).

Objetivo: distinguir un sesgo de AB *local* (artefacto de mapeo en un gen
puntual — HLA, SERPINA1, etc., ya cubierto por FLAG_NOISY_GENE) de un sesgo
*compartido por miles de variantes en el mismo cromosoma*, que es la firma
esperada de una aneuploidía/CNV grande en mosaico (igual principio que la
B-allele frequency de un arreglo citogenómico). Para esto último se necesita
poder estadístico genoma-completo: dentro de las 640 candidatas codificantes
un cromosoma entero puede tener 0-2 variantes, insuficiente para decidir nada.

Método:
  1. Para cada cromosoma, entre los SNV heterocigotos PASS con profundidad
     suficiente, se calcula qué fracción cae fuera de una ventana de AB
     "normal" (por defecto 0.35-0.65).
  2. Se construye una base de referencia robusta (mediana + MAD) a partir de
     esa tasa en los autosomas — robusta a que 1-3 cromosomas ya sean outliers
     (p. ej. los acrocéntricos), porque la mediana no se mueve por unos pocos.
  3. Un cromosoma se marca si su tasa excede la referencia por un z-score
     robusto >= --z-threshold Y el número absoluto de variantes desviadas
     alcanza --min-skewed (poder estadístico mínimo).

LIMITACIÓN IMPORTANTE, por diseño no resuelta aquí: este estadístico por sí
solo NO distingue una aneuploidía mosaica real de un sesgo técnico de mapeo
que también es específico de cromosoma (los acrocéntricos y ricos en GC
—13,14,15,21,22— son conocidos por tener más profundidad Y más desbalance con
lectura corta, por razones puramente técnicas). Por eso esta señal solo debe
REDUCIR la penalización de calidad en composite_score.py, nunca convertirse en
"evidencia positiva" de aneuploidía por sí misma. Para confirmar de verdad se
necesita: comparar contra el mismo pipeline corrido en muestras control, un
llamador de CNV dedicado (profundidad en bins + BAF), o cariotipo/FISH.

Solo biblioteca estándar; usa bcftools vía subprocess (ya está en el PATH).
"""
import argparse
import statistics
import subprocess
import sys

AUTOSOMES = [str(i) for i in range(1, 23)]
SEX = ["X", "Y"]


def robust_z(x, values):
    med = statistics.median(values)
    mad = statistics.median([abs(v - med) for v in values]) * 1.4826  # ~ SD bajo normalidad
    if mad == 0:
        return 0.0 if x == med else float("inf")
    return (x - med) / mad


def scan(vcf, min_dp, ab_lo, ab_hi):
    """bcftools query en streaming; un pase, agregación en memoria por cromosoma."""
    cmd = ["bcftools", "query", "-f", "%CHROM\t[%AD]\n",
           "-i", f'TYPE="snp" && GT="het" && FMT/DP>={min_dp}', vcf]
    n_het = {}
    n_skew = {}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        chrom, ad = line.rstrip("\n").split("\t")
        try:
            r, a = ad.split(",")[:2]
            r, a = int(r), int(a)
        except (ValueError, IndexError):
            continue
        tot = r + a
        if tot < min_dp:
            continue
        ab = a / tot
        n_het[chrom] = n_het.get(chrom, 0) + 1
        if ab < ab_lo or ab > ab_hi:
            n_skew[chrom] = n_skew.get(chrom, 0) + 1
    proc.wait()
    if proc.returncode != 0:
        sys.exit(f"bcftools terminó con código {proc.returncode}")
    return n_het, n_skew


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="VCF normalizado, GENOMA COMPLETO (no el de candidatas)")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--min-dp", type=int, default=20)
    ap.add_argument("--ab-lo", type=float, default=0.35, help="ventana de AB 'normal' (límite inferior)")
    ap.add_argument("--ab-hi", type=float, default=0.65)
    ap.add_argument("--z-threshold", type=float, default=3.0,
                    help="z-score robusto mínimo (vs. mediana de autosomas) para marcar un cromosoma")
    ap.add_argument("--min-skewed", type=int, default=500,
                    help="mínimo de variantes desviadas (poder estadístico) para poder marcar")
    args = ap.parse_args()

    sys.stderr.write(f"Escaneando {args.vcf} (DP>={args.min_dp}, ventana AB {args.ab_lo}-{args.ab_hi})...\n")
    n_het, n_skew = scan(args.vcf, args.min_dp, args.ab_lo, args.ab_hi)

    rate = {c: n_skew.get(c, 0) / n_het[c] for c in n_het if n_het[c] > 0}
    auto_rates = [rate[c] for c in AUTOSOMES if c in rate]
    if len(auto_rates) < 10:
        sys.exit("Muy pocos autosomas con datos; ¿el VCF es genoma completo y normalizado?")

    rows = []
    for c in AUTOSOMES + SEX:
        if c not in rate:
            continue
        z = robust_z(rate[c], auto_rates)
        flagged = (c in AUTOSOMES and z >= args.z_threshold and n_skew.get(c, 0) >= args.min_skewed)
        rows.append({
            "CHROM": c, "N_HET": n_het[c], "N_SKEWED": n_skew.get(c, 0),
            "SKEW_RATE": rate[c], "ROBUST_Z": z, "FLAGGED": int(flagged),
        })

    med = statistics.median(auto_rates)
    sys.stderr.write(f"Mediana de tasa de sesgo en autosomas: {med:.3f}\n\n")
    sys.stderr.write(f"{'CHR':<4}{'N_HET':>9}{'N_SESGADOS':>11}{'TASA':>8}{'Z':>7}{'MARCADO':>9}\n")
    for r in rows:
        sys.stderr.write(f"{r['CHROM']:<4}{r['N_HET']:>9}{r['N_SKEWED']:>11}"
                         f"{r['SKEW_RATE']:>8.3f}{r['ROBUST_Z']:>7.2f}{'  SI' if r['FLAGGED'] else '':>9}\n")

    flagged_chroms = [r["CHROM"] for r in rows if r["FLAGGED"]]
    sys.stderr.write(f"\nCromosomas marcados (posible señal de CNV/aneuploidía genoma-completo, "
                     f"NO necesariamente real — ver limitaciones en el docstring): "
                     f"{', '.join(flagged_chroms) or 'ninguno'}\n")

    with open(args.out, "w") as out:
        out.write("CHROM\tN_HET\tN_SKEWED\tSKEW_RATE\tROBUST_Z\tFLAGGED\n")
        for r in rows:
            out.write(f"{r['CHROM']}\t{r['N_HET']}\t{r['N_SKEWED']}\t"
                     f"{r['SKEW_RATE']:.4f}\t{r['ROBUST_Z']:.2f}\t{r['FLAGGED']}\n")
    sys.stderr.write(f"\nEscrito: {args.out}\n")


if __name__ == "__main__":
    main()
