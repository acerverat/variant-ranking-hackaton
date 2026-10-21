#!/usr/bin/env python3
"""Selecciona variantes candidatas (codificantes / splicing, raras) de un
*_multianno.txt de ANNOVAR generado con --vcfinput y agrega columnas útiles
para el ranking: genotipo, zigosidad, AF máxima en gnomAD y conteo de
variantes candidatas por gen (posibles heterocigotos compuestos).

Solo usa la biblioteca estándar para que corra en cualquier contenedor.
"""
import argparse
import csv
import re
import sys
from collections import Counter

VCF_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def af_columns(header, extra):
    """Columnas de frecuencia poblacional global de gnomAD (p. ej. gnomad41_genome_AF,
    gnomad312_AF, gnomAD_exome_ALL) más las indicadas por el usuario."""
    cols = [c for c in header
            if c.lower().startswith("gnomad") and re.search(r"_(af|all)$", c, re.I)]
    return cols + [c for c in extra if c in header and c not in cols]


def parse_sample(fmt, sample):
    fields = dict(zip(fmt.split(":"), sample.split(":")))
    gt = fields.get("GT", "./.")
    alleles = re.split(r"[/|]", gt)
    alt = [a for a in alleles if a not in ("0", ".")]
    if not alt:
        zyg = "ref"
    elif len(alleles) == 1:
        zyg = "hemi"
    elif len(alt) == len(alleles) and len(set(alt)) == 1:
        zyg = "hom"
    else:
        zyg = "het"
    return gt, zyg, fields.get("AD", "."), fields.get("DP", "."), fields.get("GQ", ".")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("multianno")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--gene-col", default="Gene.refGene")
    ap.add_argument("--func-col", default="Func.refGene")
    ap.add_argument("--exonicfunc-col", default="ExonicFunc.refGene")
    ap.add_argument("--clinsig-col", default="CLNSIG")
    ap.add_argument("--func-keep", default="exonic,splicing,exonic;splicing",
                    help="valores de Func.* a conservar (separados por coma)")
    ap.add_argument("--max-af", type=float, default=0.01,
                    help="AF máxima en gnomAD (las P/LP de ClinVar se conservan hasta --max-af-clinvar)")
    ap.add_argument("--max-af-clinvar", type=float, default=0.05)
    ap.add_argument("--keep-synonymous", action="store_true")
    ap.add_argument("--extra-af-cols", default="", help="otras columnas de frecuencia a considerar")
    args = ap.parse_args()

    func_keep = set(args.func_keep.split(","))
    csv.field_size_limit(sys.maxsize)

    # errors="replace": algunos campos de texto libre de ANNOVAR (p. ej. CLNDN de ClinVar) traen
    # bytes que no son UTF-8 válido (comillas de Windows-1252 en el volcado original). No son columnas
    # que se usan para filtrar, así que se tolera el byte suelto en vez de tronar todo el pipeline.
    with open(args.multianno, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        other = [i for i, c in enumerate(header) if c.startswith("Otherinfo")]
        if len(other) < len(VCF_COLS):
            sys.exit("No se encontraron suficientes columnas Otherinfo del VCF; ¿se usó --otherinfo o --vcfinput?")
        vcf_idx = other[-len(VCF_COLS):]  # las últimas 10 son la línea original del VCF
        base_header = [c for i, c in enumerate(header) if i not in other]
        base_idx = [i for i in range(len(header)) if i not in other]
        afc = af_columns(header, [c for c in args.extra_af_cols.split(",") if c])
        af_idx = [header.index(c) for c in afc]
        col = {c: header.index(c) for c in (args.gene_col, args.func_col, args.exonicfunc_col)}
        clin_idx = header.index(args.clinsig_col) if args.clinsig_col in header else None

        kept = []
        n_total = 0
        for row in reader:
            n_total += 1
            if row[col[args.func_col]] not in func_keep:
                continue
            efunc = row[col[args.exonicfunc_col]]
            if efunc == "synonymous SNV" and "splicing" not in row[col[args.func_col]] \
                    and not args.keep_synonymous:
                continue
            afs = [to_float(row[i]) for i in af_idx]
            afs = [a for a in afs if a is not None]
            max_af = max(afs) if afs else 0.0
            clin = row[clin_idx] if clin_idx is not None else "."
            is_plp = "pathogenic" in clin.lower() and "conflicting" not in clin.lower()
            if max_af > (args.max_af_clinvar if is_plp else args.max_af):
                continue
            vcf = dict(zip(VCF_COLS, (row[i] for i in vcf_idx)))
            gt, zyg, ad, dp, gq = parse_sample(vcf["FORMAT"], vcf["SAMPLE"])
            if zyg == "ref":
                continue
            kept.append({
                "CHROM": vcf["CHROM"], "POS": vcf["POS"], "REF": vcf["REF"], "ALT": vcf["ALT"],
                "GENE": row[col[args.gene_col]], "GT": gt, "ZYGOSITY": zyg,
                "AD": ad, "DP": dp, "GQ": gq, "QUAL": vcf["QUAL"],
                "MAX_GNOMAD_AF": f"{max_af:.6g}", "CLINVAR_PLP": int(is_plp),
                "_row": [row[i] for i in base_idx],
            })

    per_gene = Counter(k["GENE"] for k in kept)
    front = ["CHROM", "POS", "REF", "ALT", "GENE", "GT", "ZYGOSITY", "AD", "DP", "GQ", "QUAL",
             "MAX_GNOMAD_AF", "CLINVAR_PLP", "N_CANDIDATES_IN_GENE", "RECESSIVE_COMPATIBLE"]
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(front + base_header)
        for k in kept:
            n = per_gene[k["GENE"]]
            # hom/hemi, o ≥2 candidatas en el gen (compuesto posible; sin padres no hay fase)
            k["N_CANDIDATES_IN_GENE"] = n
            k["RECESSIVE_COMPATIBLE"] = int(k["ZYGOSITY"] in ("hom", "hemi") or n >= 2)
            w.writerow([k[c] for c in front] + k["_row"])

    print(f"variantes leídas: {n_total}", file=sys.stderr)
    print(f"candidatas: {len(kept)} en {len(per_gene)} genes", file=sys.stderr)
    print(f"columnas de AF usadas: {', '.join(afc) or 'ninguna'}", file=sys.stderr)


if __name__ == "__main__":
    main()
