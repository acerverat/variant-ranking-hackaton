#!/usr/bin/env python3
"""Score compuesto para priorizar variantes candidatas (top-N clínico).

Combina cinco ejes de evidencia en un score 0-100, interpretable: cada eje se
reporta como columna aparte y se genera una cadena REASON legible, para que un
analista pueda ver POR QUÉ una variante quedó arriba (requisito clínico: el
ranking tiene que ser auditable, no una caja negra).

    SCORE = 30*S_FENOTIPO + 25*S_EVIDENCIA + 25*S_IMPACTO + 10*S_FRECUENCIA
            + 10*S_HERENCIA  -  penalización_benigna

Decisiones de diseño que importan:

* **Los scores in-silico solo aplican a missense.** REVEL/AlphaMissense/MetaRNN no
  existen para stopgain/frameshift. Promediarlos a ciegas hundiría a las variantes
  truncantes, que son las de mayor impacto. Por eso S_IMPACTO se calcula por clase
  de consecuencia, y la ausencia de scores en una LoF no penaliza.
* **La restricción génica (LOEUF/pLI) NO penaliza a los candidatos recesivos.** Un gen
  recesivo como BUB1B tiene pLI~0 porque la pérdida heterocigota se tolera; usar pLI
  como filtro general descartaría justo ese tipo de diagnóstico. Aquí la restricción
  solo da un bono a la hipótesis dominante (het aislado en gen restringido).
* **Sin padres no se puede confirmar fase ni de novo.** Dos candidatas en el mismo gen
  se marcan como "posible het compuesto" (no confirmado) y un het aislado recibe
  crédito parcial, no el de un hallazgo recesivo confirmado.
* **Gen sin anotación HPO no se anula a 0**, recibe un piso bajo (`--pheno-floor`): un
  gen aún no ligado a enfermedad puede ser un hallazgo nuevo, pero es menos probable
  que explique un cuadro reconocible.
* La normalización del eje fenotípico usa una escala FIJA (`--pheno-ref`), no min-max
  sobre la muestra, para que los scores sean comparables entre corridas.

Solo biblioteca estándar.
"""
import argparse
import csv
import math
import sys
from collections import Counter

# --- clases de consecuencia -------------------------------------------------
LOF_STRONG = {"stopgain", "frameshift deletion", "frameshift insertion", "startloss"}
LOF_WEAKER = {"stoploss"}
INFRAME = {"nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution"}
MISSENSE = {"nonsynonymous SNV"}

ACMG_PATH = ["PVS1", "PS1", "PS2", "PS3", "PS4",
             "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
             "PP1", "PP2", "PP3", "PP4", "PP5"]
ACMG_BENIGN = ["BA1", "BS1", "BS2", "BS3", "BS4",
               "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]

# Genes que generan LoF espurias de forma recurrente con llamado corto contra referencia:
# región HLA y otros loci hiperpolimórficos (el haplotipo real se aparta de la referencia),
# mucinas y genes con repetidos, y genes gigantes/con alelos nulos comunes en población sana.
# NO se eliminan: se marcan y se penalizan, para que un hallazgo real en estos genes siga visible.
NOISY_PREFIX = ("HLA-", "MUC", "KRTAP", "PRAMEF", "NBPF", "ZNF", "OR1", "OR2", "OR4",
                "OR5", "OR6", "OR8", "OR9", "OR10", "OR11", "OR13", "OR51", "OR52",
                "TAS2R", "GOLGA", "SIRPB", "USP17L", "CT47", "NPIP", "FAM90A")
NOISY_EXACT = {"FLG", "FLG2", "HRNR", "AHNAK2", "TTN", "OBSCN", "NEB", "RP1L1",
               "PDE4DIP", "FCGBP", "DSPP", "MST1L", "CDC27", "PLIN4", "MUC16",
               "SYNE1", "SYNE2", "DNAH17", "CCDC168", "GPRIN2", "ANKRD36",
               "ANKRD36C", "LILRB3", "LILRA6", "KIR2DL1", "KIR3DL1", "PRSS1"}


def is_noisy_gene(gene):
    """True si el gen es de los que producen LoF artefactuales de forma recurrente."""
    for g in (gene or "").split(","):          # ANNOVAR puede listar varios genes
        g = g.strip()
        if not g:
            continue
        if g in NOISY_EXACT or g.startswith(NOISY_PREFIX):
            return True
    return False


def allele_balance(row):
    """Fracción del alelo alterno a partir de AD. None si no se puede calcular."""
    ad = (row.get("AD") or "").strip()
    if not ad or ad == ".":
        return None
    try:
        counts = [int(x) for x in ad.split(",")]
    except ValueError:
        return None
    total = sum(counts)
    if total <= 0 or len(counts) < 2:
        return None
    return sum(counts[1:]) / total


def fnum(value):
    """float tolerante: '.', '', 'NA' -> None."""
    if value is None:
        return None
    v = value.strip()
    if v in (".", "", "NA", "na", "None", "-"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def clamp(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


def load_constraint(path):
    """gene -> (LOEUF, pLI). Devuelve {} si no se da archivo."""
    out = {}
    if not path:
        return out
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["gene"]] = (fnum(r.get("LOEUF")), fnum(r.get("pLI")))
    return out


# --- ejes de evidencia ------------------------------------------------------
def score_phenotype(row, ref, floor):
    sim = fnum(row.get("HPO_SIM_P2G"))
    if sim is None:
        return floor, "gen sin anotación HPO"
    s = clamp(sim / ref)
    best = row.get("HPO_BEST_MATCH", ".")
    return s, f"HPO sim={sim:.2f} ({best})"


def clinvar_class(clnsig):
    """Normaliza CLNSIG a una etiqueta gruesa."""
    c = (clnsig or ".").lower()
    if c in (".", "", "na"):
        return None
    if "conflicting" in c:
        return "conflicting"
    # 'pathogenic' aparece dentro de 'likely_pathogenic'; mirar el más fuerte primero
    if "pathogenic" in c and "likely_pathogenic" not in c.replace("pathogenic/likely_pathogenic", ""):
        if c.startswith("pathogenic") or "/likely_pathogenic" in c:
            return "P"
    if "likely_pathogenic" in c:
        return "LP"
    if "pathogenic" in c:
        return "P"
    if "benign" in c and "likely_benign" not in c.replace("benign/likely_benign", ""):
        if c.startswith("benign") or "/likely_benign" in c:
            return "B"
    if "likely_benign" in c:
        return "LB"
    if "benign" in c:
        return "B"
    if "uncertain" in c:
        return "VUS"
    return "other"


def score_evidence(row):
    """ClinVar + InterVar + conteo de criterios ACMG. Devuelve (score, penalización, notas)."""
    notes = []
    cv = clinvar_class(row.get("CLNSIG"))
    iv = (row.get("InterVar_automated") or ".").strip().lower()

    cv_score = {"P": 1.00, "LP": 0.90, "conflicting": 0.40,
                "VUS": 0.30, "LB": 0.10, "B": 0.00}.get(cv)
    iv_score = {"pathogenic": 0.85, "likely pathogenic": 0.70,
                "uncertain significance": 0.30,
                "likely benign": 0.10, "benign": 0.00}.get(iv)

    n_path = sum(1 for c in ACMG_PATH if (fnum(row.get(c)) or 0) >= 1)
    n_ben = sum(1 for c in ACMG_BENIGN if (fnum(row.get(c)) or 0) >= 1)
    hit_path = [c for c in ACMG_PATH if (fnum(row.get(c)) or 0) >= 1]

    parts = [s for s in (cv_score, iv_score) if s is not None]
    base = max(parts) if parts else 0.30            # sin datos -> neutral, no cero
    # pequeño refuerzo por cantidad de criterios ACMG patogénicos disparados
    base = clamp(base + min(0.10, 0.025 * n_path))

    if cv in ("P", "LP"):
        notes.append(f"ClinVar {row.get('CLNSIG')}")
    elif cv == "conflicting":
        notes.append("ClinVar conflicting")
    if iv in ("pathogenic", "likely pathogenic"):
        notes.append(f"InterVar {row.get('InterVar_automated')}")
    if hit_path:
        notes.append("ACMG+: " + ",".join(hit_path))

    # Penalización: evidencia benigna explícita hunde la variante por debajo de un VUS.
    penalty = 0.0
    if cv == "B":
        penalty = 25.0
        notes.append("ClinVar BENIGNO")
    elif cv == "LB":
        penalty = 15.0
        notes.append("ClinVar probablemente benigno")
    elif iv == "benign":
        penalty = 12.0
        notes.append("InterVar benigno")
    elif iv == "likely benign":
        penalty = 7.0
        notes.append("InterVar probablemente benigno")
    if n_ben and cv not in ("P", "LP"):
        penalty += min(6.0, 2.0 * n_ben)

    return base, penalty, notes


def score_impact(row):
    """Severidad por clase de consecuencia; in-silico solo donde aplica."""
    func = (row.get("Func.refGeneWithVer") or "").lower()
    efunc = (row.get("ExonicFunc.refGeneWithVer") or ".").strip()
    notes = []

    # ensemble in-silico (missense)
    ens = []
    for col, scale in (("REVEL_score", 1.0), ("AlphaMissense_score", 1.0),
                       ("MetaRNN_score", 1.0)):
        v = fnum(row.get(col))
        if v is not None:
            ens.append(clamp(v / scale))
    cadd = fnum(row.get("CADD_phred"))
    if cadd is not None:
        ens.append(clamp(cadd / 35.0))
    ens_mean = sum(ens) / len(ens) if ens else None

    # splicing (dbscSNV)
    spl = [v for v in (fnum(row.get("dbscSNV_ADA_SCORE")),
                       fnum(row.get("dbscSNV_RF_SCORE"))) if v is not None]
    spl_max = max(spl) if spl else None

    s = None
    if efunc in LOF_STRONG:
        s = 0.95
        notes.append(f"truncante ({efunc})")
    elif efunc in LOF_WEAKER:
        s = 0.75
        notes.append(f"{efunc}")
    elif "splicing" in func:
        if spl_max is not None:
            s = clamp(0.45 + 0.55 * spl_max)
            notes.append(f"splicing dbscSNV={spl_max:.2f}")
        else:
            s = 0.55
            notes.append("región de splicing (sin dbscSNV)")
    elif efunc in INFRAME:
        s = 0.50 if ens_mean is None else clamp(0.40 + 0.40 * ens_mean)
        notes.append(f"indel en marco ({efunc})")
    elif efunc in MISSENSE:
        if ens_mean is None:
            s = 0.40
            notes.append("missense sin predictores")
        else:
            s = ens_mean
            bits = []
            for lbl, col in (("REVEL", "REVEL_score"), ("AlphaMissense", "AlphaMissense_score"),
                             ("CADD", "CADD_phred")):
                v = fnum(row.get(col))
                if v is not None:
                    bits.append(f"{lbl}={v:g}")
            notes.append("missense " + " ".join(bits))
    else:
        s = 0.35 if ens_mean is None else clamp(0.30 + 0.40 * ens_mean)
        notes.append(f"consecuencia {efunc or func or '?'}")

    # un sitio de splicing con señal fuerte sobre una missense: tomar lo peor de los dos
    if spl_max is not None and efunc in MISSENSE:
        s = max(s, clamp(0.45 + 0.55 * spl_max))
    return clamp(s), notes


def score_frequency(row):
    af = fnum(row.get("MAX_GNOMAD_AF"))
    if af is None or af <= 0:
        return 1.0, "ausente en gnomAD"
    s = clamp((-math.log10(af) - 2.0) / 3.0)
    return s, f"AF={af:.3g}"


def score_inheritance(row, constraint, male):
    """Compatibilidad con un modo de herencia, siendo explícito sobre lo no verificable."""
    chrom = (row.get("CHROM") or "").replace("chr", "")
    zyg = (row.get("ZYGOSITY") or "").lower()
    n_gene = fnum(row.get("N_CANDIDATES_IN_GENE")) or 1
    gene = row.get("GENE", "")
    loeuf, pli = constraint.get(gene, (None, None))
    notes = []

    sex_chrom = chrom in ("X", "Y")
    if sex_chrom and male and zyg == "hom":
        s = 1.00
        notes.append(f"hemicigoto en {chrom} (varón) — compatible con ligado al {chrom}")
    elif zyg == "hom":
        s = 1.00
        notes.append("homocigoto — compatible con recesivo")
    elif n_gene >= 2:
        s = 0.80
        notes.append(f"{int(n_gene)} candidatas en el gen — posible het compuesto (FASE NO CONFIRMADA, faltan padres)")
    else:
        s = 0.40
        notes.append("het aislado — requeriría de novo/dominante (no verificable sin padres)")
        # la restricción génica solo refuerza la hipótesis dominante; nunca penaliza recesivos
        if loeuf is not None and loeuf < 0.6:
            s = clamp(s + 0.25)
            notes.append(f"gen restringido (LOEUF={loeuf:.2f}) — apoya dominante")
        elif pli is not None and pli > 0.9:
            s = clamp(s + 0.20)
            notes.append(f"pLI={pli:.2f} — apoya dominante")
    return s, notes


def load_chrom_skew(path):
    """gen -> set de cromosomas con exceso de hets desviados genoma-completo
    (salida de chrom_ab_scan.py). Vacío si no se da el archivo."""
    flagged = set()
    if not path:
        return flagged
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("FLAGGED") == "1":
                flagged.add(r["CHROM"])
    return flagged


def quality_flags(row, ab_lo, ab_hi, chrom_skew, chrom_skew_factor):
    """Penalizaciones por artefacto probable. Devuelve (penalización, flag_gen, flag_ab, flag_chrom, notas)."""
    penalty = 0.0
    notes = []
    flag_gene = "1" if is_noisy_gene(row.get("GENE", "")) else "0"
    if flag_gene == "1":
        penalty += 12.0
        notes.append("gen hiperpolimórfico/propenso a artefacto")

    ab = allele_balance(row)
    flag_ab = "0"
    chrom = (row.get("CHROM") or "").replace("chr", "")
    flag_chrom = "1" if chrom in chrom_skew else "0"
    zyg = (row.get("ZYGOSITY") or "").lower()
    if ab is not None and zyg == "het":
        if ab >= 0.98 or ab <= 0.02:
            base_pen = 15.0
            reason = "contradice el genotipo het (probable error de llamado)"
        elif ab < ab_lo or ab > ab_hi:
            base_pen = 10.0
            reason = "sesgado para un het"
        else:
            base_pen = 0.0
            reason = None
        if base_pen > 0:
            flag_ab = "1"
            if flag_chrom == "1":
                # el cromosoma completo tiene un exceso de hets desviados (ver chrom_ab_scan.py):
                # podría ser aneuploidía/CNV grande en mosaico y no un artefacto LOCAL de esta
                # variante. No se puede distinguir con certeza de un sesgo técnico de mapeo
                # (los acrocéntricos/ricos en GC también lo producen) -> se reduce, no se anula.
                pen = base_pen * chrom_skew_factor
                notes.append(f"AB={ab:.2f} {reason} — PERO chr{chrom} tiene exceso de hets "
                             f"desviados genoma-completo (posible aneuploidía/CNV grande; "
                             f"penalización reducida a {pen:.1f}, no es necesariamente artefacto local)")
            else:
                pen = base_pen
                notes.append(f"AB={ab:.2f} {reason}")
            penalty += pen
    return penalty, flag_gene, flag_ab, flag_chrom, ("." if ab is None else f"{ab:.2f}"), notes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True,
                    help="TSV con similitud HPO ya anotada (salida de hpo_gene_similarity.py)")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--constraint", default=None, help="refs/gene_constraint.tsv (LOEUF/pLI)")
    ap.add_argument("--top", type=int, default=25, help="cuántas filas imprimir como resumen")
    ap.add_argument("--pheno-ref", type=float, default=2.0,
                    help="escala fija del eje fenotípico (HPO_SIM_P2G que vale 1.0)")
    ap.add_argument("--pheno-floor", type=float, default=0.15,
                    help="valor del eje fenotípico cuando el gen no tiene anotación HPO")
    ap.add_argument("--male", dest="male", action="store_true", default=True,
                    help="tratar X/Y homocigoto como hemicigoto (por defecto: sí)")
    ap.add_argument("--not-male", dest="male", action="store_false")
    ap.add_argument("--w-pheno", type=float, default=30.0)
    ap.add_argument("--w-evidence", type=float, default=25.0)
    ap.add_argument("--w-impact", type=float, default=25.0)
    ap.add_argument("--w-freq", type=float, default=10.0)
    ap.add_argument("--w-inherit", type=float, default=10.0)
    ap.add_argument("--ab-lo", type=float, default=0.25,
                    help="balance alélico mínimo aceptable en un het")
    ap.add_argument("--ab-hi", type=float, default=0.75,
                    help="balance alélico máximo aceptable en un het")
    ap.add_argument("--no-quality-flags", action="store_true",
                    help="no penalizar genes propensos a artefacto ni balance alélico sesgado")
    ap.add_argument("--chrom-ab-summary", default=None,
                    help="salida de chrom_ab_scan.py (genoma completo); reduce la penalización de "
                         "AB en cromosomas con exceso de hets desviados (posible aneuploidía/CNV)")
    ap.add_argument("--chrom-skew-factor", type=float, default=0.3,
                    help="factor que multiplica la penalización de AB en esos cromosomas (0=anula, 1=sin cambio)")
    args = ap.parse_args()

    constraint = load_constraint(args.constraint)
    chrom_skew = load_chrom_skew(args.chrom_ab_summary)
    if args.chrom_ab_summary:
        sys.stderr.write(f"Cromosomas con exceso de hets desviados genoma-completo: "
                         f"{', '.join(sorted(chrom_skew)) or 'ninguno'} "
                         f"(penalización de AB x{args.chrom_skew_factor} ahí)\n")
    if args.constraint:
        sys.stderr.write(f"Restricción génica: {len(constraint)} genes\n")

    with open(args.candidates, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames
    sys.stderr.write(f"Candidatas de entrada: {len(rows)}\n")

    for row in rows:
        s_phe, n_phe = score_phenotype(row, args.pheno_ref, args.pheno_floor)
        s_evi, penalty, n_evi = score_evidence(row)
        s_imp, n_imp = score_impact(row)
        s_frq, n_frq = score_frequency(row)
        s_inh, n_inh = score_inheritance(row, constraint, args.male)

        if args.no_quality_flags:
            q_pen, f_gene, f_ab, f_chrom, ab_str, n_qual = 0.0, "0", "0", "0", ".", []
        else:
            q_pen, f_gene, f_ab, f_chrom, ab_str, n_qual = quality_flags(
                row, args.ab_lo, args.ab_hi, chrom_skew, args.chrom_skew_factor)
        penalty += q_pen

        total = (args.w_pheno * s_phe + args.w_evidence * s_evi
                 + args.w_impact * s_imp + args.w_freq * s_frq
                 + args.w_inherit * s_inh - penalty)
        loeuf, pli = constraint.get(row.get("GENE", ""), (None, None))

        row.update({
            "SCORE": f"{max(0.0, total):.2f}",
            "AB": ab_str,
            "FLAG_NOISY_GENE": f_gene,
            "FLAG_AB": f_ab,
            "FLAG_CHROM_SKEW": f_chrom,
            "S_PHENO": f"{s_phe:.3f}",
            "S_EVIDENCE": f"{s_evi:.3f}",
            "S_IMPACT": f"{s_imp:.3f}",
            "S_FREQ": f"{s_frq:.3f}",
            "S_INHERIT": f"{s_inh:.3f}",
            "PENALTY": f"{penalty:.1f}",
            "LOEUF": "." if loeuf is None else f"{loeuf:.3f}",
            "pLI": "." if pli is None else f"{pli:.3g}",
            "REASON": " | ".join([n_phe] + n_evi + n_imp + [n_frq] + n_inh + n_qual),
        })

    rows.sort(key=lambda r: float(r["SCORE"]), reverse=True)
    for i, row in enumerate(rows, 1):
        row["RANK"] = i

    new_cols = ["RANK", "SCORE", "S_PHENO", "S_EVIDENCE", "S_IMPACT", "S_FREQ",
                "S_INHERIT", "PENALTY", "AB", "FLAG_NOISY_GENE", "FLAG_AB", "FLAG_CHROM_SKEW",
                "LOEUF", "pLI", "REASON"]
    front = ["RANK", "SCORE", "CHROM", "POS", "REF", "ALT", "GENE", "ZYGOSITY",
             "AAChange.refGeneWithVer", "ExonicFunc.refGeneWithVer", "MAX_GNOMAD_AF",
             "CLNSIG", "InterVar_automated", "S_PHENO", "S_EVIDENCE", "S_IMPACT",
             "S_FREQ", "S_INHERIT", "PENALTY", "AB", "FLAG_NOISY_GENE", "FLAG_AB", "FLAG_CHROM_SKEW",
             "LOEUF", "pLI", "REASON"]
    out_header = front + [c for c in (header + new_cols) if c not in front]

    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    sys.stderr.write(f"Escrito: {args.out}\n\n")
    sys.stderr.write(f"TOP {args.top}\n")
    sys.stderr.write(f"{'#':>3} {'SCORE':>6}  {'GEN':<12} {'ZIG':<5} {'CONSECUENCIA':<22} {'CLINVAR':<34} CAMBIO\n")
    for row in rows[:args.top]:
        aa = (row.get("AAChange.refGeneWithVer") or ".").split(",")[0]
        sys.stderr.write(f"{row['RANK']:>3} {row['SCORE']:>6}  {row['GENE'][:12]:<12} "
                         f"{row['ZYGOSITY'][:5]:<5} {(row.get('ExonicFunc.refGeneWithVer') or '.')[:22]:<22} "
                         f"{(row.get('CLNSIG') or '.')[:34]:<34} {aa[:46]}\n")
    cnt = Counter(r["GENE"] for r in rows[:args.top])
    multi = {g: c for g, c in cnt.items() if c > 1}
    if multi:
        sys.stderr.write("\nGenes con >1 candidata en el top (posible het compuesto, sin fase): "
                         + ", ".join(f"{g} x{c}" for g, c in multi.items()) + "\n")


if __name__ == "__main__":
    main()
