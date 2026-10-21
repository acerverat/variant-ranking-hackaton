#!/usr/bin/env python3
"""Similitud fenotipo-gen basada en HPO (Resnik + best-match average).

Para cada gen candidato, compara el conjunto de términos HPO del paciente
contra el conjunto de términos HPO asociado a ese gen (genes_to_phenotype.txt
de HPO, que agrega las anotaciones de todas las enfermedades ligadas al gen).

Método (estándar en priorización de variantes, el mismo principio que usan
Exomiser/Phenolyzer):
  1. IC(t) = -log(frecuencia(t)) donde frecuencia(t) = fracción de genes
     anotados (en genes_to_phenotype.txt) a t o a un descendiente de t.
     Términos raros y específicos => IC alto; términos genéricos => IC ~0.
  2. sim(t1, t2) = IC del ancestro común más informativo (MICA) entre t1 y t2
     (similitud de Resnik).
  3. Para cada término del paciente, se toma la mejor coincidencia entre los
     términos del gen (best match), y se promedia (BMA). Esa es la dirección
     clínicamente relevante: "qué tan bien explica este gen lo que tiene el
     paciente". También se reporta la versión simétrica (BMA gen<->paciente).

Solo biblioteca estándar (sin numpy/scipy) para correr en cualquier entorno.
"""
import argparse
import csv
import math
import sys
from collections import defaultdict

ROOT = "HP:0000001"


def parse_obo(path):
    """Devuelve (is_a, alt2primary, names): grafo mínimo de HPO desde hp.obo."""
    is_a = defaultdict(set)
    alt2primary = {}
    names = {}
    cur_id = None
    cur_alts = []
    obsolete = False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur_id, cur_alts, obsolete = None, [], False
            elif line.startswith("id: "):
                cur_id = line[4:].strip()
            elif line.startswith("name: ") and cur_id:
                names[cur_id] = line[6:].strip()
            elif line.startswith("is_a: ") and cur_id:
                parent = line[6:].split("!")[0].strip()
                is_a[cur_id].add(parent)
            elif line.startswith("alt_id: ") and cur_id:
                cur_alts.append(line[8:].strip())
            elif line.startswith("is_obsolete: true"):
                obsolete = True
            elif line == "" and cur_id and not obsolete:
                for a in cur_alts:
                    alt2primary[a] = cur_id
    return is_a, alt2primary, names


def build_ancestors(is_a):
    """Devuelve una función term_id -> set(ancestros incluyendo el término mismo), memoizada."""
    cache = {}

    def climb(term, seen):
        if term in cache:
            return cache[term]
        if term in seen:  # ciclo defensivo, no debería pasar en HPO
            return {term}
        seen = seen | {term}
        acc = {term}
        for parent in is_a.get(term, ()):
            acc |= climb(parent, seen)
        cache[term] = acc
        return acc

    return lambda term: climb(term, frozenset())


def normalize(term_id, alt2primary, is_a):
    t = alt2primary.get(term_id, term_id)
    if t not in is_a and t != ROOT and t not in alt2primary.values():
        return None
    return t


def load_gene_terms(path, alt2primary, is_a, valid_terms):
    gene_terms = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = row.get("gene_symbol")
            term = row.get("hpo_id")
            if not gene or not term:
                continue
            t = alt2primary.get(term, term)
            if t in valid_terms:
                gene_terms[gene].add(t)
    return gene_terms


def compute_ic(gene_terms, ancestors_of):
    """IC(t) = -log(#genes anotados a t o descendiente / #genes totales)."""
    term_gene_count = defaultdict(int)
    for terms in gene_terms.values():
        seen = set()
        for t in terms:
            seen |= ancestors_of(t)
        for a in seen:
            term_gene_count[a] += 1
    n_genes = len(gene_terms)
    ic = {}
    for t, count in term_gene_count.items():
        ic[t] = -math.log(count / n_genes)
    return ic


def resnik(t1, t2, ancestors_of, ic):
    common = ancestors_of(t1) & ancestors_of(t2)
    if not common:
        return 0.0
    return max(ic.get(a, 0.0) for a in common)


def gene_similarity(patient_terms, gene_terms, ancestors_of, ic):
    if not gene_terms:
        return None
    sim_matrix = {(p, g): resnik(p, g, ancestors_of, ic)
                  for p in patient_terms for g in gene_terms}
    p2g = [max(sim_matrix[(p, g)] for g in gene_terms) for p in patient_terms]
    g2p = [max(sim_matrix[(p, g)] for p in patient_terms) for g in gene_terms]
    score_p2g = sum(p2g) / len(p2g)
    score_bma = (sum(p2g) / len(p2g) + sum(g2p) / len(g2p)) / 2
    # mejor pareja (término paciente, término gen) para explicar el score
    best_pair = max(sim_matrix, key=sim_matrix.get)
    return {
        "HPO_SIM_P2G": score_p2g,
        "HPO_SIM_BMA": score_bma,
        "N_GENE_HPO_TERMS": len(gene_terms),
        "HPO_BEST_MATCH": f"{best_pair[0]}~{best_pair[1]} (IC={sim_matrix[best_pair]:.2f})",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True, help="candidates.tsv (debe tener columna GENE)")
    ap.add_argument("--hpo-terms", required=True,
                     help="términos HPO del paciente, separados por coma (HP:xxxxxxx,...)")
    ap.add_argument("--obo", default="refs/hp.obo")
    ap.add_argument("--g2p", default="refs/genes_to_phenotype.txt")
    ap.add_argument("--gene-col", default="GENE")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    sys.stderr.write("Cargando ontología HPO...\n")
    is_a, alt2primary, names = parse_obo(args.obo)
    ancestors_of = build_ancestors(is_a)
    valid_terms = set(is_a.keys()) | {ROOT}

    patient_terms_raw = [t.strip() for t in args.hpo_terms.split(",") if t.strip()]
    patient_terms = []
    for t in patient_terms_raw:
        norm = alt2primary.get(t, t)
        if norm not in valid_terms:
            sys.stderr.write(f"AVISO: término HPO no encontrado en la ontología: {t}\n")
            continue
        patient_terms.append(norm)
    if not patient_terms:
        sys.exit("Ningún término HPO del paciente es válido; nada que comparar.")
    sys.stderr.write("Términos del paciente:\n")
    for t in patient_terms:
        sys.stderr.write(f"  {t}  {names.get(t, '?')}\n")

    sys.stderr.write("Cargando anotaciones gen -> HPO...\n")
    gene_terms_all = load_gene_terms(args.g2p, alt2primary, is_a, valid_terms)
    sys.stderr.write(f"  {len(gene_terms_all)} genes con anotación HPO\n")

    sys.stderr.write("Calculando contenido de información (IC) por término...\n")
    ic = compute_ic(gene_terms_all, ancestors_of)

    with open(args.candidates, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames

    sys.stderr.write(f"Anotando {len(rows)} variantes candidatas...\n")
    cache = {}
    n_no_data = 0
    for row in rows:
        gene = row.get(args.gene_col, "")
        if gene not in cache:
            gt = gene_terms_all.get(gene, set())
            cache[gene] = gene_similarity(patient_terms, gt, ancestors_of, ic)
        result = cache[gene]
        if result is None:
            n_no_data += 1
            row.update({"HPO_SIM_P2G": "NA", "HPO_SIM_BMA": "NA",
                        "N_GENE_HPO_TERMS": 0, "HPO_BEST_MATCH": "."})
        else:
            row.update({k: (f"{v:.4f}" if isinstance(v, float) else v)
                        for k, v in result.items()})

    out_header = header + ["HPO_SIM_P2G", "HPO_SIM_BMA", "N_GENE_HPO_TERMS", "HPO_BEST_MATCH"]
    rows.sort(key=lambda r: float(r["HPO_SIM_P2G"]) if r["HPO_SIM_P2G"] != "NA" else -1,
              reverse=True)
    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    sys.stderr.write(f"Listo: {args.out}\n")
    sys.stderr.write(f"  {len(rows) - n_no_data} variantes con gen anotado en HPO, "
                      f"{n_no_data} sin anotación (gen no ligado a enfermedad conocida)\n")


if __name__ == "__main__":
    main()
