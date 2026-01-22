#!/usr/bin/env python3
"""Build GO-enriched GenBank table (no raw DNA) for SAE concept mapping.

Outputs:
- ecoli_go_enriched_table.csv
- stats/feature_type_counts.csv
- stats/normalized_label_counts.csv
- stats/go_counts_by_aspect.csv
- stats/go_support_tiers.csv
- stats/goslim_counts.csv (now includes name/aspect metadata)
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO


def normalize_text(text: str) -> str:
    import re
    t = text.lower()
    t = re.sub(r"[^a-z0-9\s\-]", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def load_rules(path: Path):
    rules = json.loads(path.read_text())
    import re
    low_info = [re.compile(pat) for pat in rules.get("low_info_patterns", [])]
    compiled_rules = []
    for rule in rules.get("rules", []):
        pats = [re.compile(pat) for pat in rule.get("patterns", [])]
        compiled_rules.append({"label": rule["label"], "quality": rule["quality"], "patterns": pats})
    return low_info, compiled_rules


def classify_label(label: str, low_info, rules):
    if not label:
        return "unknown", "low"
    norm = normalize_text(label)
    for pat in low_info:
        if pat.search(norm):
            return "low_info", "low"
    for rule in rules:
        if any(p.search(norm) for p in rule["patterns"]):
            return rule["label"], rule["quality"]
    return "other", "low"


def load_go_map(go_annot_csv: Path):
    # geneid -> dict of sets
    mapping = defaultdict(lambda: {"go_terms": set(), "go_bp": set(), "go_mf": set(), "go_cc": set()})
    with go_annot_csv.open() as f:
        r = csv.DictReader(f)
        for row in r:
            geneid = row.get("geneid", "")
            if not geneid:
                continue
            for field in ["go_terms", "go_bp", "go_mf", "go_cc"]:
                terms = row.get(field, "")
                if not terms:
                    continue
                for t in terms.split(";"):
                    if t:
                        mapping[geneid][field].add(t)
    return mapping


def load_goslim_map(path: Path):
    mapping = defaultdict(set)
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "," in line:
                parts = line.split(",")
            else:
                parts = line.split("\t")
            if len(parts) < 2:
                continue
            go_id = parts[0].strip()
            slim_id = parts[1].strip()
            if go_id and slim_id:
                mapping[go_id].add(slim_id)
    return mapping


def load_go_obo_terms(path: Path):
    # Returns go_id -> rich metadata from go-basic.obo
    if not path or not path.exists():
        return {}
    terms = {}
    current = {}
    in_term = False
    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if line == "[Term]":
                if in_term and current.get("id") and not current.get("is_obsolete"):
                    terms[current["id"]] = current
                current = {}
                in_term = True
                continue
            if not in_term:
                continue
            if line == "":
                if current.get("id") and not current.get("is_obsolete"):
                    terms[current["id"]] = current
                current = {}
                in_term = False
                continue
            if line.startswith("id:"):
                current["id"] = line.split("id:", 1)[1].strip()
            elif line.startswith("name:"):
                current["name"] = line.split("name:", 1)[1].strip()
            elif line.startswith("namespace:"):
                current["namespace"] = line.split("namespace:", 1)[1].strip()
            elif line.startswith("def:"):
                start = line.find('"')
                end = line.rfind('"')
                if start != -1 and end != -1 and end > start:
                    current["def"] = line[start + 1 : end]
            elif line.startswith("synonym:"):
                start = line.find('"')
                end = line.rfind('"')
                if start != -1 and end != -1 and end > start:
                    current.setdefault("synonyms", []).append(line[start + 1 : end])
            elif line.startswith("is_a:"):
                token = line.split("is_a:", 1)[1].strip().split()[0]
                if token:
                    current.setdefault("is_a", []).append(token)
            elif line.startswith("alt_id:"):
                token = line.split("alt_id:", 1)[1].strip()
                if token:
                    current.setdefault("alt_id", []).append(token)
            elif line.startswith("subset:"):
                token = line.split("subset:", 1)[1].strip()
                if token:
                    current.setdefault("subset", []).append(token)
            elif line.startswith("xref:"):
                token = line.split("xref:", 1)[1].strip()
                if token:
                    current.setdefault("xref", []).append(token)
            elif line.startswith("comment:"):
                current["comment"] = line.split("comment:", 1)[1].strip()
            elif line.startswith("is_obsolete:"):
                current["is_obsolete"] = line.split("is_obsolete:", 1)[1].strip() == "true"
    if in_term and current.get("id") and not current.get("is_obsolete"):
        terms[current["id"]] = current
    return terms


def extract_geneid(feature) -> str | None:
    xrefs = feature.qualifiers.get("db_xref", [])
    for x in xrefs:
        if x.startswith("GeneID:"):
            return x.split(":", 1)[1]
    return None


def feature_start_end_strand(feature):
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else ""
    return start, end, strand


def make_seq_id(
    record_id: str,
    feature_type: str,
    locus_tag: str,
    geneid: str,
    protein_id: str,
    start: int,
    end: int,
    strand: str,
):
    # Prefer stable identifiers; fall back to genomic coordinates.
    if feature_type == "CDS":
        if protein_id:
            return protein_id
        if locus_tag:
            return locus_tag
    elif feature_type == "gene":
        if locus_tag:
            return locus_tag
        if geneid:
            return geneid
    if locus_tag:
        return locus_tag
    if geneid:
        return geneid
    if protein_id:
        return protein_id
    return f"{record_id}:{feature_type}:{start}-{end}:{strand or '.'}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genbank", required=True)
    parser.add_argument("--go-annot", required=True)
    parser.add_argument("--goslim-map", required=True)
    parser.add_argument("--rules", required=True)
    default_obo = (Path(__file__).resolve().parents[1] / "solve_problem_data_codes/data_downloads/go-basic.obo")
    parser.add_argument("--go-obo", default=str(default_obo))
    parser.add_argument("--out", default="ecoli_go_enriched_table.csv")
    parser.add_argument("--stats-dir", default="stats")
    args = parser.parse_args()

    low_info, rules = load_rules(Path(args.rules))
    go_map = load_go_map(Path(args.go_annot))
    goslim_map = load_goslim_map(Path(args.goslim_map))
    go_terms_map = load_go_obo_terms(Path(args.go_obo)) if args.go_obo else {}
    namespace_to_aspect = {
        "biological_process": "BP",
        "molecular_function": "MF",
        "cellular_component": "CC",
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    stats_dir = Path(args.stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)

    # stats
    feature_counts = Counter()
    normalized_counts = Counter()
    go_counts_by_aspect = defaultdict(Counter)
    go_support = Counter()
    goslim_counts = Counter()

    record = SeqIO.read(args.genbank, "genbank")

    # Index gene/CDS features by locus_tag to avoid duplicated gene/CDS rows.
    genes_by_locus = defaultdict(list)
    cds_by_locus = defaultdict(list)
    other_features = []
    rna_locus_tags = set()
    rna_coords = set()
    for feature in record.features:
        if feature.type not in {"gene", "CDS", "rRNA", "tRNA", "ncRNA", "mobile_element", "misc_feature"}:
            continue
        locus_tag = ";".join(feature.qualifiers.get("locus_tag", []))
        if feature.type in {"tRNA", "rRNA", "ncRNA"}:
            if locus_tag:
                rna_locus_tags.add(locus_tag)
            start, end, strand = feature_start_end_strand(feature)
            rna_coords.add((start, end, strand))
        if feature.type == "gene":
            if locus_tag:
                genes_by_locus[locus_tag].append(feature)
            else:
                other_features.append(feature)
        elif feature.type == "CDS":
            if locus_tag:
                cds_by_locus[locus_tag].append(feature)
            else:
                other_features.append(feature)
        else:
            other_features.append(feature)

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "seq_id",
                "feature_type",
                "start",
                "end",
                "strand",
                "geneid",
                "locus_tag",
                "gene",
                "protein_id",
                "product",
                "normalized_label",
                "normalized_quality",
                "go_terms",
                "go_bp",
                "go_mf",
                "go_cc",
                "goslim_terms",
            ]
        )

        used_seq_ids = Counter()

        def is_gene_rna_duplicate(feature) -> bool:
            locus_tag = ";".join(feature.qualifiers.get("locus_tag", []))
            if locus_tag and locus_tag in rna_locus_tags:
                return True
            start, end, strand = feature_start_end_strand(feature)
            return (start, end, strand) in rna_coords

        def write_row(feature, feature_type_override=None, gene_fallback=None):
            feature_type = feature_type_override or feature.type
            geneid = extract_geneid(feature) or ""
            locus_tag = ";".join(feature.qualifiers.get("locus_tag", []))
            gene = ";".join(feature.qualifiers.get("gene", []))
            protein_id = ";".join(feature.qualifiers.get("protein_id", []))
            product = ";".join(feature.qualifiers.get("product", []))

            if gene_fallback is not None:
                gene = gene or gene_fallback.get("gene", "")
                if not product:
                    product = gene_fallback.get("product", "")

            normalized_label, quality = classify_label(product, low_info, rules)

            go_terms = go_map.get(geneid, {}).get("go_terms", set())
            go_bp = go_map.get(geneid, {}).get("go_bp", set())
            go_mf = go_map.get(geneid, {}).get("go_mf", set())
            go_cc = go_map.get(geneid, {}).get("go_cc", set())

            goslim_terms = set()
            for t in go_terms:
                goslim_terms.update(goslim_map.get(t, []))

            start, end, strand = feature_start_end_strand(feature)
            seq_id = make_seq_id(record.id, feature_type, locus_tag, geneid, protein_id, start, end, strand)
            used_seq_ids[seq_id] += 1
            if used_seq_ids[seq_id] > 1:
                seq_id = f"{seq_id}|{used_seq_ids[seq_id]}"

            feature_counts[feature_type] += 1
            normalized_counts[normalized_label] += 1
            for t in go_terms:
                go_support[t] += 1
            for t in go_bp:
                go_counts_by_aspect["BP"][t] += 1
            for t in go_mf:
                go_counts_by_aspect["MF"][t] += 1
            for t in go_cc:
                go_counts_by_aspect["CC"][t] += 1
            for t in goslim_terms:
                goslim_counts[t] += 1

            writer.writerow(
                [
                    seq_id,
                    feature_type,
                    start,
                    end,
                    strand,
                    geneid,
                    locus_tag,
                    gene,
                    protein_id,
                    product,
                    normalized_label,
                    quality,
                    ";".join(sorted(go_terms)),
                    ";".join(sorted(go_bp)),
                    ";".join(sorted(go_mf)),
                    ";".join(sorted(go_cc)),
                    ";".join(sorted(goslim_terms)),
                ]
            )

        # 1) CDS + gene: prefer CDS rows, enrich with gene fields
        all_locus = sorted(set(cds_by_locus.keys()) | set(genes_by_locus.keys()))
        for locus_tag in all_locus:
            gene_features = genes_by_locus.get(locus_tag, [])
            cds_features = cds_by_locus.get(locus_tag, [])
            gene_fallback = None
            if gene_features:
                gene_feature = gene_features[0]
                gene_fallback = {
                    "gene": ";".join(gene_feature.qualifiers.get("gene", [])),
                    "product": ";".join(gene_feature.qualifiers.get("product", [])),
                }
            if cds_features:
                for cds_feature in cds_features:
                    write_row(cds_feature, feature_type_override="CDS", gene_fallback=gene_fallback)
            else:
                # Gene-only (no CDS)
                for gene_feature in gene_features:
                    if is_gene_rna_duplicate(gene_feature):
                        continue
                    write_row(gene_feature, feature_type_override="gene")

        # 2) Other feature types (and gene/CDS without locus_tag)
        for feature in other_features:
            if feature.type == "gene" and is_gene_rna_duplicate(feature):
                continue
            write_row(feature)
    # write stats
    with (stats_dir / "feature_type_counts.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["feature_type", "count"])
        for k, v in feature_counts.most_common():
            writer.writerow([k, v])

    with (stats_dir / "normalized_label_counts.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["normalized_label", "count"])
        for k, v in normalized_counts.most_common():
            writer.writerow([k, v])

    with (stats_dir / "go_counts_by_aspect.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aspect", "go_id", "count"])
        for aspect, counter in go_counts_by_aspect.items():
            for k, v in counter.most_common():
                writer.writerow([aspect, k, v])

    # support tiers
    tiers = [20, 50, 100, 200, 500, 1000]
    with (stats_dir / "go_support_tiers.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["tier", "num_terms"])
        for t in tiers:
            writer.writerow([t, sum(1 for c in go_support.values() if c >= t)])

    with (stats_dir / "goslim_counts.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "goslim_id",
                "goslim_name",
                "aspect",
                "definition",
                "synonyms",
                "is_a",
                "alt_id",
                "subset",
                "xref",
                "comment",
                "count",
            ]
        )
        for k, v in goslim_counts.most_common():
            meta = go_terms_map.get(k, {})
            aspect = namespace_to_aspect.get(meta.get("namespace", ""), "")
            writer.writerow(
                [
                    k,
                    meta.get("name", ""),
                    aspect,
                    meta.get("def", ""),
                    ";".join(meta.get("synonyms", [])),
                    ";".join(meta.get("is_a", [])),
                    ";".join(meta.get("alt_id", [])),
                    ";".join(meta.get("subset", [])),
                    ";".join(meta.get("xref", [])),
                    meta.get("comment", ""),
                    v,
                ]
            )

    print(f"Wrote table: {out_path}")
    print(f"Wrote stats under: {stats_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
