"""
Gene-locus dataset preparation pipeline.

Creates gene-level locus sequences and GO labels in a DNACLIP-style format.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from dataclasses import dataclass
from statistics import median
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.parse import unquote


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


@dataclass
class GeneRecord:
    gene_key: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    source: str
    tax_id: str
    entrez_gene_id: str


def parse_gff3_attrs(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_str.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = unquote(value)
        else:
            attrs[part] = ""
    return attrs


def parse_gtf_attrs(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, value = part.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def infer_annotation_format(path: str) -> Optional[str]:
    base = path[:-3] if path.endswith(".gz") else path
    lower = base.lower()
    if lower.endswith(".gff3") or lower.endswith(".gff"):
        return "gff3"
    if lower.endswith(".gtf"):
        return "gtf"
    return None


def _extract_geneid_from_dbxref(dbxref: str) -> Optional[str]:
    for item in dbxref.split(","):
        if item.startswith("GeneID:"):
            return item.split(":", 1)[1]
    return None


def parse_gff3_genes(
    gff3_path: str,
    tax_id: str,
    chrom_rename: Optional[Dict[str, str]] = None,
) -> List[GeneRecord]:
    genes: List[GeneRecord] = []
    with open_text(gff3_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = fields
            if feature != "gene":
                continue
            if chrom_rename and chrom in chrom_rename:
                chrom = chrom_rename[chrom]
            ad = parse_gff3_attrs(attrs)

            gene_id = None
            dbx = ad.get("Dbxref") or ad.get("db_xref") or ""
            if dbx:
                gene_id = _extract_geneid_from_dbxref(dbx)

            fallback_key = ad.get("ID") or ad.get("Name") or ad.get("gene") or ad.get("locus_tag")
            gene_key = gene_id or fallback_key
            if not gene_key:
                gene_key = f"{chrom}:{start}-{end}:{strand}"

            gene_name = ad.get("Name") or ad.get("gene") or ad.get("locus_tag") or gene_key

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            genes.append(
                GeneRecord(
                    gene_key=gene_key,
                    gene_name=gene_name,
                    chrom=chrom,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    source=source or "gff3",
                    tax_id=tax_id,
                    entrez_gene_id=gene_id or "",
                )
            )
    return genes


def parse_gtf_genes(
    gtf_path: str,
    tax_id: str,
    chrom_rename: Optional[Dict[str, str]] = None,
) -> List[GeneRecord]:
    genes: List[GeneRecord] = []
    with open_text(gtf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = fields
            if feature != "gene":
                continue
            if chrom_rename and chrom in chrom_rename:
                chrom = chrom_rename[chrom]
            ad = parse_gtf_attrs(attrs)

            gene_key = ad.get("gene_id")
            if not gene_key:
                continue

            gene_name = ad.get("gene_name") or gene_key

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            genes.append(
                GeneRecord(
                    gene_key=gene_key,
                    gene_name=gene_name,
                    chrom=chrom,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    source=source or "gtf",
                    tax_id=tax_id,
                    entrez_gene_id="",
                )
            )
    return genes


def load_genome_fasta(fasta_path: str) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    current_seq_name: Optional[str] = None
    seq_chunks: List[str] = []
    with open_text(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq_name is not None:
                    sequences[current_seq_name] = "".join(seq_chunks).upper()
                current_seq_name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if current_seq_name is not None:
            sequences[current_seq_name] = "".join(seq_chunks).upper()
    return sequences


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def extract_loci(
    genes: List[GeneRecord],
    genome_sequences: Dict[str, str],
) -> Tuple[List[Tuple[GeneRecord, str]], int, int]:
    loci: List[Tuple[GeneRecord, str]] = []
    missing_contig = 0
    invalid_coords = 0
    for gene in genes:
        seq = genome_sequences.get(gene.chrom)
        if seq is None:
            missing_contig += 1
            continue
        if gene.start < 1 or gene.end > len(seq) or gene.start > gene.end:
            invalid_coords += 1
            continue
        subseq = seq[gene.start - 1 : gene.end]
        if gene.strand == "-":
            subseq = reverse_complement(subseq)
        loci.append((gene, subseq))
    return loci, missing_contig, invalid_coords


def load_gene2go(
    gene2go_path: str,
    tax_id: str,
    drop_not: bool = True,
) -> Dict[str, Dict[str, set[str]]]:
    gene_go: Dict[str, Dict[str, set[str]]] = {}
    with open_text(gene2go_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            row_tax_id, gene_id, go_id, evidence, qualifier, go_term, pubmed, category = parts[:8]
            if row_tax_id != tax_id:
                continue
            if not go_id.startswith("GO:"):
                continue
            if drop_not and qualifier:
                if "NOT" in qualifier.split("|"):
                    continue
            category_lower = category.strip().lower()
            if category_lower in {"function", "molecular_function", "mf"}:
                aspect = "MF"
            elif category_lower in {"process", "biological_process", "bp"}:
                aspect = "BP"
            elif category_lower in {"component", "cellular_component", "cc"}:
                aspect = "CC"
            else:
                continue

            bucket = gene_go.setdefault(gene_id, {"all": set(), "BP": set(), "MF": set(), "CC": set()})
            bucket["all"].add(go_id)
            bucket[aspect].add(go_id)
    return gene_go


def load_gene2ensembl_map(gene2ensembl_path: str, tax_id: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open_text(gene2ensembl_path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            row_tax_id, gene_id, ensg = parts[0], parts[1], parts[2]
            if row_tax_id != tax_id:
                continue
            if not ensg or ensg == "-":
                continue
            mapping[ensg] = gene_id
    return mapping


def apply_entrez_mapping(genes: List[GeneRecord], mapping: Dict[str, str]) -> None:
    for gene in genes:
        if gene.entrez_gene_id:
            continue
        mapped = mapping.get(gene.gene_key)
        if mapped:
            gene.entrez_gene_id = mapped


def wrap_fasta(seq: str, width: int = 60) -> List[str]:
    return [seq[i : i + width] for i in range(0, len(seq), width)]


def write_genes_tsv(path: str, genes: List[GeneRecord]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "gene_key",
            "gene_name",
            "chrom",
            "start",
            "end",
            "strand",
            "source",
            "tax_id",
            "entrez_gene_id",
        ])
        for gene in genes:
            writer.writerow([
                gene.gene_key,
                gene.gene_name,
                gene.chrom,
                gene.start,
                gene.end,
                gene.strand,
                gene.source,
                gene.tax_id,
                gene.entrez_gene_id,
            ])


def write_loci_outputs(
    out_fasta: str,
    out_tsv: str,
    loci: List[Tuple[GeneRecord, str]],
    species: str,
) -> None:
    with open(out_fasta, "w") as fasta_f, open(out_tsv, "w", newline="") as tsv_f:
        writer = csv.writer(tsv_f, delimiter="\t")
        writer.writerow([
            "gene_key",
            "gene_name",
            "chrom",
            "start",
            "end",
            "strand",
            "length",
            "tax_id",
            "entrez_gene_id",
        ])
        for gene, seq in loci:
            header = f">{gene.gene_key}|{species}|{gene.chrom}:{gene.start}-{gene.end}({gene.strand})|name={gene.gene_name}"
            fasta_f.write(header + "\n")
            for line in wrap_fasta(seq):
                fasta_f.write(line + "\n")
            writer.writerow([
                gene.gene_key,
                gene.gene_name,
                gene.chrom,
                gene.start,
                gene.end,
                gene.strand,
                len(seq),
                gene.tax_id,
                gene.entrez_gene_id,
            ])


def write_labels_go(
    path: str,
    loci: List[Tuple[GeneRecord, str]],
    gene_go: Dict[str, Dict[str, set[str]]],
) -> Tuple[int, int]:
    num_entrez_mapped = 0
    num_go_labeled = 0
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "gene_key",
            "entrez_gene_id",
            "go_terms",
            "go_bp",
            "go_mf",
            "go_cc",
        ])
        for gene, _seq in loci:
            gene_id = gene.entrez_gene_id or (gene.gene_key if gene.gene_key.isdigit() else "")
            if gene_id:
                num_entrez_mapped += 1
            go_bucket = gene_go.get(gene_id, {"all": set(), "BP": set(), "MF": set(), "CC": set()})
            go_terms = ";".join(sorted(go_bucket["all"]))
            go_bp = ";".join(sorted(go_bucket["BP"]))
            go_mf = ";".join(sorted(go_bucket["MF"]))
            go_cc = ";".join(sorted(go_bucket["CC"]))
            if go_terms:
                num_go_labeled += 1
            writer.writerow([
                gene.gene_key,
                gene_id,
                go_terms,
                go_bp,
                go_mf,
                go_cc,
            ])
    return num_entrez_mapped, num_go_labeled


def write_summary(
    path: str,
    species: str,
    tax_id: str,
    genes: List[GeneRecord],
    loci: List[Tuple[GeneRecord, str]],
    missing_contig: int,
    invalid_coords: int,
    num_entrez_mapped: int,
    num_go_labeled: int,
) -> None:
    lengths = [len(seq) for _gene, seq in loci]
    summary = {
        "species": species,
        "tax_id": tax_id,
        "num_genes_total": len(genes),
        "num_loci_extracted": len(loci),
        "num_missing_contig_in_fasta": missing_contig,
        "num_invalid_coordinates": invalid_coords,
        "num_entrez_mapped": num_entrez_mapped,
        "num_go_labeled": num_go_labeled,
        "locus_length_min": min(lengths) if lengths else 0,
        "locus_length_median": median(lengths) if lengths else 0,
        "locus_length_max": max(lengths) if lengths else 0,
    }
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)


def load_chrom_rename_map(path: Optional[str]) -> Optional[Dict[str, str]]:
    if not path:
        return None
    with open(path) as f:
        return json.load(f)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare gene-locus dataset with GO labels.")
    parser.add_argument("--species", required=True, help="Species or dataset label.")
    parser.add_argument("--tax_id", required=True, help="NCBI taxonomy ID.")
    parser.add_argument("--genome_fasta", required=True, help="Genome FASTA (can be .gz).")
    parser.add_argument("--annotation", required=True, help="Gene annotation GFF3/GTF (can be .gz).")
    parser.add_argument("--annotation_format", choices=["gff3", "gtf"], help="Annotation format.")
    parser.add_argument("--gene2go", required=True, help="NCBI gene2go TSV (can be .gz).")
    parser.add_argument("--gene2ensembl", help="NCBI gene2ensembl TSV (optional).")
    parser.add_argument("--chrom_rename_map", help="JSON mapping of contig names (optional).")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--keep_not_qualifier", action="store_true", help="Keep GO annotations with NOT qualifier.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotation_format = args.annotation_format or infer_annotation_format(args.annotation)
    if not annotation_format:
        raise SystemExit("Could not infer annotation format; specify --annotation_format.")

    chrom_rename = load_chrom_rename_map(args.chrom_rename_map)

    if annotation_format == "gff3":
        genes = parse_gff3_genes(args.annotation, args.tax_id, chrom_rename)
    else:
        genes = parse_gtf_genes(args.annotation, args.tax_id, chrom_rename)

    genes.sort(key=lambda g: (g.chrom, g.start, g.end, g.gene_key))

    if args.gene2ensembl:
        mapping = load_gene2ensembl_map(args.gene2ensembl, args.tax_id)
        apply_entrez_mapping(genes, mapping)

    genome = load_genome_fasta(args.genome_fasta)
    loci, missing_contig, invalid_coords = extract_loci(genes, genome)

    gene_go = load_gene2go(args.gene2go, args.tax_id, drop_not=not args.keep_not_qualifier)

    import os
    os.makedirs(args.out_dir, exist_ok=True)

    genes_tsv = os.path.join(args.out_dir, "genes.tsv")
    loci_fasta = os.path.join(args.out_dir, "loci.fasta")
    loci_tsv = os.path.join(args.out_dir, "loci.tsv")
    labels_go_tsv = os.path.join(args.out_dir, "labels_go.tsv")
    summary_json = os.path.join(args.out_dir, "summary.json")

    write_genes_tsv(genes_tsv, genes)
    write_loci_outputs(loci_fasta, loci_tsv, loci, args.species)
    num_entrez_mapped, num_go_labeled = write_labels_go(labels_go_tsv, loci, gene_go)
    write_summary(
        summary_json,
        args.species,
        args.tax_id,
        genes,
        loci,
        missing_contig,
        invalid_coords,
        num_entrez_mapped,
        num_go_labeled,
    )

    print(f"Wrote {len(genes)} genes to {genes_tsv}")
    print(f"Wrote {len(loci)} loci to {loci_fasta} and {loci_tsv}")
    print(f"Wrote GO labels to {labels_go_tsv}")
    print(f"Summary: {summary_json}")


if __name__ == "__main__":
    main()
