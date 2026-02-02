# DNACLIP-style Gene-Locus Data Pipeline (for SAE Study)

This document specifies a **DNACLIP-style** data preparation pipeline to build a **gene-locus DNA dataset** (gene start→end, strand-aware) and attach **GO concept labels** for downstream **SAE** experiments. d
Target organisms for initial validation: **E. coli K-12 MG1655** and **Mus musculus (mouse, GRCm39)**.
You can extend the same pipeline to other organisms by swapping input FASTA/annotation files and choosing the correct GO/ID mapping strategy.

---

## 0) Goal and Definitions

### What “gene-locus” means (DNACLIP-style)

For each **gene**, represent it by the **contiguous genomic locus** on the reference assembly:

- DNA sequence = `reference_genome[gene_start : gene_end]` (1-based inclusive coordinates in GFF/GTF)
- If gene is on the **minus strand**, take the **reverse complement**.
- This locus includes **exons + introns + UTRs** in eukaryotes (contiguous span), and is typically a single contiguous region in prokaryotes.

### What “instance” means

An **instance** is one gene identifier:

- For NCBI RefSeq GFF3: use **Entrez GeneID** when available (`Dbxref=GeneID:...`)
- For Ensembl GTF: use **Ensembl gene_id** (e.g., `ENSMUSG...`) and optionally map to **Entrez GeneID** for GO labels.

---

## 1) Required Data Inputs

### 1.1 Reference genome sequence

- **FASTA** for the reference assembly (chromosome/contig sequences).

### 1.2 Gene annotation

One of:

- **GFF3** (common for NCBI/RefSeq genomes)
- **GTF** (common for Ensembl genomes)

We only need **gene-level records** (`feature == "gene"` in GFF3; `feature == "gene"` in GTF).

### 1.3 GO concept labels (for SAE evaluation / grounding)

Use **NCBI gene2go**:

- `gene2go.tsv` provides `tax_id`, `GeneID`, `GO_ID`, `Aspect (BP/MF/CC)`, evidence, etc.

If your genes are **Ensembl IDs** (mouse), also download:

- `gene2ensembl.tsv` to map **Ensembl gene ID → Entrez GeneID**.

> Notes:
>
> - GO labels are **gene-level**, not nucleotide-level.
> - You can keep full GO or compress to **GO-slim** (optional post-processing step).

---

## 2) Download Commands (E. coli, mouse, GO files)

These are **examples** that you can adapt. Prefer pinning to a specific build/release in your project.

### 2.1 E. coli K-12 MG1655 (NCBI RefSeq)

Download FASTA + GFF3 (gzip), then decompress:

```bash
mkdir -p data/ecoli
ECOLI_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"

curl -L "${ECOLI_BASE}/GCF_000005845.2_ASM584v2_genomic.fna.gz" -o data/ecoli/ecoli.fna.gz
curl -L "${ECOLI_BASE}/GCF_000005845.2_ASM584v2_genomic.gff.gz" -o data/ecoli/ecoli.gff.gz

gunzip -c data/ecoli/ecoli.fna.gz > data/ecoli/ecoli.fna
gunzip -c data/ecoli/ecoli.gff.gz > data/ecoli/ecoli.gff3
```

### 2.2 Mouse (Ensembl GRCm39, example release)

Download FASTA + GTF (gzip), then decompress:

```bash
mkdir -p data/mouse
REL=115
EN_BASE="https://ftp.ensembl.org/pub/release-${REL}"

curl -L "${EN_BASE}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
  -o data/mouse/mouse.fa.gz

curl -L "${EN_BASE}/gtf/mus_musculus/Mus_musculus.GRCm39.${REL}.gtf.gz" \
  -o data/mouse/mouse.gtf.gz

gunzip -c data/mouse/mouse.fa.gz  > data/mouse/mouse.fa
gunzip -c data/mouse/mouse.gtf.gz > data/mouse/mouse.gtf
```

### 2.3 GO label files (NCBI)

```bash
mkdir -p data/ncbi

curl -L "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"      -o data/ncbi/gene2go.gz
curl -L "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz" -o data/ncbi/gene2ensembl.gz

gunzip -c data/ncbi/gene2go.gz      > data/ncbi/gene2go.tsv
gunzip -c data/ncbi/gene2ensembl.gz > data/ncbi/gene2ensembl.tsv
```

---

## 3) Pipeline Steps (Implementation Spec)

### Step A — Parse gene coordinates (gene-level only)

**Inputs**

- `annotation.gff3` OR `annotation.gtf`

**Output**

- `genes.tsv` containing (at minimum):

| column         | meaning                                                             |
| -------------- | ------------------------------------------------------------------- |
| gene_key       | primary gene identifier in this dataset (GeneID or Ensembl gene_id) |
| gene_name      | symbol/name if available                                            |
| chrom          | contig/chromosome id (must match FASTA headers)                     |
| start          | 1-based inclusive start                                             |
| end            | 1-based inclusive end                                               |
| strand         | `+` or `-`                                                      |
| source         | annotation source tag                                               |
| tax_id         | NCBI taxonomy id (optional but recommended)                         |
| entrez_gene_id | resolved Entrez GeneID (optional but recommended)                   |

**GFF3 parsing rules (NCBI)**

- Keep rows with `feature == "gene"`
- Parse attribute column (9th column) into dict (split by `;`, then `=`)
- Gene identifier preference:
  1) `Dbxref` contains `GeneID:<id>` → use `<id>` as `gene_key` and `entrez_gene_id`
  2) Else fallback to `ID` or `Name`

**GTF parsing rules (Ensembl)**

- Keep rows with `feature == "gene"`
- Parse attribute field: key-value pairs (e.g., `gene_id "ENSMUSG..."`)
- Use `gene_id` as `gene_key`
- `gene_name` from `gene_name` if present
- Keep (chrom, start, end, strand)

**Common pitfall**

- Chromosome naming must match between FASTA and annotation (e.g., `1` vs `chr1`). Provide a renaming map if needed.

### Step B — Load genome FASTA

**Input**

- `genome.fa` (FASTA)

**Output**

- `genome_dict[chrom] = sequence_string` (uppercase recommended)

**Pitfalls**

- FASTA headers: take the first token before spaces as the contig key.
- Some genomes have many contigs; keep memory usage in mind.

### Step C — Extract gene-locus sequence (start→end, strand-aware)

**Input**

- `genes.tsv` + `genome.fa`

**Output**

- `loci.fasta` (one record per gene instance)
- `loci.tsv` (metadata including length)

**Extraction rules**

- Coordinates are 1-based inclusive in GFF/GTF.
- Python slicing uses 0-based half-open: `seq[start-1 : end]`.
- If `strand == '-'`: reverse-complement the extracted locus.

**FASTA record naming**
Use a stable structured header, e.g.:

```
>{gene_key}|{species}|{chrom}:{start}-{end}({strand})|name={gene_name}
ACGT...
```

### Step D — Attach GO labels (extension for SAE)

**Goal**
Create `labels_go.tsv` mapping each gene instance to a set of GO terms.

#### D1) If you already have Entrez GeneID (e.g., E. coli)

- Join directly:
  - `gene2go.tsv`: filter by `tax_id`, collect GO_IDs per GeneID.

Output columns example:

| gene_key | entrez_gene_id | go_terms | go_bp | go_mf | go_cc |
| -------- | -------------: | -------- | ----- | ----- | ----- |

Where `go_terms` is a `;`-joined list, and aspect-specific columns are optional.

#### D2) If you have Ensembl gene IDs (e.g., mouse)

- Map **Ensembl gene ID → Entrez GeneID** using `gene2ensembl.tsv`.
- Then join Entrez GeneID to `gene2go.tsv`.

**Mapping rules**

- `gene2ensembl.tsv` columns include: `tax_id`, `GeneID`, `Ensembl_gene_identifier`, ...
- Filter rows by `tax_id == 10090` for mouse.
- Build dict: `ensg -> GeneID`.
- For each gene: `entrez_gene_id = dict.get(gene_key)`

**Coverage expectation**

- Some Ensembl genes may not map to an Entrez GeneID; record `entrez_gene_id = ""` and `go_terms = ""`.

### Step E — (Optional) GO-slim compression

If you prefer a smaller concept space:

- Use a GO-slim mapping file (external) to convert GO_IDs → GO-slim terms.
- Keep both full GO and GO-slim columns for flexibility.

---

## 4) Outputs (What downstream agents should expect)

After running the pipeline for an organism, the output directory should contain:

- `genes.tsv` — gene coordinate manifest (one row per gene)
- `loci.fasta` — gene-locus sequences (one record per gene)
- `loci.tsv` — locus metadata (length, coordinates, id fields)
- `labels_go.tsv` — GO labels per gene (possibly empty for unmapped genes)
- `summary.json` — counts and basic QC stats

**Recommended stats in summary.json**

- num_genes_total
- num_loci_extracted
- num_missing_contig_in_fasta
- num_entrez_mapped
- num_go_labeled
- distribution of locus lengths (min/median/max)

## 5) Example Commands (expected usage)

Assuming you implement a CLI script `prepare_gene_locus_dataset.py`:

```bash
# E. coli
python prepare_gene_locus_dataset.py \
  --species ecoli_k12 \
  --tax_id 511145 \
  --genome_fasta data/ecoli/ecoli.fna \
  --annotation data/ecoli/ecoli.gff3 \
  --annotation_format gff3 \
  --gene2go data/ncbi/gene2go.tsv \
  --out_dir out/ecoli_gene_locus

# Mouse
python prepare_gene_locus_dataset.py \
  --species mouse_grcm39 \
  --tax_id 10090 \
  --genome_fasta data/mouse/mouse.fa \
  --annotation data/mouse/mouse.gtf \
  --annotation_format gtf \
  --gene2go data/ncbi/gene2go.tsv \
  --gene2ensembl data/ncbi/gene2ensembl.tsv \
  --out_dir out/mouse_gene_locus
```

---

## 6) Notes for SAE Integration (why this matters)

This dataset supports **sequence-level SAE** experiments where:

- input = frozen DNA LM embedding of gene-locus sequence
- SAE learns a sparse dictionary over embeddings
- GO labels provide **concept grounding** (post hoc mapping / evaluation)

**Important:** GO is not guaranteed to be fully determined by DNA sequence alone (especially for regulatory BP terms). For best signal, start with **GO-slim** and/or restrict to subsets (e.g., MF-heavy terms) during early prototyping.

---

## 7) Extending to other organisms (agent instructions)

To add a new organism `X`:

1) Decide data source (NCBI RefSeq vs Ensembl). Download FASTA + annotation.
2) Implement/choose the correct parser (GFF3 or GTF).
3) Determine the identifier that will serve as `gene_key`.
4) Choose GO label strategy:
   - If gene_key is Entrez GeneID → direct gene2go join
   - If gene_key is Ensembl gene ID → map using gene2ensembl → then join gene2go
   - If neither is available → add an additional mapping file (organism-specific)
5) Run acceptance tests above.
6) Commit the `summary.json` + counts so you can track coverage changes over time.
