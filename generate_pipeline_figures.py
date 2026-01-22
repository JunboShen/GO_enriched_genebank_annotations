#!/usr/bin/env python3
"""Generate pipeline figures as PNGs for the annotation preprocessing report."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def add_box(ax, xy, text, width=2.6, height=0.6, fc="#E3F2FD", ec="#1E88E5", text_color="#0D47A1", fontsize=9):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.2,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=text_color,
        wrap=True,
    )
    return patch


def add_arrow(ax, start, end):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=1.2,
        color="#546E7A",
    )
    ax.add_patch(arrow)

def add_group_box(ax, xy, title, items, width=3.4, height=3.2, fc="#EEF3F8", ec="#90A4AE"):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.04,rounding_size=0.12",
        linewidth=1.2,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.18,
        y + height - 0.35,
        title,
        ha="left",
        va="center",
        fontsize=10,
        color="#263238",
        weight="bold",
    )
    text = "\n".join([f"• {item}" for item in items])
    ax.text(
        x + 0.18,
        y + height - 0.8,
        text,
        ha="left",
        va="top",
        fontsize=8.6,
        color="#37474F",
    )
    return patch

def fig_data_lineage(out_path: Path):
    fig, ax = plt.subplots(figsize=(12.5, 4.6))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")
    ax.set_axis_off()
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 5)

    g1 = add_group_box(
        ax,
        (0.4, 0.9),
        "Inputs",
        [
            "GenBank: NC_000913.gb",
            "NCBI gene2go (gene2go.gz)",
            "GO ontology: go-basic.obo",
            "GO-slim subset: goslim_generic.obo",
            "Normalization rules (JSON)",
        ],
        width=3.2,
        height=3.2,
        fc="#E3F2FD",
        ec="#64B5F6",
    )
    g2 = add_group_box(
        ax,
        (4.3, 0.9),
        "Preprocessing / Mapping",
        [
            "filter_gene2go_by_taxon.py",
            "link_genbank_to_go_gene2go.py",
            "make_goslim_map.py",
        ],
        width=3.8,
        height=3.2,
        fc="#E8F5E9",
        ec="#81C784",
    )
    g3 = add_group_box(
        ax,
        (8.7, 0.9),
        "Table Construction",
        [
            "build_ecoli_go_enriched_table.py",
            "feature extraction + QC",
            "attach GO/GO-slim + labels",
        ],
        width=3.2,
        height=3.2,
        fc="#FFF3E0",
        ec="#FFB74D",
    )
    g4 = add_group_box(
        ax,
        (12.4, 0.9),
        "Outputs",
        [
            "ecoli_go_enriched_table.csv",
            "goslim_map_ecoli.tsv",
            "genbank_go_annotation…csv",
            "stats/*.csv",
        ],
        width=3.0,
        height=3.2,
        fc="#F3E5F5",
        ec="#BA68C8",
    )

    add_arrow(ax, (3.6, 2.5), (4.3, 2.5))
    add_arrow(ax, (8.1, 2.5), (8.7, 2.5))
    add_arrow(ax, (11.9, 2.5), (12.4, 2.5))

    fig.suptitle("End-to-end data lineage (files + scripts)", fontsize=12, color="#263238")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=240)
    plt.close(fig)


def fig_label_hierarchy(out_path: Path):
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")
    ax.set_axis_off()
    ax.set_xlim(0, 6.4)
    ax.set_ylim(0, 6.4)

    y_positions = [4.9, 3.5, 2.1, 0.7]
    labels = [
        "L0: Feature type\n(CDS, tRNA, mobile_element, …)",
        "L1: Normalized product concept\n(transporter, kinase, efflux_pump, …)",
        "L2: GO terms\n(4,161 unique in E. coli test)",
        "L3: GO-slim terms\n(86 unique in E. coli test)",
    ]
    colors = [
        ("#E3F2FD", "#64B5F6"),
        ("#E8F5E9", "#81C784"),
        ("#FFF3E0", "#FFB74D"),
        ("#F3E5F5", "#BA68C8"),
    ]
    for y, text, (fc, ec) in zip(y_positions, labels, colors):
        add_box(ax, (0.7, y), text, width=5.0, height=0.78, fc=fc, ec=ec, text_color="#263238", fontsize=9.5)

    for i in range(len(y_positions) - 1):
        add_arrow(ax, (3.2, y_positions[i]), (3.2, y_positions[i] - 0.4))

    fig.suptitle("Label hierarchy used for SAE concept mapping", fontsize=12, color="#263238")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=240)
    plt.close(fig)


def fig_goslim_mapping(out_path: Path):
    fig, ax = plt.subplots(figsize=(10.2, 3.4))
    fig.patch.set_facecolor("#FAFAFA")
    ax.set_facecolor("#FAFAFA")
    ax.set_axis_off()
    ax.set_xlim(0, 12.5)
    ax.set_ylim(0, 3.6)

    add_box(ax, (0.4, 1.2), "Annotated GO term\n(e.g., GO:0022857)", width=2.6, height=0.9,
            fc="#E3F2FD", ec="#64B5F6", text_color="#263238")
    add_box(ax, (3.2, 1.2), "Traverse parents\n(is_a, part_of)\nusing go-basic.obo", width=2.8, height=0.9,
            fc="#E8F5E9", ec="#81C784", text_color="#263238")
    add_box(ax, (6.2, 1.2), "Intersect with GO-slim set\n(goslim_generic.obo)", width=2.8, height=0.9,
            fc="#FFF3E0", ec="#FFB74D", text_color="#263238")
    add_box(ax, (9.4, 1.2), "GO-slim term(s)\n(e.g., GO:0005215)", width=2.6, height=0.9,
            fc="#F3E5F5", ec="#BA68C8", text_color="#263238")

    add_arrow(ax, (3.0, 1.65), (3.2, 1.65))
    add_arrow(ax, (6.0, 1.65), (6.2, 1.65))
    add_arrow(ax, (8.9, 1.65), (9.4, 1.65))

    fig.suptitle("GO → GO-slim mapping intuition", fontsize=12, color="#263238")
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(out_path, dpi=240)
    plt.close(fig)


def main() -> int:
    fig_data_lineage(OUT_DIR / "fig1_data_lineage.png")
    fig_label_hierarchy(OUT_DIR / "fig2_label_hierarchy.png")
    fig_goslim_mapping(OUT_DIR / "fig3_goslim_mapping.png")
    print(f"Wrote figures to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
