#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.image as mpimg
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib_venn import venn2, venn3  # venn3 pour le Venn à 3 ensembles
from upsetplot import UpSet, from_memberships


# ==========================================================
# 1) Parse the .clstr file
# ==========================================================
def parse_clstr(filename: str) -> pd.DataFrame:
    """Parse a CD-HIT .clstr file into a DataFrame, incl. RepBase as source."""
    data = []
    cluster_id = None
    cluster_pattern = re.compile(r'^>Cluster\s+(\d+)')
    seq_pattern = re.compile(r'^\s*\d+\s+(\d+)nt, >([^\.]+)\.\.\. (.*)$')

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Start of cluster
            m = cluster_pattern.match(line)
            if m:
                cluster_id = int(m.group(1))
                continue

            # Sequence line
            m = seq_pattern.match(line)
            if not m:
                continue

            length = int(m.group(1))
            seq_info = m.group(2)
            end_info = m.group(3)

            # Source and name
            src_m = re.match(r'^(RM2|HC|HLC)\|(.+)$', seq_info)
            if src_m:
                source = src_m.group(1)
                name = src_m.group(2)
            else:
                # Pas de préfixe => RepBase
                source = "RepBase"
                name = seq_info

            # Classification (si présent après un #)
            classification = None
            if "#" in name:
                name, classification = name.split("#", 1)

            # Representative / similarity
            if "*" in end_info:
                representative = True
                similarity = None
            else:
                sim_m = re.search(r'/([\d\.]+)%', end_info)
                similarity = float(sim_m.group(1)) if sim_m else None
                representative = False

            data.append({
                "cluster": cluster_id,
                "length": length,
                "source": source,
                "name": name,
                "classification": classification,
                "similarity": similarity,
                "representative": representative
            })

    return pd.DataFrame(data)


# ==========================================================
# 2) Cross-tool matches
# ==========================================================
def add_cross_tool_matches(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add boolean column telling if a sequence's cluster mixes tools.
    - Pour RM2 : True si le cluster contient HC, HLC ou RepBase
    - Pour HLC : True si le cluster contient RM2 ou RepBase
    - Pour HC  : conservé de toute façon, valeur None
    - Pour RepBase : retiré au filtrage, valeur None
    """
    cluster_sources = df.groupby("cluster")["source"].apply(set).to_dict()
    def has_other(row):
        srcset = cluster_sources[row["cluster"]]
        if row["source"] == "RM2":
            return any(x in srcset for x in ["HC", "HLC", "RepBase"])
        if row["source"] == "HLC":
            return any(x in srcset for x in ["RM2", "RepBase"])
        # HC conservé de toute façon ; RepBase retiré ensuite
        return None

    df = df.copy()
    df["has_other_tool_match"] = df.apply(has_other, axis=1)
    return df


# ==========================================================
# 3) Filter sequences
# ==========================================================
def filter_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep RM2/HLC only if cross-tool match; always keep HC.
    Remove all RepBase sequences (RepBase is only used for overlap/confirmation).
    """
    mask = (
        ((df["source"].isin(["RM2", "HLC"])) & (df["has_other_tool_match"] == True))
        | (df["source"] == "HC")
    )
    # Retirer RepBase explicitement
    mask &= (df["source"] != "RepBase")
    return df.loc[mask].reset_index(drop=True)


# ==========================================================
# 4) UpSet (no title in the file)
# ==========================================================
def plot_upset_from_clusters(df: pd.DataFrame, outfile: str):
    """Save an UpSet plot showing cluster overlaps (no title)."""
    rm2_classified = set(df[(df["source"] == "RM2") & (df["classification"] != "Unknown")]["cluster"])
    rm2_unclassified = set(df[(df["source"] == "RM2") & (df["classification"] == "Unknown")]["cluster"])
    hite_confident = set(df[df["source"] == "HC"]["cluster"])
    hite_low_confident = set(df[df["source"] == "HLC"]["cluster"])

    all_clusters = rm2_classified | rm2_unclassified | hite_confident | hite_low_confident

    memberships = [
        [
            grp for grp, cluster_set in {
                "RM2_classified": rm2_classified,
                "RM2_unclassified": rm2_unclassified,
                "HiTE_confident": hite_confident,
                "HiTE_low_confident": hite_low_confident
            }.items() if cl in cluster_set
        ]
        for cl in all_clusters
    ]

    upset_data = from_memberships(memberships)
    plt.figure(figsize=(8, 5))
    UpSet(upset_data, subset_size='count', show_counts=True).plot()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()


# ==========================================================
# 5) Venn RM2 vs HiTE vs RepBase (no title in the file)
# ==========================================================
def plot_venn_rm2_vs_hite(df: pd.DataFrame, outfile: str):
    """Save a 3-set Venn (RM2, HiTE, RepBase; no title)."""
    # NB: on lit dans le DF PARSÉ (non filtré) si on veut inclure RepBase
    rm2_clusters = set(df[df["source"] == "RM2"]["cluster"])
    hite_clusters = set(df[df["source"].isin(["HC", "HLC"])]["cluster"])
    repbase_clusters = set(df[df["source"] == "RepBase"]["cluster"])

    plt.figure(figsize=(7.2, 6.2))
    venn3(
        [rm2_clusters, hite_clusters, repbase_clusters],
        set_labels=("RepeatModeler2 (RM2)", "HiTE (HC + HLC)", "RepBase"),
        # couleurs douces, alpha comme avant
        set_colors=("skyblue", "lightgreen", "orange"),
        alpha=0.7
    )
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()


# ==========================================================
# 6) Cluster size distribution (bigger & clearer)
# ==========================================================
def plot_cluster_size_distribution(df: pd.DataFrame, outfile: str):
    """
    Plot cluster size distribution with larger text and bars.
    (No title in the file to avoid duplicates on the A4 page.)
    """
    cluster_sizes = df.groupby("cluster").size()

    plt.figure(figsize=(9.5, 6.5))
    plt.hist(
        cluster_sizes,
        bins=range(1, cluster_sizes.max() + 2),
        edgecolor='black',
        alpha=0.85,
        rwidth=0.92
    )

    plt.xlabel("Number of sequences per cluster", fontsize=24, weight="bold")
    plt.ylabel("Number of clusters", fontsize=24, weight="bold")

    plt.yscale("log")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.tick_params(axis='both', which='major', labelsize=22)
    plt.grid(alpha=0.45, which="both", linestyle="--", linewidth=0.6)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()


# ==========================================================
# 7) Confirmation pies — 2×2, with connector lines (no donut)
# ==========================================================
def plot_confirmation_pies(df: pd.DataFrame, outfile_prefix: str = "confirmation_pies"):
    """
    Build four pies on a 2x2 grid with bigger text.
    Percent labels are placed outside with connector lines (no donut).
    No global title; captions will be added on the A4 page.

    Ajoute 'confirmed by RepBase' si une autre séquence du même cluster a source RepBase.
    Priorité: RM2 classified > RM2 unclassified > HC > HLC > RepBase > without confirmation.
    """
    cluster_groups = df.groupby("cluster")

    def get_confirmation(row):
        cl = row["cluster"]
        others = cluster_groups.get_group(cl)
        others = others[others["name"] != row["name"]]
        if len(others) == 0:
            return "without confirmation"
        if any((others["source"] == "RM2") & (others["classification"] != "Unknown")):
            return "confirmed by RM2 classified"
        if any((others["source"] == "RM2") & (others["classification"] == "Unknown")):
            return "confirmed by RM2 unclassified"
        if any(others["source"] == "HC"):
            return "confirmed by HC"
        if any(others["source"] == "HLC"):
            return "confirmed by HLC"
        if any(others["source"] == "RepBase"):
            return "confirmed by RepBase"
        return "without confirmation"

    df = df.copy()
    df["confirmation"] = df.apply(get_confirmation, axis=1)

    subsets = {
        "RepeatModeler classified": df[(df["source"] == "RM2") & (df["classification"] != "Unknown")],
        "RepeatModeler unclassified": df[(df["source"] == "RM2") & (df["classification"] == "Unknown")],
        "HiTE confident": df[df["source"] == "HC"],
        "HiTE low confident": df[df["source"] == "HLC"]
    }

    categories = [
        "without confirmation",
        "confirmed by RM2 classified",
        "confirmed by RM2 unclassified",
        "confirmed by HC",
        "confirmed by HLC",
        "confirmed by RepBase"
    ]
    colors = {
        "without confirmation": "lightgray",
        "confirmed by RM2 classified": "#1f77b4",
        "confirmed by RM2 unclassified": "#aec7e8",
        "confirmed by HC": "#2ca02c",
        "confirmed by HLC": "#98df8a",
        "confirmed by RepBase": "orange"
    }

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
    axes = axes.flatten()

    for ax, (panel_title, subset) in zip(axes, subsets.items()):
        counts = subset["confirmation"].value_counts().reindex(categories, fill_value=0)
        total = int(counts.sum())

        # Classic full pie (no autopct)
        wedges, texts = ax.pie(
            counts,
            labels=None,
            colors=[colors[c] for c in counts.index],
            startangle=90,
            counterclock=False
        )

        # Add percentage labels outside with connectors
        for idx, w in enumerate(wedges):
            if total == 0:
                continue
            pct = 100.0 * counts.iloc[idx] / total
            if pct < 1.0:
                continue  # skip tiny slices

            # Mid-angle of the wedge
            ang = (w.theta2 + w.theta1) / 2.0
            x = np.cos(np.deg2rad(ang))
            y = np.sin(np.deg2rad(ang))

            # Outside label position
            text_x = 1.22 * x
            text_y = 1.22 * y

            ax.annotate(
                f"{pct:.1f}%",
                xy=(x, y),
                xytext=(text_x, text_y),
                ha="left" if x >= 0 else "right",
                va="center",
                fontsize=14,
                fontweight="bold",
                arrowprops=dict(
                    arrowstyle="-",
                    color="black",
                    connectionstyle=f"arc3,rad={0.15*np.sign(x)}",
                    shrinkA=0, shrinkB=0, lw=1.0
                )
            )

        # Caption under each pie
        ax.text(
            0.5, -0.18, f"{panel_title} (n={total})",
            ha="center", va="top", transform=ax.transAxes,
            fontsize=18, weight="bold"
        )

    # Shared legend, placed to the right
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', label=cat,
                   markerfacecolor=colors[cat], markersize=10)
        for cat in categories
    ]
    fig.legend(
        handles=handles,
        loc='center left',
        bbox_to_anchor=(1.01, 0.5),
        fontsize=18,
        title="Confirmation type",
        title_fontsize=22
    )

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # place pour la légende
    plt.savefig(f"{outfile_prefix}_confirmation_pies.png", dpi=300, bbox_inches='tight')
    plt.close()


# ==========================================================
# 8) A4 summary page (titles live here)
# ==========================================================
def make_a4_summary_page(species_name: str, out_prefix: str,
                         outfile_png: str = None, outfile_pdf: str = None):
    if outfile_png is None:
        outfile_png = f"{out_prefix}_summary_A4.png"
    if outfile_pdf is None:
        outfile_pdf = f"{out_prefix}_summary_A4.pdf"

    files = {
        "UpSet": f"{out_prefix}_upset_TE_clusters.png",
        "Venn": f"{out_prefix}_venn_TE_clusters.png",
        "Size": f"{out_prefix}_cluster_size_distribution.png",
        "Pies": f"{out_prefix}_confirmation_pies.png",
    }

    fig = plt.figure(figsize=(8.27, 11.69), constrained_layout=True)
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.05, 1.0, 1.35])

    def put_image(ax, path, caption: str):
        try:
            img = mpimg.imread(path)
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(caption, fontsize=13, weight="bold")
        except FileNotFoundError:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{caption}\n(file not found)\n{path}",
                    ha="center", va="center", fontsize=10.5, color="red")

    put_image(fig.add_subplot(gs[0, :]), files["UpSet"], "UpSet — cluster overlap (4 groups)")
    put_image(fig.add_subplot(gs[1, 0]), files["Venn"], "Venn — RM2 vs HiTE vs RepBase")
    put_image(fig.add_subplot(gs[1, 1]), files["Size"], "Cluster size distribution")
    put_image(fig.add_subplot(gs[2, :]), files["Pies"], "Sequence confirmation per source")

    fig.suptitle(f"{species_name} — CD-HIT clustering summary", fontsize=18, weight="bold")

    fig.savefig(outfile_png, dpi=300, bbox_inches="tight")
    with PdfPages(outfile_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ==========================================================
# 9) Main
# ==========================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse and visualize CD-HIT clustering results.")
    parser.add_argument("-i", "--input", required=True, help="Input .clstr file path")
    parser.add_argument("-n", "--name", required=True, help="Species name (used on the A4 summary)")
    parser.add_argument("-o", "--output", required=True, help="Output file prefix (without extension)")
    args = parser.parse_args()

    infile = args.input
    species_name = args.name
    out_prefix = args.output

    # Parse & process
    df = parse_clstr(infile)
    df = add_cross_tool_matches(df)
    df_filtered = filter_sequences(df)

    # Tables
    df.to_csv(f"{out_prefix}_parsed.tsv", sep="\t", index=False)
    df_filtered.to_csv(f"{out_prefix}_filtered.tsv", sep="\t", index=False)

    # Individual figures (WITHOUT titles)
    # UpSet et distribution utilisent df_filtered (comme avant)
    plot_upset_from_clusters(df_filtered, f"{out_prefix}_upset_TE_clusters.png")
    plot_cluster_size_distribution(df_filtered, f"{out_prefix}_cluster_size_distribution.png")

    # Venn 3-ensembles : utilise df PARSÉ pour inclure RepBase
    plot_venn_rm2_vs_hite(df, f"{out_prefix}_venn_TE_clusters.png")

    # Pies (sur df_filtered) avec confirmation RepBase détectée dans le cluster (même si RepBase supprimée)
    # NB: la confirmation RepBase se base sur la présence d'une séquence RepBase dans df (parse),
    # mais ici on calcule la confirmation à partir de df_filtered ; si tu préfères que la
    # confirmation « voie » aussi RepBase, passe df complet ci-dessous (au prix d'inclure RepBase
    # dans les subsets — non désiré). La fonction actuelle cherche RepBase dans "others" du cluster,
    # donc il faut qu'elles existent dans df donné. On utilise donc df complet ici.
    plot_confirmation_pies(df, out_prefix)

    # A4 summary (with panel captions and global title)
    make_a4_summary_page(species_name, out_prefix)

    print("✅ Analysis complete. Results saved:")
    print(f" - Parsed table: {out_prefix}_parsed.tsv")
    print(f" - Filtered table: {out_prefix}_filtered.tsv")
    print(f" - UpSet plot: {out_prefix}_upset_TE_clusters.png")
    print(f" - Venn diagram: {out_prefix}_venn_TE_clusters.png")
    print(f" - Cluster size distribution: {out_prefix}_cluster_size_distribution.png")
    print(f" - Confirmation pies: {out_prefix}_confirmation_pies.png")
    print(f" - A4 summary (PNG): {out_prefix}_summary_A4.png")
    print(f" - A4 summary (PDF): {out_prefix}_summary_A4.pdf")
