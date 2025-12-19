#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyse des clusters CD-HIT pour produire :
- diagramme de Venn (clusters only de_novo (5-letter_prefix_TE...), only RepBase, both)
- pie chart (parmi clusters contenant au moins un de_novo: single-species vs multi-species)
- heatmap de co-occurrence des codes espèces (5-letter prefix)

Sorties : plusieurs PNG et CSV avec préfixe donné.

Usage:
  cluster_species_overlap.py -i file.clstr -o out_prefix

Le script est volontairement autonome et léger (pas de dépendances externes
non présentes dans le dépôt à l'exception de matplotlib et pandas).
"""

import argparse
import re
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def parse_clstr_minimal(path):
    """Parse un fichier .clstr et retourne une DataFrame minimale.
    Colonnes : cluster, raw_name, name, source ('de_novo' ou 'RepBase'), species_code (pour 'de_novo' si trouvée)
    """
    data = []
    cluster_id = None
    cluster_pattern = re.compile(r'^>Cluster\s+(\d+)')
    seq_pattern = re.compile(r'^\s*\d+\s+(\d+)nt,\s+>([^\.]+)\.\.\.\s*(.*)$')

    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            m = cluster_pattern.match(line)
            if m:
                cluster_id = int(m.group(1))
                continue
            m = seq_pattern.match(line)
            if not m:
                continue
            length = int(m.group(1))
            seq_info = m.group(2)

            # If seq_info contains a prefix like RM2|..., keep last part
            if '|' in seq_info:
                raw_name = seq_info.split('|', 1)[1]
            else:
                raw_name = seq_info

            # classification after # (keep raw_name before '#')
            name = raw_name.split('#', 1)[0]

            # Source detection: 'de novo' if name matches the pattern
            #   5-letter prefix followed by '_TE' and a number, e.g. AGAAC_TE507...,
            # otherwise consider it RepBase.
            m_denovo = re.match(r'^([A-Za-z]{5})_TE\d+', name)
            if m_denovo:
                src = 'de_novo'
                species_code = m_denovo.group(1)
            else:
                src = 'RepBase'

            data.append({
                'cluster': cluster_id,
                'length': length,
                'raw_name': raw_name,
                'name': name,
                'source': src,
                'species': species_code
            })

    df = pd.DataFrame(data)
    return df


def compute_cluster_sets(df: pd.DataFrame):
    """Retourne ensembles de clusters: only_de_novo, only_RepBase, both_sources"""
    clusters = df.groupby('cluster')['source'].apply(set)
    only_denovo = set(clusters[clusters.apply(lambda s: s == {'de_novo'})].index)
    only_rep = set(clusters[clusters.apply(lambda s: s == {'RepBase'})].index)
    both = set(clusters[clusters.apply(lambda s: len(s) > 1)].index)
    return only_denovo, only_rep, both


def plot_venn_clusters(only_denovo, only_rep, both, out_png):
    # venn2 expects sets A, B where intersection is provided by set operations
    denovo_all = only_denovo | both
    rep_all = only_rep | both
    plt.figure(figsize=(6, 5))
    v = venn2([denovo_all, rep_all], set_labels=('de_novo (5-letter prefix)', 'RepBase'))
    # add counts as integers
    for text in v.set_labels:
        if text:
            text.set_fontsize(12)
    plt.title('Clusters: de_novo vs RepBase', fontsize=14)
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()


def pie_single_vs_multi(df: pd.DataFrame, out_png, out_csv):
    # consider clusters that contain at least one de_novo sequence
    grp = df[df['source'] == 'de_novo'].groupby('cluster')
    single = 0
    multi = 0
    unknown = 0
    rows = []
    for cl, sub in grp:
        species = set([s for s in sub['species'] if s is not None])
        if len(species) == 0:
            unknown += 1
            rows.append((cl, 'unknown', 0))
        elif len(species) == 1:
            single += 1
            rows.append((cl, list(species)[0], 1))
        else:
            multi += 1
            rows.append((cl, ';'.join(sorted(species)), 2))

    # save CSV of clusters -> species membership
    with open(out_csv, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['cluster', 'denovo_species_codes', 'category'])
        for r in rows:
            w.writerow(r)

    labels = []
    sizes = []
    if single > 0:
        labels.append('single-species (de_novo)')
        sizes.append(single)
    if multi > 0:
        labels.append('multi-species (de_novo)')
        sizes.append(multi)
    if unknown > 0:
        labels.append('unknown-species')
        sizes.append(unknown)

    plt.figure(figsize=(6, 6))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#2ca02c', '#1f77b4', 'lightgray'])
    plt.axis('equal')
    plt.title('Clusters with at least one de_novo sequence\n(single vs multi species)')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()


def build_cooccurrence_matrix(df: pd.DataFrame):
    # consider only de_novo sequences and clusters that have at least one de_novo
    denovo = df[df['source'] == 'de_novo']
    cluster_to_species = denovo.groupby('cluster')['species'].apply(lambda s: set([x for x in s if x is not None]))
    all_species = sorted({sp for s in cluster_to_species for sp in s})
    idx = {s: i for i, s in enumerate(all_species)}

    mat = np.zeros((len(all_species), len(all_species)), dtype=int)
    for species_set in cluster_to_species:
        sset = set([x for x in species_set if x is not None])
        for a in sset:
            for b in sset:
                mat[idx[a], idx[b]] += 1

    dfm = pd.DataFrame(mat, index=all_species, columns=all_species).astype(float)
    # Set diagonal entries to NA (when row name == column name)
    idx = np.diag_indices_from(dfm.values)
    dfm.values[idx] = np.nan
    return dfm


def plot_heatmap(dfm: pd.DataFrame, out_png):
    plt.figure(figsize=(8, 7))
    im = plt.imshow(dfm.values, interpolation='nearest', cmap='viridis')
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.xticks(range(len(dfm.columns)), dfm.columns, rotation=90)
    plt.yticks(range(len(dfm.index)), dfm.index)
    plt.xlabel('Species code (5-letter prefix)')
    plt.ylabel('Species code (5-letter prefix)')
    plt.title('Co-occurrence of species within the same clusters (count of clusters)')
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Cluster-species overlap analysis (de_novo vs RepBase)')
    parser.add_argument('-i', '--input', required=True, help='Input .clstr file')
    parser.add_argument('-o', '--output', required=True, help='Output prefix (without extension)')
    args = parser.parse_args()

    df = parse_clstr_minimal(args.input)
    print(df)
    df.to_csv(f'{args.output}_parsed.tsv', sep='\t', index=False)

    only_rm, only_rep, both = compute_cluster_sets(df)

    # Venn (clusters)
    plot_venn_clusters(only_rm, only_rep, both, f'{args.output}_venn_RM_repbase.png')

    # Pie: among clusters with at least one RM
    pie_single_vs_multi(df, f'{args.output}_rm_single_vs_multi_pie.png', f'{args.output}_rm_clusters_species.tsv')

    # Heatmap co-occurrence species0
    dfm = build_cooccurrence_matrix(df)
    if dfm.shape[0] > 0:
        dfm.to_csv(f'{args.output}_cooccurrence_matrix.csv')
        plot_heatmap(dfm, f'{args.output}_species_cooccurrence_heatmap.png')
    else:
        print('No RM species codes found: heatmap not produced.')

    # Summary print
    print('Outputs:')
    print(f' - parsed TSV: {args.output}_parsed.tsv')
    print(f' - venn PNG: {args.output}_venn_RM_repbase.png')
    print(f' - pie PNG: {args.output}_rm_single_vs_multi_pie.png')
    print(f' - clusters->species TSV: {args.output}_rm_clusters_species.tsv')
    print(f' - cooccurrence CSV: {args.output}_cooccurrence_matrix.csv')
    print(f' - heatmap PNG: {args.output}_species_cooccurrence_heatmap.png')


if __name__ == '__main__':
    main()
