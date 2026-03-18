#!/usr/bin/env python
"""
Plot BUSCO Q-vecs (Single / Duplicated / Fragmented / Missing) for all BUSCO JSONs
found directly in an input directory (non-recursive).

Usage:
    python bin/plot_busco_qvecs.py -i input_dir -o out.png

The script is permissive when parsing BUSCO JSONs: it looks for common keys
such as `counts`, `total_buscos`, `percentages`, or explicit keys like
`complete_single_copy`, `complete_duplicated`, `fragmented`, `missing`.
"""

import os
import re
import json
import argparse
from collections import defaultdict
import pandas as pd
from matplotlib import pyplot as plt


def build_arg_parser():
    p = argparse.ArgumentParser(description='Plot BUSCO qvecs for all JSONs in a folder')
    p.add_argument('-i', '--input', required=True, help='Input folder containing BUSCO .json files')
    p.add_argument('-o', '--output', default='busco_qvecs.png', help='Output PNG file')
    p.add_argument('--width', type=float, default=None, help='Figure width in inches')
    p.add_argument('--height', type=float, default=None, help='Figure height in inches')
    p.add_argument('--fontsize', type=int, default=10, help='Font size for labels/legend')
    p.add_argument('--no-labels', action='store_true', help='Hide x tick labels')
    return p


def parse_busco_json(path):
    """Parse a BUSCO JSON file and return a dict with percentages for
    Single, Duplicated, Fragmented and Missing, and total_buscos when available.
    Returns None if parsing fails or required fields missing.
    """
    with open(path) as fh:
        try:
            data = json.load(fh)
        except Exception:
            return None

    # First, prefer explicit values in `results` (common in BUSCO JSONs)
    results = data.get('results', {}) if isinstance(data.get('results'), dict) else {}

    def pick(*names):
        for n in names:
            if n in results:
                return results[n]
        return None

    # Common BUSCO keys (as in the example):
    # "Single copy percentage", "Multi copy percentage", "Fragmented percentage", "Missing percentage", "n_markers" or "n"
    single_p = pick('Single copy percentage', 'Single copy %', 'Single copy percentage (%)', 'Single copy')
    dup_p = pick('Multi copy percentage', 'Multi copy %', 'Multi copy percentage (%)', 'Multi copy', 'Multi copy BUSCOs')
    frag_p = pick('Fragmented percentage', 'Fragmented %', 'Fragmented percentage (%)', 'Fragmented')
    miss_p = pick('Missing percentage', 'Missing %', 'Missing percentage (%)', 'Missing')
    total = pick('n_markers', 'n', 'n_markers') or data.get('lineage_dataset', {}).get('number_of_buscos') or results.get('n_markers')
    one_line = results.get('one_line_summary') or data.get('one_line_summary')

    # If the fields exist in `results`, convert to floats and return
    try:
        if single_p is not None or dup_p is not None or frag_p is not None or miss_p is not None:
            sp = float(single_p) if single_p is not None else None
            dp = float(dup_p) if dup_p is not None else None
            fp = float(frag_p) if frag_p is not None else 0.0
            mp = float(miss_p) if miss_p is not None else 0.0
            # If single provided and dup missing but Complete provided, try to compute
            if sp is None and dp is not None and results.get('Complete percentage') is not None:
                sp = float(results.get('Complete percentage')) - float(dp)
            # Fill missing with 0.0
            sp = sp if sp is not None else 0.0
            dp = dp if dp is not None else 0.0
            total_int = int(total) if total not in (None, '') else None
            return {'Single': sp, 'Duplicated': dp, 'Fragmented': fp, 'Missing': mp, 'Total': total_int, 'one_line': one_line}
    except Exception:
        pass

    # Fallback: reuse previous more generic parsing behaviour (counts/percentages elsewhere)
    # Attempt older logic: look for counts or nested percentages
    # helpers to get keys flexibly
    def get(d, *keys):
        for k in keys:
            if k in d:
                return d[k]
        return None

    counts = get(data, 'counts') or get(results, 'counts')
    total = get(data, 'total_buscos') or get(data, 'total') or get(data, 'n_total_buscos') or total

    if isinstance(counts, dict):
        c_single = (counts.get('complete_single_copy') or counts.get('single_copy') or counts.get('single') or counts.get('complete_single'))
        c_dup = (counts.get('complete_duplicated') or counts.get('duplicated') or counts.get('complete_duplicate'))
        c_frag = (counts.get('fragmented') or counts.get('fragged') or counts.get('fragments'))
        c_missing = (counts.get('missing') or counts.get('miss'))
        if not total:
            vals = [v for v in (c_single, c_dup, c_frag, c_missing) if isinstance(v, (int, float))]
            if vals:
                total = sum(vals)
        if total and total > 0:
            try:
                single_p = (float(c_single) / float(total)) * 100 if c_single is not None else 0.0
                dup_p = (float(c_dup) / float(total)) * 100 if c_dup is not None else 0.0
                frag_p = (float(c_frag) / float(total)) * 100 if c_frag is not None else 0.0
                miss_p = (float(c_missing) / float(total)) * 100 if c_missing is not None else 0.0
            except Exception:
                return None
            return {'Single': single_p, 'Duplicated': dup_p, 'Fragmented': frag_p, 'Missing': miss_p, 'Total': int(total), 'one_line': one_line}

    # Generic percentages dict fallback
    percentages = get(data, 'percentages') or get(results, 'percentages')
    if isinstance(percentages, dict):
        comp = (percentages.get('complete') or percentages.get('C'))
        single = (percentages.get('single_copy') or percentages.get('S') or percentages.get('single'))
        dup = (percentages.get('duplicated') or percentages.get('D') or percentages.get('duplicate'))
        frag = (percentages.get('fragmented') or percentages.get('F') or percentages.get('fragment'))
        miss = (percentages.get('missing') or percentages.get('M'))
        if single is not None and dup is not None and frag is not None and miss is not None:
            return {'Single': float(single), 'Duplicated': float(dup), 'Fragmented': float(frag), 'Missing': float(miss), 'Total': total, 'one_line': one_line}

    return None


def clean_label(name):
    """Clean a BUSCO-derived filename for display on the plot.
    Removes common prefixes such as anything before 'busco_' and leading
    'short_summary' fragments.
    """
    # remove common 'short_summary...' prefixes
    name = re.sub(r'^short[_\.-]?summary\.[^.]*\.', '', name)
    name = re.sub(r'^short[_\.-]?summary[_\.-]?', '', name)
    # remove everything up to and including 'busco_'
    name = re.sub(r'.*busco_', '', name)
    # also remove common leading tokens like 'short_' or 'short-'
    name = re.sub(r'^short[_\.-]?', '', name)
    return name


def collect_busco_jsons(input_dir):
    """Return a list of JSON file paths directly inside `input_dir` (non-recursive)."""
    if not os.path.isdir(input_dir):
        raise ValueError(f"Input path is not a directory: {input_dir}")
    files = [os.path.join(input_dir, x) for x in os.listdir(input_dir) if x.lower().endswith('.json') and os.path.isfile(os.path.join(input_dir, x))]
    return sorted(files)


def plot_busco_df(df, outpath, width=None, height=None, fontsize=10, no_labels=False):
    if width is None:
        width = max(len(df) * 0.3, 6)
    if height is None:
        # Make plot height proportional to number of samples so bars are readable
        # Use at least 8 inches and ~0.7 inch per sample
        height = max(len(df) * 0.3, 8)
    fig, ax = plt.subplots(1, 1, figsize=(width, height))

    # ensure numeric and fill NaNs
    for col in ['Single', 'Duplicated', 'Fragmented', 'Missing']:
        if col not in df.columns:
            df[col] = 0.0
    df = df.fillna(0.0)

    # plotting order and colors (BUSCO-like)
    cols = ['Single', 'Duplicated', 'Fragmented', 'Missing']
    colors = ['#5DADE2', '#2E86C1', '#F7DC6F', '#F8766D']

    y = list(range(len(df)))
    # reverse so first entry is top
    y_rev = y[::-1]

    # draw stacked horizontal bars manually to allow annotations
    for i, idx in enumerate(df.index[::-1]):
        left = 0.0
        row = df.loc[idx]
        total_width = 0.0
        for col, color in zip(cols, colors):
            width_val = float(row[col]) if col in row else 0.0
            total_width += width_val
        for col, color in zip(cols, colors):
            width_val = float(row[col]) if col in row else 0.0
            if width_val <= 0:
                continue
            ax.barh(i, width_val, left=left, color=color, edgecolor='white', height=0.7)
            # annotate percentage inside segment if it is wide enough
            if width_val >= 3.0:
                # choose text color for contrast
                text_color = 'black' if color in ['#F7DC6F', '#F8766D'] else 'white' if color in ['#2E86C1'] else 'black'
                ax.text(left + width_val / 2.0, i, f"{width_val:.1f}%", va='center', ha='center', fontsize=max(fontsize-1, 6), color=text_color)
            left += width_val

        # (removed BUSCO-style one-line text to keep plot clean)

    ax.set_yticks(list(range(len(df))))
    ax.set_yticklabels(df.index[::-1], fontsize=fontsize)
    ax.set_xlim(0, 100)
    ax.set_xlabel('% BUSCOs', fontsize=fontsize)
    ax.invert_yaxis()
    if no_labels:
        ax.get_yaxis().set_visible(False)

    # Legend
    legend_patches = [plt.Rectangle((0, 0), 1, 1, facecolor=c) for c in colors]
    ax.legend(legend_patches, ['Complete (S)', 'Complete (D)', 'Fragmented (F)', 'Missing (M)'], loc='lower center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=fontsize)

    ax.set_xlim(0, 100)
    ax.grid(axis='x', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches='tight')


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    input_dir = args.input
    out = args.output

    files = collect_busco_jsons(input_dir)
    if not files:
        print(f'No .json BUSCO files found in {input_dir}')
        return

    records = []
    names = []
    for j in files:
        parsed = parse_busco_json(j)
        name = os.path.splitext(os.path.basename(j))[0]
        name = clean_label(name)
        if parsed is None:
            print(f'Warning: could not parse {j}, skipping')
            continue
        records.append(parsed)
        names.append(name)

    if not records:
        print('No parsable BUSCO JSONs found.')
        return

    df = pd.DataFrame(records, index=names)
    # sort alphabetically for plotting (we sort descending so the plotting routine
    # that reverses the index produces top-to-bottom alphabetical order)
    df = df.sort_index(ascending=False)
    # ensure columns exist
    for col in ['Single', 'Duplicated', 'Fragmented', 'Missing']:
        if col not in df.columns:
            df[col] = 0.0

    plot_busco_df(df, out, width=args.width, height=args.height, fontsize=args.fontsize, no_labels=args.no_labels)
    print(f'Wrote plot to {out}')


if __name__ == '__main__':
    main()
