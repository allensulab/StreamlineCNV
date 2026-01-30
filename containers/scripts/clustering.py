#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNV clustering on a GLOBAL 500kb grid built from chromosome sizes.

Usage:
  python clustering_global_grid_sizes.py -c chrom.sizes -i seg_anno -o cnv_heatmap.pdf --bin-size 500000 --neutral 3 --anchor 0 --coords one_based_inclusive

Inputs:
  chrom.sizes (TSV/CSV): columns: chr, size
  seg_anno (TSV):        columns: sample, chr, start, end, state, label

Key behavior:
  • Build a global window list per chromosome using chr sizes and bin size.
  • Map each segment to all overlapping windows via index math (no per-base loops).
  • Pivot to sample x window matrix; fill missing with neutral (default 3).
  • Cluster rows (samples); show CNV states with colored rectangle patches.
"""

import argparse
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from sklearn.preprocessing import LabelEncoder
from typing import List


def parse_args():
    p = argparse.ArgumentParser(description="Global-grid CNV clustering heatmap (from chrom sizes)")
    p.add_argument("-c", "--chrom-sizes", required=True,
                   help="Chromosome sizes file (TSV/CSV) with columns: chr,size")
    p.add_argument("-i", "--input", default="seg_anno",
                   help="Segment TSV with columns: sample, chr, start, end, state, label (default: seg_anno)")
    p.add_argument("-o", "--output", required=True, help="Output PDF filename")
    p.add_argument("--bin-size", type=int, default=500000,
                   help="Bin size in bp (default: 500000)")
    p.add_argument("--neutral", type=float, default=3,
                   help="Neutral state to fill missing bins (default: 3)")
    p.add_argument("--drop-chr", nargs="*", default=None,
                   help="Chromosomes to drop (e.g., --drop-chr Y)")
    p.add_argument("--anchor", type=int, default=0, choices=[0, 1],
                   help="Grid anchor: 0 or 1 (default 0).")
    p.add_argument("--coords", choices=["one_based_inclusive", "zero_based_half_open"],
                   default="one_based_inclusive",
                   help="Coordinate convention for input segments (default: one_based_inclusive)")
    p.add_argument("--figsize", type=float, nargs=2, default=(18, 10),
                   help="Figure size width height (default: 18 10)")
    p.add_argument("--palette", default="Set3",
                   help="Seaborn palette for labels (default: Set3)")
    p.add_argument("--no-sampleLabel", type=lambda x: str(x).lower() == 'true', default=False,
               help="If True, disable Y-axis sample labeling (default: False)")
    return p.parse_args()


def read_chrom_sizes(path: str) -> pd.DataFrame:
    # Try TSV then CSV
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.read_csv(path)
    # normalize columns
    cols = {c.lower(): c for c in df.columns}
    if "chr" not in cols and "chromosome" in cols:
        df.rename(columns={cols["chromosome"]: "chr"}, inplace=True)
        cols = {c.lower(): c for c in df.columns}
    if "chr" not in cols or "size" not in cols:
        raise ValueError("chrom sizes file must have columns 'chr' and 'size'")
    df = df[["chr", "size"]].copy()
    df["chr"] = df["chr"].astype(str)
    df["size"] = df["size"].astype(int)
    return df


def build_global_windows_zero_based(chrom_df, bin_size, drop=None):
    rows = []
    chroms = chrom_df.copy()
    if drop:
        chroms = chroms[~chroms["chr"].astype(str).isin(map(str, drop))]
    for _, r in chroms.iterrows():
        chr_, L = str(r["chr"]), int(r["size"])
        start = 0
        while start < L:
            end = min(start + bin_size, L)  # half-open [start, end)
            rows.append({
                "chr": chr_,
                "start": start,
                "end": end,
                "segment": f"{chr_}:{start}-{end}"
            })
            start = end
    return pd.DataFrame(rows)


def segments_to_windows(df_seg: pd.DataFrame, windows: pd.DataFrame,
                        bin_size: int, anchor: int, coords: str) -> pd.DataFrame:
    # Efficient mapping: compute window indices using integer division
    out_rows = []
    # Pre-build per-chrom windows start list for bounds checking
    win_dict = {}
    for chr_, sub in windows.groupby("chr"):
        sub = sub.reset_index(drop=True)
        win_dict[chr_] = {
            "df": sub,
            "start0": sub.loc[0, "start"],
            "n": len(sub)
        }

    for _, row in df_seg.iterrows():
        chr_ = str(row["chr"])
        if chr_ not in win_dict:
            continue
        s = int(row["start"])
        e = int(row["end"])
        # Normalize to 0-based half-open coordinates matching windows definition
        if coords == "one_based_inclusive":
            s0 = s - 1
            e0 = e
        else:  # zero_based_half_open
            s0 = s
            e0 = e
        # Map to window indices
        start0 = win_dict[chr_]["start0"]
        n = win_dict[chr_]["n"]
        # Convert to grid-relative positions
        k_start = (s0 - anchor) // bin_size
        k_end = (max(e0 - 1, s0) - anchor) // bin_size
        k_start = max(k_start, 0)
        k_end = min(k_end, n - 1)
        # Emit overlaps
        for k in range(k_start, k_end + 1):
            win = win_dict[chr_]["df"].iloc[k]
            new_row = row.copy()
            new_row["start"] = int(win["start"])
            new_row["end"] = int(win["end"])
            new_row["segment"] = str(win["segment"])
            out_rows.append(new_row)
    if not out_rows:
        return pd.DataFrame(columns=list(df_seg.columns) + ["segment"])
    return pd.DataFrame(out_rows)


def main():
    args = parse_args()

    # -------------------------
    # Drop list as strings
    # -------------------------
    drop_list = [str(c) for c in args.drop_chr] if args.drop_chr else []

    # -------------------------
    # 1) Read chromosome sizes FIRST and define canonical order
    # -------------------------
    try:
        chrom_df = read_chrom_sizes(args.chrom_sizes)
    except Exception as e:
        print(f"Failed to read chrom sizes: {e}", file=sys.stderr)
        sys.exit(1)

    chrom_df["chr"] = chrom_df["chr"].astype(str)

    # Drop 'genome' pseudo-chromosome (this was creating the big white block)
    chrom_df = chrom_df[chrom_df["chr"].str.lower() != "genome"].copy()

    # Apply user-specified drops to chrom sizes
    if drop_list:
        chrom_df = chrom_df[~chrom_df["chr"].isin(drop_list)].copy()

    if chrom_df.empty:
        print("No chromosomes left in chrom sizes after filtering.", file=sys.stderr)
        sys.exit(1)

    # 👉 THIS is the chromosome order, exactly as in the chromosome size file (minus genome / drops)
    chrom_order = chrom_df["chr"].tolist()

    # -------------------------
    # 2) Read seg_anno (segments)
    # -------------------------
    try:
        df = pd.read_csv(args.input, sep="\t", dtype={"chr": str})
    except Exception as e:
        print(f"Failed to read {args.input}: {e}", file=sys.stderr)
        sys.exit(1)

    required_cols = {"sample", "chr", "start", "end", "state", "label"}
    missing = required_cols - set(df.columns)
    if missing:
        print(f"Input is missing required columns: {sorted(missing)}", file=sys.stderr)
        sys.exit(1)

    df["chr"] = df["chr"].astype(str)

    # Apply same drop_chr to seg_anno
    if drop_list:
        df = df[~df["chr"].isin(drop_list)].copy()

    # Keep only chromosomes that exist in the chromosome size file (and keep their order from chrom_order)
    df = df[df["chr"].isin(chrom_order)].copy()
    if df.empty:
        print("No segments left after matching to chromosome sizes / drop list.", file=sys.stderr)
        sys.exit(1)

    # -------------------------
    # 3) Build global windows using chrom_df in chrom_size-file order
    # -------------------------
    win_df = build_global_windows_zero_based(chrom_df, args.bin_size, drop_list)
    if win_df.empty:
        print("No windows constructed from chrom sizes.", file=sys.stderr)
        sys.exit(1)

    # -------------------------
    # 4) Map segments -> windows
    # -------------------------
    df_cnv = segments_to_windows(df, win_df, args.bin_size, args.anchor, args.coords)
    if df_cnv.empty:
        print("No overlaps between segments and windows. Check coords/anchor/bin-size.", file=sys.stderr)
        sys.exit(1)

    # -------------------------
    # 5) Pivot to sample × window matrix
    #    Columns are reindexed to match win_df order ⇒ chromosome order follows chrom sizes file
    # -------------------------
    matrix = (
        df_cnv.pivot_table(index="sample", columns="segment", values="state", aggfunc="first")
        .reindex(columns=win_df["segment"].tolist(), fill_value=np.nan)
        .fillna(args.neutral)
    )
    matrix.index.name = None

    # -------------------------
    # 6) Sample labels and colors
    # -------------------------
    sample_labels = (
        df[["sample", "label"]]
        .drop_duplicates()
        .set_index("sample")
        .reindex(matrix.index)
    )
    sample_labels["label_str"] = sample_labels["label"].astype(str)
    label_encoder = LabelEncoder()
    sample_labels["label_code"] = label_encoder.fit_transform(sample_labels["label_str"])
    unique_labels = list(label_encoder.classes_)
    palette = sns.color_palette(args.palette, len(unique_labels))
    label_to_color = dict(zip(unique_labels, palette))
    row_colors = sample_labels["label_str"].map(label_to_color)
    row_colors.index.name = None
    row_colors.name = None

    # -------------------------
    # 7) Plot clustermap (columns are NOT clustered)
    # -------------------------
    sns.set(style="white")
    plt.rcParams.update({'font.family': 'DejaVu Sans'})

    g = sns.clustermap(
        matrix,
        cmap=["white"],          # patches will provide colors
        row_cluster=True,
        col_cluster=False,       # 👈 important to preserve column order
        row_colors=row_colors,
        xticklabels=False,
        yticklabels=False,
        figsize=tuple(args.figsize),
        dendrogram_ratio=(0.1, 0.01),
        cbar_pos=None
    )

    # Reorder samples according to row dendrogram
    clustered_order = g.dendrogram_row.reordered_ind
    clustered_samples = matrix.index[clustered_order]
    df_cnv = df_cnv[df_cnv["sample"].isin(clustered_samples)].copy()

    # Optional Y-axis labeling
    if not args.no_sampleLabel:
        g.ax_heatmap.set_yticks(np.arange(len(clustered_samples)) + 0.5)
        g.ax_heatmap.set_yticklabels(clustered_samples, fontsize=8)

    # CN state legend
    state_labels = [
        'homozygous deletion',
        'heterozygous deletion',
        'neutral',
        'gain',
        'amplification',
        'high-level amplification'
    ]
    state_colors = ['#77dd77', '#c2f0c2', 'white', '#f7c2c2', '#f77e7e', '#ff0000']
    cmap_custom = mcolors.ListedColormap(state_colors)
    norm_custom = mcolors.BoundaryNorm(boundaries=[1, 2, 3, 4, 5, 6, 7], ncolors=6)

    cb_ax = g.fig.add_axes([1.07, 0.60, 0.2, 0.02])
    cb = plt.colorbar(
        plt.cm.ScalarMappable(cmap=cmap_custom, norm=norm_custom),
        cax=cb_ax,
        orientation="horizontal"
    )
    cb.ax.tick_params(axis='x', length=0)
    for i, label in enumerate(state_labels):
        xpos = (i + 0.5) / 6
        if label == "homozygous deletion":
            xpos += 0.11
        elif label == "heterozygous deletion":
            xpos += 0.11
        elif label == "neutral":
            xpos += 0.02
        elif label == "amplification":
            xpos += 0.08
        elif label == "high-level amplification":
            xpos += 0.12
        cb_ax.text(xpos, 1.2, label, ha="center", va="bottom", fontsize=8,
                   rotation=45, transform=cb_ax.transAxes)

    # Sample label legend
    handles = [mpatches.Patch(color=palette[i], label=lab)
               for i, lab in enumerate(unique_labels)]
    g.ax_heatmap.legend(handles=handles,
                        loc="upper left",
                        bbox_to_anchor=(1.10, 0.5),
                        fontsize=12,
                        frameon=False)

    # Patch rectangles for each overlapping window
    col_index = {c: i for i, c in enumerate(matrix.columns)}
    sample_index = {s: i for i, s in enumerate(clustered_samples)}
    df_cnv["segment_id"] = df_cnv["segment"].astype(str)
    df_cnv.drop_duplicates(subset=["sample", "segment_id"], inplace=True)

    for _, r in df_cnv.iterrows():
        seg = r["segment_id"]
        j = col_index.get(seg, None)
        i = sample_index.get(r["sample"], None)
        if j is None or i is None:
            continue
        try:
            color = cmap_custom.colors[int(r["state"]) - 1]
        except Exception:
            color = "white"
        g.ax_heatmap.add_patch(Rectangle((j, i), 1, 1, color=color, linewidth=0))

    # Chromosome dividers and labels
    # counts/mids follow the order of win_df, which follows chrom_df, which follows chrom sizes file
    counts = win_df.groupby("chr", sort=False).size()
    cum = counts.cumsum()
    starts = cum.shift(fill_value=0)
    mids = (starts + cum) / 2

    for chr_name, pos in cum.items():
        g.ax_heatmap.axvline(x=pos, color="black", linewidth=0.3, zorder=100)

    g.ax_heatmap.set_xlabel("chromosome", labelpad=25)
    g.ax_heatmap.set_ylabel("")

    for pos, label in zip(mids, mids.index):
        g.ax_heatmap.text(pos, matrix.shape[0] + 0.3,
                          label, ha="center", va="top", fontsize=8)

    plt.savefig(args.output, bbox_inches="tight")
    print(f"[Saved] {args.output}")


if __name__ == "__main__":
    main()

