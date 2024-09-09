#!/usr/bin/env python3
# File: collect_tre.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 20, 2024
# Updated: Jul 28, 2024
import os
import math
from multiprocessing import Pool
from pathlib import Path
from argparse import ArgumentParser
import anndata as ad

import pysam as ps
import pandas as pd
import polars as pl
import genomicranges as gr
import matplotlib.pyplot as plt


REQUIRED_COLS = ["seqnames", "starts", "ends"]


def fetch_region_info(peak_db, region):
    for k, v in peak_db.items():
        raw, full, info = [str(v[x]) for x in ["raw", "full", "info"]]
        raw_regions = ps.TabixFile(raw, parser=ps.asBed()).fetch(region)
        sig_regions = ps.TabixFile(full, parser=ps.asBed()).fetch(region)
        info_regions = ps.TabixFile(info, parser=ps.asBed()).fetch(region)

    return raw_regions, sig_regions, info_regions


def get_cli():
    parser = ArgumentParser("")
    parser.add_argument("-d", "--tre-dir", required=True, help="The folder storing TRE results.")
    parser.add_argument("-g", "--gff-file", required=False, help="The GFF file")
    parser.add_argument("-f", "--fasta-file", required=False, help="The FASTA file.")

    return parser.parse_args()


def read_gff(infile):
    gtf_gr = gr.read_gtf(infile)

    return gtf_gr

# Utility functions
def parse_region(region):
    chrom, start, end = None, None, None
    if ":" in region:
        chrom, start_end = region.replace(",", "").split(":")
        if "-" in start_end:
            start, end = start_end.split("-")
            start, end = int(start), int(end)
        else:
            start, end = int(start_end), int(math.inf)

    return chrom, start, end


def plot_region(peak, region, save_to, label_col="SRR_id", fig_size=(10, 6)):
    seqnames, starts, ends = parse_region(region)
    selected_peaks = (peak.filter(pl.col("seqnames") == seqnames, pl.col("starts") >= starts, pl.col("ends") <= ends)
                      .with_columns(pl.int_range(pl.len()).alias("y_position")))

    y_lab = selected_peaks.select(pl.col(label_col)).to_numpy().squeeze()
    y_pos = selected_peaks.select(pl.col("y_position")).to_numpy().squeeze()
    x_pos_min = selected_peaks.select(pl.col("starts")).to_numpy().squeeze()
    x_pos_max = selected_peaks.select(pl.col("ends")).to_numpy().squeeze()

    fig, axe = plt.subplots(1, 1, figsize=fig_size)
    axe.hlines(y_pos, x_pos_min, x_pos_max, color="0", linewidth=1)
    for _y, _x, _xx in zip(y_pos, x_pos_min, x_pos_max):
        axe.scatter(_x, _y, color="0", marker="|")
        axe.scatter(_xx, _y, color="0", marker="|")
    axe.set_title(region)
    axe.set_yticks(y_pos)
    axe.set_yticklabels(y_lab)

    fig.savefig(save_to)
    fig.clear()
    plt.close(fig)


# 1. Load all full peaks and merge just base on the positions
def read_peak_bed(infile, new_columns=None, extra_info={}):
    pk_pl = pl.read_csv(infile, separator="\t", has_header=False, new_columns=new_columns)

    if isinstance(extra_info, dict):
        for k, v in extra_info.items():
            pk_pl = pk_pl.with_columns(pl.lit(v).alias(k))
    elif isinstance(extra_info, pl.Series):
        pass
    elif isinstance(extra_info, pd.Series):
        pass

    return pk_pl


def merge_intervals(intervals, max_gap=0):
    groups = [0]
    for idx, item in enumerate(intervals[1:]):
        last_chrom, _, last_end = intervals[idx]
        cur_chrom, cur_start, _ = item
        on_same_chrom = cur_chrom == last_chrom
        is_overlapped = cur_start <= last_end - max_gap

        if on_same_chrom and is_overlapped:
            groups.append(groups[-1])
        else:
            groups.append(groups[-1] + 1)

    return groups, intervals


# 2. Estimate how complex of the merged regions
# 2.1 More than one significant peak from the same sample.
def complexity(peak_tab, id_col="interval"):
    peak_tab.group_by(id_col)


# 3. Refine the merged regions based on the peak information
def harmonize_peaks(peak_list, method="expand", min_overlap=100, max_gap=0, bin_size=100):
    """
    Harmonize peaks by merging peaks from multiple samples.

    Args:
      peak_list: list, a list of polars.DataFrame.
      method: str, method to merge peaks.

    Returns:
      GenomicRanges
    """
    inter_cols = ["seqnames", "starts", "ends"]
    peak_pldf = pl.concat(peak_list).sort(pl.col(["seqnames", "starts", "ends"]))
    intervals = peak_pldf.select(pl.col(inter_cols)).to_numpy().tolist()
    groups, _ = merge_intervals(intervals, max_gap=max_gap)
    peak_pldf = peak_pldf \
        .with_columns(group=pl.Series(values=groups)) \
        .with_columns(interval=pl.col("seqnames") + ":" + pl.col("starts").min().over("group").cast(str) + "-" + pl.col("ends").max().over("group").cast(str))

    return peak_pldf


def main():
    pass

proj_dir = Path("/home/zzhang/Documents/projects/wp_ipf/outputs/analysis/tre_identification/tre_results")
sig_peak_col_names = REQUIRED_COLS + ["score", "p_value", "centroid"]
samples = [x for x in os.listdir(proj_dir) if x.startswith("SRR")]
sig_peak_file_path = [proj_dir / f'{x}/{x}.dREG.peak.full.bed.gz' for x in samples if (proj_dir / f'{x}/{x}.dREG.peak.raw.bed.gz').exists()]
sig_peak_list = [read_peak_bed(p, new_columns=sig_peak_col_names, extra_info={"SRR_id": p.parent.stem}) for p in sig_peak_file_path]
sig_peaks = harmonize_peaks(sig_peak_list)
sig_peaks.group_by("interval").agg(pl.col("SRR_id").count().alias("n_samples"), pl.col("SRR_id").sort().str.concat(";").alias("samples"), pl.col("seqnames").first(), pl.col("starts").min(), pl.col("ends").max()).sort(pl.col(["seqnames", "starts", "ends"]))
sig_peaks.group_by("interval").agg((pl.col("SRR_id").count() != pl.col("SRR_id").unique().count()).alias("Fusion_peak")).select("Fusion_peak").sum()
sig_peaks.group_by("interval").agg((pl.col("SRR_id").count() != pl.col("SRR_id").unique().count()).alias("Fusion_peak"))

infp_tab_file_path = [proj_dir / f'{x}/{x}.dREG.infp.bed.gz' for x in samples if (proj_dir / f'{x}/{x}.dREG.infp.bed.gz').exists()]
infp_tab_list = [pl.read_csv(p, separator="\t", new_columns=["seqnames", "starts", "ends", "score", "Unknown"], has_header=False) for p in sig_peak_file_path[:10]]
infp_tabs = pl.concat(infp_tab_list).sort(pl.col(["seqnames", "starts", "ends"]))

bed_dir = [proj_dir / f'{x}' for x in samples if (proj_dir / f'{x}/{x}.dREG.peak.raw.bed.gz').exists()]
peak_db = {
    x: {
        "raw": proj_dir / f'{x}/{x}.dREG.peak.raw.bed.gz',
        "full": proj_dir / f'{x}/{x}.dREG.peak.full.bed.gz',
        "info": proj_dir / f'{x}/{x}.dREG.infp.bed.gz'
    }
    for x in samples if (proj_dir / f'{x}/{x}.dREG.peak.raw.bed.gz').exists()
}



if __name__ == "__main__":
    main()
