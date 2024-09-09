#!/usr/bin/env python3
import time
from pathlib import Path
from multiprocessing import Pool
import random

from tqdm import tqdm
import pysam as psm
import pyarrow as paw
import h5py as hpy
import sqlite3 as sql
import polars as pls
import pyarrow.parquet as par
import pandas as pds
import numpy as npy

def chunknize(start: int, end: int, step: int = 500000):
    starts = list(range(start, end + 1, step))
    ends = [x - 1 for x in starts[1:]] + [end]
    return zip(starts, ends)


def create_indice_tab(step: int = 500000):
    CHROM_LENGTH = {
         "chr1": 247199719,  "chr2": 242751149,  "chr3": 199446827,  "chr4": 191263063,  "chr5": 180837866,
         "chr6": 170896993,  "chr7": 158821424,  "chr8": 146274826,  "chr9": 140442298, "chr10": 135374737,
        "chr11": 134452384, "chr12": 132289534, "chr13": 114127980, "chr14": 106360585, "chr15": 100338915,
        "chr16": 88822254,  "chr17": 78654742,  "chr18": 76117153,  "chr19": 63806651,  "chr20": 62435965,
        "chr21": 46944323,  "chr22": 49528953,   "chrX": 154913754,  "chrY": 57741652,
    }

    indice_list = []
    for k, v in CHROM_LENGTH.items():
        per_chunk = [[k, s, e, i + 1] for i, (s, e) in enumerate(chunknize(1, v, step))]
        indice_list.extend(per_chunk)

    new_names = {"column_0": "seqnames", "column_1": "chunk_start", "column_2": "chunk_end", "column_3": "chunk_index"}
    return pls.DataFrame(indice_list).rename(new_names)


def merge_intervals(intervals, max_gap=0):
    groups = [0]
    for idx, item in enumerate(intervals[1:]):
        last_chrom, _, last_end = intervals[idx]
        cur_chrom, cur_start, _ = item
        on_same_chrom = cur_chrom == last_chrom
        is_overlapped = cur_start <= last_end + max_gap

        if on_same_chrom and is_overlapped:
            groups.append(groups[-1])
        else:
            groups.append(groups[-1] + 1)

    return groups, intervals


def bed_to_any(bed_file: Path, out_file: Path, source: str = "raw", chunk_size=500000):
    grange_cols = ["seqnames", "starts", "ends"]
    if source in ["raw", "full"]:
        new_cols = grange_cols + ["score", "p_value"]
    elif source == "infp":
        new_cols = grange_cols + ["score"]
    else:
        raise ValueError(f"Unsupported source: {source}")

    # Create indice table
    bed_idx = create_indice_tab(chunk_size)

    # Data with spliced into chunks which matches indice table.
    bed_tab = pls.read_csv(
        bed_file, has_header=False, new_columns=new_cols, separator="\t"
    ).with_columns(# Splite data into chunks by floor-division
        pls.col("starts").floordiv(chunk_size).mul(chunk_size).add(1).alias("chunk_start"),
        pls.col("ends").floordiv(chunk_size).add(1).mul(chunk_size).alias("chunk_end"),
    ).join(bed_idx, on=["seqnames", "chunk_start", "chunk_end"], how="left").sort(grange_cols)

    out_type = out_file.suffix # Decide the output type
    if out_type == ".sql":
        # Create SQL database. Including a table of indice, and a hand of tables of data.
        print(f"[I]: Creating database file: {out_file}")
        with sql.connect(out_file) as db:
            db.execute("DROP TABLE IF EXISTS indices")
            db.execute("CREATE TABLE indices(seqnames VARCHAR(8), starts INT, ends INT, chunk_index INT)")
            db.executemany("INSERT INTO indices VALUES(?, ?, ?, ?);", bed_idx.unique().to_numpy())
            db.execute("DROP TABLE IF EXISTS sample_1")

            if source in ["raw", "full"]:
                pc_data = bed_tab.select("seqnames", "starts", "ends", "chunk_index", "score", "p_value").to_numpy()
                db.execute("CREATE TABLE sample_1(seqnames VARCHAR(8), starts INT, ends INT, chunk_index INT, score REAL, p_value REAL)")
                db.executemany(f"INSERT INTO sample_1 VALUES(?, ?, ?, ?, ?, ?)", pc_data)
            elif source == "infp":
                pc_data = bed_tab.select("seqnames", "starts", "ends", "chunk_index", "score").to_numpy()
                db.execute("CREATE TABLE sample_1(seqnames VARCHAR(8), starts INT, ends INT, chunk_index INT, score REAL)")
                db.executemany(f"INSERT INTO sample_1 VALUES(?, ?, ?, ?, ?)", pc_data)
            db.commit()
        print(f"[I]: Check the database file: {out_file}")
    elif out_type == ".hdf5":
        # Create HDF5 file. Including a group of indice and a group of data.
        # Both indice and data were then splited into dataset per chromsome.
        print(f"[I]: Creating hdf5 file: {out_file}")
        all_seqnames = bed_tab["seqnames"].unique().to_list()
        with hpy.File(out_file, "w") as hdf5_file:
            idx_grp = hdf5_file.create_group("indices")
            dat_grp = hdf5_file.create_group("sample_1")
            for pc in all_seqnames:
                pc_idx = bed_idx.filter(pls.col("seqnames") == pc).select("chunk_start", "chunk_end", "chunk_index")
                idx_grp.create_dataset(pc, data=pc_idx, dtype=npy.int32)

                if source in ["raw", "full"]:
                    pc_dat = bed_tab.filter(pls.col("seqnames") == pc).select("starts", "ends", "chunk_index", "score", "p_value")
                    dat_grp.create_dataset(f"{pc}_index", data=pc_dat.select("starts", "ends", "chunk_index"), dtype="int32")
                    dat_grp.create_dataset(f"{pc}_data", data=pc_dat.select("score", "p_value"), dtype="f")
                else:
                    pc_dat = bed_tab.filter(pls.col("seqnames") == pc).select("starts", "ends", "chunk_index", "score")
                    dat_grp.create_dataset(f"{pc}_index", data=pc_dat.select("starts", "ends", "chunk_index"), dtype="int32")
                    dat_grp.create_dataset(f"{pc}_data", data=pc_dat.select("score"), dtype="f")
        print(f"[I]: Check the hdf5 file: {out_file}")
    else:
        raise ValueError(f"Unsupported output type: {out_type}")


def test_random_access_bed(in_file, n_access=1000):
    """Test the random access speed via pysam on BED files."""
    m1, m2 = None, None
    print("[I]: testing random access opening and closing for each access ...")
    for idx in tqdm(range(n_access)):
        random.seed(idx)
        end = random.randint(1, 1000000)
        with psm.TabixFile(in_file, parser=psm.asBed()) as f:
            m1 = [x for x in f.fetch("chr10", 1, end)]

    print("[I]: testing random access opening at the beginning and closing at the end...")
    with psm.TabixFile(in_file, parser=psm.asBed()) as f:
        for idx in tqdm(range(n_access)):
            random.seed(idx)
            end = random.randint(1, 1000000)
            m2 = [x for x in f.fetch("chr10", 1, end)]

    return m1, m2


def test_random_access_hdf5(in_file, n_access=100000):
    """Test the random access speed via polars on hdf5 files."""
    seqname = "chr1"
    chunk_tab = None
    with hpy.File(in_file, "r") as hdf5_file:
        chunk_tab = hdf5_file.get(f"indices/{seqname}")[:, :]

        for idx in tqdm(range(n_access)):
            random.seed(idx)
            end = random.randint(1, 1000000)
            pls.DataFrame(chunk_tab).filter(1 <= pls.col("column_0"), pls.col("column_1") <= end).select("column_2")

        return chunk_tab


def test_random_access_sql(in_file, n_access=100000):
    """Test the random access speed via polars on sql files."""
    pass


def main():
    pass

proj_dir = Path("~/Documents/projects/wp_ipf").expanduser()

bed_file = proj_dir / "outputs/analysis/tre_identification/tre_results/SRR13523284/SRR13523284.dREG.peak.raw.bed.gz"
out_file = proj_dir / "temps/SRR13523284.raw.sql"
bed_to_any(bed_file, out_file, "raw", chunk_size=1024 * 1024 * 2)
out_file = proj_dir / "temps/SRR13523284.raw.hdf5"
bed_to_any(bed_file, out_file, "raw", chunk_size=1024 * 1024 * 2)

bed_file = proj_dir / "outputs/analysis/tre_identification/tre_results/SRR13523284/SRR13523284.dREG.peak.full.bed.gz"
out_file = proj_dir / "temps/SRR13523284.full.sql"
bed_to_any(bed_file, out_file, "full", chunk_size=1024 * 1024 * 2)
out_file = proj_dir / "temps/SRR13523284.full.hdf5"
bed_to_any(bed_file, out_file, "full", chunk_size=1024 * 1024 * 2)

bed_file = proj_dir / "outputs/analysis/tre_identification/tre_results/SRR13523284/SRR13523284.dREG.infp.bed.gz"
out_file = proj_dir / "temps/SRR13523284.infp.sql"
bed_to_any(bed_file, out_file, "infp", chunk_size=1024 * 1024 * 2)
out_file = proj_dir / "temps/SRR13523284.infp.hdf5"
bed_to_any(bed_file, out_file, "infp", chunk_size=1024 * 1024 * 2)

m1, m2 = test_random_access_bed(str(bed_file))

bed_file = proj_dir / "outputs/analysis/tre_identification/tre_results/SRR13523284/SRR13523284.dREG.peak.raw.bed.gz"
out_file = proj_dir / "temps/SRR13523284.raw.hdf5"
bed_to_any(bed_file, out_file, "raw", chunk_size=1024 * 1024 * 2)

chunk_tab = test_random_access_hdf5(str(proj_dir / "temps/SRR13523284.raw.hdf5"))

hdf5_file = hpy.File(str(proj_dir / "temps/SRR13523284.raw.hdf5"), "r")
chunk_tab = hdf5_file.get("indices/chr1")
pls.DataFrame(chunk_tab[:, :])
dat_tab = hdf5_file.get("sample_1/chr1_data")
pls.DataFrame(dat_tab[:200, :])
idx_tab = hdf5_file.get("sample_1/chr1_index")
pls.DataFrame(idx_tab[:200, :])
hdf5_file.close()


if __name__ == "__main__":
    main()
