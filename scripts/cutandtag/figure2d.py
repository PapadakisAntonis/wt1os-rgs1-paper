#!/usr/bin/env python3
# Figure 2d: average the per-sample H3K27ac gene-body matrix into CTRL and KO
# and render the heatmap.
# Input:  all_genes.computeMatrix.mat.gz  (nf-core/cutandrun computeMatrix output)
# Output: matrix_aggregated.gz, heatmap_aggregated.pdf

import gzip
import json
import subprocess
import numpy as np

infile    = "all_genes.computeMatrix.mat.gz"
outmatrix = "matrix_aggregated.gz"
outpdf    = "heatmap_aggregated.pdf"

# read header params, per-row gene annotation (first 6 cols) and values
with gzip.open(infile, "rt") as f:
    params = json.loads(f.readline().strip()[1:])
    annotations, values = [], []
    for line in f:
        parts = line.rstrip("\n").split("\t")
        annotations.append(parts[:6])
        values.append([np.nan if x == "nan" else float(x) for x in parts[6:]])

values = np.array(values)

# CTRL/KO column blocks and bins-per-sample, taken from the header
labels   = params["sample_labels"]
ctrl_idx = [i for i, l in enumerate(labels) if "CTRL" in l]
ko_idx   = [i for i, l in enumerate(labels) if "KO"   in l]
n_bins   = params["sample_boundaries"][1] - params["sample_boundaries"][0]

def group_mean(idx):
    return np.nanmean([values[:, i * n_bins:(i + 1) * n_bins] for i in idx], axis=0)

aggregated = np.hstack([group_mean(ctrl_idx), group_mean(ko_idx)])

params["sample_labels"]     = ["CTRL", "KO"]
params["sample_boundaries"] = [0, n_bins, 2 * n_bins]

with gzip.open(outmatrix, "wt") as f:
    f.write("@" + json.dumps(params) + "\n")
    for anno, row in zip(annotations, aggregated):
        f.write("\t".join(anno + [str(x) for x in row]) + "\n")

subprocess.run([
    "plotHeatmap",
    "-m", outmatrix,
    "-o", outpdf,
    "--colorMap", "RdBu_r",
    "--whatToShow", "heatmap and colorbar",
    "--heatmapHeight", "14",
    "--heatmapWidth", "8",
    "--samplesLabel", "CTRL", "KO",
], check=True)
