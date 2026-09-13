# CUT&Tag processing — nf-core/cutandrun

H3K27ac CUT&Tag on isolated glomeruli, *Wt1os*-/- vs WT, with IgG controls and
an *E. coli* spike-in for normalisation.

## Version

The original run used `-r master`, resolved on 2024-06-28 to commit
`2727de382fb7df8f37fc4c6dbde99c979be04e77`, which is tag **3.2.1**.
`master` has moved since, so reproduce with the tag, not the branch:

```bash
nextflow run nf-core/cutandrun -r 3.2.1 ...
```

| | |
|---|---|
| Pipeline | nf-core/cutandrun 3.2.1 (`2727de38`) |
| Nextflow | 23.04.1 build 5866 |
| Profile | singularity |
| Run | 2024-06-28 to 2024-06-30 (2d 1h 21m) |

## Command

Paths below are relative to the repository; the original run used absolute
paths under the operator's scratch directory.

```bash
nextflow run nf-core/cutandrun -r 3.2.1 -profile singularity \
  --input metadata/samples_cutandrun.csv \
  --bowtie2 data/genome/index/bowtie2 \
  --gtf data/genome/Mus_musculus.GRCm38.102.gtf \
  --fasta data/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
  --spikein_bowtie2 data/spikein/Escherichia_coli_K_12_MG1655/Bowtie2Index \
  --spikein_fasta data/spikein/Escherichia_coli_K_12_MG1655/genome.fa \
  --outdir data/nf_CUTouts
```

## Reference

| | |
|---|---|
| Genome | GRCm38 / mm10, `dna_sm` primary assembly |
| Annotation | Ensembl release 102 GTF |
| Spike-in | *E. coli* K-12 MG1655, NCBI 2001-10-15 (from the mm10-2020-A reference bundle) |

The same GRCm38.102 annotation is used for the RNA-seq (see the top-level
README), so gene identifiers are consistent across the two assays.

## Non-default parameters

| Parameter | Value |
|---|---|
| `minimum_alignment_q_score` | 20 |
| `normalisation_binsize` | 50 |
| `macs_gsize` | 2.7e9 |
| `replicate_threshold` | 1 |
| `min_frip_overlap` | 0.20 |
| `min_peak_overlap` | 0.20 |

`replicate_threshold 1` means a peak is retained if it is called in a single
replicate, rather than requiring reproducibility across replicates. This is
permissive; the downstream DESeq2 step in `scripts/cutandtag/` is what
establishes differential enrichment.

## Downstream

The consensus peak counts from this run, annotated with ChIPseeker, become
`processed/cuttag_deseq.csv.gz`. The deepTools `computeMatrix` output
(`all_genes.computeMatrix.mat.gz`) is the input to `figure2d.py`.
