# Wt1os / Rgs1 in focal segmental glomerulosclerosis

Analysis code for the *Wt1os* glomerular RNA-seq and H3K27ac CUT&Tag work
(manuscript in revision).

The RNA-seq analysis was written by **Amina Lemsara** and originally released at
[AminaLEM/Wt1os_FSGS](https://github.com/AminaLEM/Wt1os_FSGS), from which this
repository is derived. Licensed GPL-3.0, inherited from the original.

## Layout

```
scripts/rnaseq/       DESeq2 differential expression, PCA, GSEA
scripts/cutandrun/    H3K27ac differential peaks, PCA, GO, figures
nextflow/             pipeline configuration
metadata/             sample sheets
processed/            committed intermediates (consensus peak table, TF list)
envs/                 package versions
```

Raw data and large external files are not in this repository; `data/` is
git-ignored. Derived tables small enough to ship live in `processed/`.

## Processing

RNA-seq quantification, nf-core/rnaseq:

```bash
nextflow run nf-core/rnaseq -profile singularity --input samples.csv \
  --gtf data/Mus_musculus/mm10-2020-A-reference-sources/Mus_musculus.GRCm38.102.gtf.gz \
  --fasta data/Mus_musculus/mm10-2020-A-reference-sources/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz \
  --outdir data/nf_RNAouts_salmon --pseudo_aligner salmon \
  --extra_salmon_quant_args '--numBootstraps 100 --validateMappings --gcBias --seqBias' \
  --aligner star_salmon --extra_star_align_args '--outFilterScoreMin 20'
```

Downstream analysis is R 4.2.3; see `scripts/rnaseq/analysis.ipynb`.

## Licence

GPL-3.0, see `LICENSE`.
