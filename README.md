# lncRNA in pathogenesis of focal segmental glomerulosclerosis -Wt1os

## Quality Control and Transcript Quantification

```bash
nextflow run nf-core/rnaseq -profile singularity --input samples.csv
--gtf data/Mus_musculus/mm10-2020-A-reference-sources/Mus_musculus.GRCm38.102.gtf.gz
--fasta data/Mus_musculus/mm10-2020-A-reference-sources/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
--outdir data/nf_RNAouts_salmon --pseudo_aligner salmon --extra_salmon_quant_args '--numBootstraps 100
--validateMappings --gcBias --seqBias' --aligner star_salmon --extra_star_align_args '--outFilterScoreMin 20'
```


## Reproducing Paper Results 

To reproduce results please use notebook `analysis.ipynb`


