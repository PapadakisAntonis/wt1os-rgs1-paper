# H3K27ac CUT&Tag: differential peaks and figures (glomeruli, GRCm38/mm10)
# Figures 2c, 4e and Supplementary Figure S2b-d.
# Inputs: cuttag.csv (consensus peak counts + peak annotation),
#         ReMap2022 mm10 BED.

library(tidyverse)
library(data.table)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)

cuttag_csv <- "/cellfile/cellnet/Antonis/AminaData/data/nf_CUTouts/cuttag_deseq.csv"
remap_bed  <- "/cellfile/cellnet/Antonis/remap2022_all_macs2_mm10_v1_0.bed.gz"
outdir     <- "/cellfile/cellnet/Antonis/Amina_results/CUTTAG"

sample_ids <- c("H3K27AC_CTRL_R1", "H3K27AC_CTRL_R2", "H3K27AC_CTRL_R3",
                "H3K27AC_KO_R1",   "H3K27AC_KO_R2",   "H3K27AC_KO_R3")

# --- Differential analysis from the consensus peak counts ---------------
csv <- read.csv(cuttag_csv)
csv$peakID <- as.character(csv$peakID)

counts <- csv %>% dplyr::select(all_of(sample_ids))
rownames(counts) <- csv$peakID

coldata <- data.frame(
  condition = factor(c(rep("wt", 3), rep("ko", 3)), levels = c("wt", "ko")),
  row.names = sample_ids
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

res <- results(dds, contrast = c("condition", "ko", "wt"),
               independentFiltering = FALSE, cooksCutoff = FALSE) %>%
  as.data.frame() %>%
  rownames_to_column("peakID")

# peak table
deseq_file <- csv %>%
  dplyr::select(peakID, seqnames, start, end, width, strand, annotation,
                distanceToTSS, geneId, transcriptId, ENSEMBL, SYMBOL, GENENAME) %>%
  left_join(dplyr::select(res, peakID, baseMean, log2FoldChange,
                          lfcSE, stat, pvalue, padj),
            by = "peakID")

# --- PCA (Fig 2c) -------------------------------------------------------
pca <- plotPCA(vsd, intgroup = "condition") +
  theme_minimal() +
  ggtitle("PCA of H3K27ac Cut&Tag Samples")
pca

# --- GO/KEGG of genes at differential peaks (Fig S2b-d) -----------------
sig_peaks <- deseq_file %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
genes <- sig_peaks$SYMBOL %>% unique() %>% na.omit()
gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

go  <- enrichGO(gene_ids, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go1 <- enrichGO(gene_ids, OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go2 <- enrichGO(gene_ids, OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05)
#kegg <- enrichKEGG(gene_ids, organism = "mmu", pAdjustMethod = "BH", qvalueCutoff = 0.05)

dotplot(go,  showCategory = 10, title = "GO Enrichment for H3K27ac DARs") + theme_minimal()
dotplot(go1, showCategory = 10, title = "GO Enrichment for H3K27ac DARs, MF") + theme_minimal()
dotplot(go2, showCategory = 10, title = "GO Enrichment for H3K27ac DARs, CC") + theme_minimal()
#dotplot(kegg, showCategory = 10, title = "KEGG Enrichment for H3K27ac DARs") + theme_minimal()

# --- TFs with a ReMap peak at the Rgs1 locus ----------------------------
remap <- fread(remap_bed, header = FALSE,
               col.names = c("chr", "start", "end", "tf_info", "score",
                             "strand", "thick_start", "thick_end", "RGB"))
rgs1_tfs <- remap %>%
  filter(chr == "chr1", end >= 144248519, start <= 144250448) %>%
  mutate(tf_name = str_split_fixed(tf_info, "\\.", 3)[, 2])
rgs1_binding_tfs <- unique(str_to_title(rgs1_tfs$tf_name))

# --- H3K27ac at Rgs1-regulating TFs (Fig 4e) ----------------------------
tf_gene_peaks_sig <- deseq_file %>%
  filter(SYMBOL %in% rgs1_binding_tfs, padj < 0.05, abs(log2FoldChange) > 1)

vst_tf_peaks <- assay(vsd)[tf_gene_peaks_sig$peakID, ]
rownames(vst_tf_peaks) <- paste0(tf_gene_peaks_sig$SYMBOL, "_", tf_gene_peaks_sig$peakID)
scaled_vst <- t(scale(t(vst_tf_peaks)))

sample_annotation <- data.frame(
  Condition = factor(c(rep("WT", 3), rep("KO", 3)), levels = c("WT", "KO")),
  row.names = colnames(scaled_vst)
)

pheatmap(scaled_vst,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = sample_annotation,
         annotation_colors = list(Condition = c("WT" = "#2E86C1", "KO" = "#E74C3C")),
         main = "H3K27ac at TF genes that regulate Rgs1",
         fontsize_row = 5,
         #filename = file.path(outdir, "heatmap_rgs1_tf_genes.png"),
         width = 8,
         height = 10)
