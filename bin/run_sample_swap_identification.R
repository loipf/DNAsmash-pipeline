#!/usr/bin/env Rscript

############################################
### gather all files

pacman::p_load(data.table, maftools, pheatmap, R.utils)

args = commandArgs(trailingOnly=TRUE)
NUM_THREADS = args[1]

MAPPED_FILES = list.files(path=".", pattern="*.bam$", full.names=T,recursive=T)


############################################
### run sample matching based on SMaSH

sample_swap_results = sampleSwaps(bams = MAPPED_FILES, build = "hg38", add=F, prefix="chr", min_depth=20, ncores=NUM_THREADS)

### output
saveRDS(sample_swap_results, "sample_swap_results.rds")

sink("matching_samples.txt")
if(length(sample_swap_results$BAM_matches) >0){
  writeLines(unlist(lapply(sample_swap_results$BAM_matches, paste, collapse=" ")), sep="\n\n")
}
sink()

### plot

# RColorBrewer::brewer.pal(n = 7, name = "Blues")
BLUE_COLOR_PALETTE = c("#EFF3FF","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#084594")

tmp_AF_table = sample_swap_results$AF_table
colnames(tmp_AF_table) = substr(colnames(tmp_AF_table), start=1, stop=20 )
pdf("sample_swap_corr_heatmap.pdf", width = 12.1, height = 12)
print(pheatmap(cor(tmp_AF_table), breaks = seq(0, 1, 0.01),
               main="AF correlation heatmap",
               cluster_col=T, cluster_row=T,
               color = colorRampPalette(BLUE_COLOR_PALETTE)(100),
               fontsize_col = 7, fontsize_row = 7,
               show_rownames=T, show_colnames=T))
dev.off()

