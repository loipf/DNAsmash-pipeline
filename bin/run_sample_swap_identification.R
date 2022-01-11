#!/usr/bin/env Rscript

############################################
### gather all files

pacman::p_load(data.table, maftools, pheatmap, R.utils)

args = commandArgs(trailingOnly=TRUE)
NUM_THREADS = args[1]
MIN_SNP_READ_DEPTH = args[2]

MAPPED_FILES = list.files(path=".", pattern="*.bam$", full.names=T,recursive=T)


############################################
### run sample matching based on SMaSH

### bugfixed from:
### https://github.com/PoisonAlien/maftools/blob/master/R/sample_swap.R
sampleSwaps_edited = function (bams = NULL, build = "hg19", prefix = NULL, add = TRUE, 
          min_depth = 30, ncores = 4, ...) 
{
  if (length(bams) < 2) {
    stop("Needs 2 or more BAM files!")
  }
  build = match.arg(arg = build, choices = c("hg19", "hg38"))
  if (build == "hg19") {
    snps = system.file("extdata", "hg19_smash_snps.tsv.gz", 
                       package = "maftools")
  }
  else {
    snps = system.file("extdata", "hg38_smash_snps.tsv.gz", 
                       package = "maftools")
  }
  snps = data.table::fread(input = snps, sep = "\t")
  if (!is.null(prefix)) {
    if (add) {
      snps$chr = paste(prefix, snps$chr, sep = "")
    }
    else {
      snps$chr = gsub(pattern = prefix, replacement = "", 
                      x = snps$chr, fixed = TRUE)
    }
  }
  rc = bamreadcounts(bam = bams, loci = snps, nthreads = ncores ) #, ...)
  snps[, `:=`(id, paste0(chr, ":", start))]
  cat("Summarizing allele frequncy table..\n")
  rc_af = lapply(rc, function(x) {
    xrc = merge(x, snps[, .(id, ref, alt)], by.x = "loci", 
                by.y = "id")
    xrc_acounts = apply(xrc, 1, function(x) {
      ref_allele = x["ref"]
      ref_rc = 0
      ref_rc = switch(ref_allele, A = x["A"], T = x["T"], 
                      G = x["G"], C = x["C"])
      alt_allele = x["alt"]
      alt_rc = 0
      alt_rc = switch(alt_allele, A = x["A"], T = x["T"], 
                      G = x["G"], C = x["C"])
      vaf_tbl = data.table::data.table(ref_rc = as.numeric(ref_rc), 
                                       alt_rc = as.numeric(alt_rc), loci = x["loci"])
      vaf_tbl[, `:=`(vaf, alt_rc/(ref_rc + alt_rc))]
      vaf_tbl
    })
    xrc = merge(xrc, data.table::rbindlist(l = xrc_acounts), 
                by = "loci")
    xrc[, .(loci, ref_rc, alt_rc, vaf)]
  })
  rc_bind = data.table::rbindlist(l = rc_af, idcol = "sample")
  rc_bind = rc_bind[!is.nan(vaf)]
  if (nrow(rc_bind) == 0) {
    stop("Zero SNPs to analyze!")
  }
  rc_bind = rc_bind[, `:=`(total, ref_rc + alt_rc)][total > 
                                                      min_depth]
  rc_df = data.table::dcast(data = rc_bind, loci ~ sample, 
                            value.var = "vaf", fill = NA)
  data.table::setDF(x = rc_df, rownames = rc_df$loci)
  rc_df$loci = NULL
  rc_df = rc_df[complete.cases(rc_df), ]
  cat("Performing pairwise comparison..\n")
  samples_no_snps_found = setdiff(names(rc_af), colnames(rc_df)) ## edit
  rc_af_snps_found = rc_af[colnames(rc_df)] ## edit
  sample_matches = parallel::mclapply(seq_along(rc_af_snps_found), function(idx) {
    print(idx)
    x = rc_af_snps_found[[idx]]
    if (idx == length(rc_af_snps_found)) {
      return(NULL)
    }
    else {
      rest_samps = seq_along(rc_af_snps_found)[(idx + 1):(length(rc_af_snps_found))]
    }
    cat("Comparing", names(rc_af_snps_found)[idx], "against rest..\n")
    x_compare = lapply(rest_samps, function(rest_idx) {
      print(rest_idx)
      y = rc_af_snps_found[[rest_idx]]
      concordant_snps = lapply(seq_along(1:nrow(y)), function(row_idx) {
        fisher.test(x = matrix(c(x[row_idx, ref_rc], 
                                 x[row_idx, alt_rc], y[row_idx, ref_rc], y[row_idx, 
                                                                           alt_rc]), nrow = 2, ncol = 2))$p.value
      })
      concordant_snps = table(unlist(lapply(concordant_snps, 
                                            function(x) ifelse(x < 0.01, no = "concordant", 
                                                               yes = "discordant"))))
      cor_coef = cor.test(x$vaf, y$vaf)$estimate
      data.table::data.table(X_bam = names(rc_af_snps_found)[idx], 
                             Y_bam = names(rc_af_snps_found)[rest_idx], concordant_snps = concordant_snps[[1]], 
                             discordant_snps = concordant_snps[[2]], fract_concordant_snps = prop.table(concordant_snps)[[1]], 
                             cor_coef = cor_coef)
    })
    x_compare = data.table::rbindlist(l = x_compare)
    x_compare[, `:=`(XY_possibly_paired, ifelse(test = fract_concordant_snps >= 
                                                  0.8 & cor_coef >= 0.9, yes = "Yes", no = "No"))]
    x_compare
  }, mc.cores = ncores)
  sample_matches = sample_matches[1:(length(sample_matches) - 
                                       1)]
  pos_mathces = lapply(sample_matches, function(sample_pair) {
    if (nrow(sample_pair) > 0) {
      unique(unlist(sample_pair[XY_possibly_paired == "Yes", 
                                .(X_bam, Y_bam)], use.names = FALSE))
    }
  })
  pos_mathces = pos_mathces[which(lapply(pos_mathces, function(x) length(x) > 
                                           0) == TRUE)]
  cat("Done!\n")
  list(AF_table = invisible(rc_df), SNP_readcounts = data.table::rbindlist(l = rc_af, 
                                                                           idcol = "BAM", use.names = TRUE, fill = TRUE), pairwise_comparison = data.table::rbindlist(sample_matches)[order(XY_possibly_paired, 
                                                                                                                                                                                            decreasing = TRUE)],
       BAM_matches = pos_mathces, samples_no_snps_found=samples_no_snps_found)
}


# sample_swap_results = sampleSwaps(bams = MAPPED_FILES, build = "hg38", add=F, prefix="chr", min_depth=MIN_SNP_READ_DEPTH, ncores=NUM_THREADS)
sample_swap_results = sampleSwaps_edited(bams = MAPPED_FILES, build = "hg38", add=F, prefix="chr", min_depth=MIN_SNP_READ_DEPTH, ncores=NUM_THREADS)

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


