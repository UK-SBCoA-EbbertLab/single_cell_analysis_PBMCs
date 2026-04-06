#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SoupX)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: soupx_pipseq.R <in_dir> <out_prefix>")

in_dir <- args[1]
out_prefix <- args[2]

# === SETUP LOGGING ===
log_file <- paste0(out_prefix, "_soupx.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("SoupX log will be saved to:", log_file, "\n")
cat(strrep("=", 60), "\n")
cat("SoupX Ambient RNA Correction\n")
cat(strrep("=", 60), "\n\n")

cat("Loading data from:", in_dir, "\n")

# Load matrices (genes × barcodes)
toc <- readMM(file.path(in_dir, "toc.mtx"))
tod <- readMM(file.path(in_dir, "tod.mtx"))

genes <- readLines(file.path(in_dir, "genes.tsv"))
cell_barcodes <- readLines(file.path(in_dir, "barcodes.tsv"))
empty_barcodes <- readLines(file.path(in_dir, "ambient_barcodes.tsv"))

# Set row/column names
rownames(toc) <- genes
colnames(toc) <- cell_barcodes
rownames(tod) <- genes
colnames(tod) <- empty_barcodes

cat("\nMatrix dimensions:\n")
cat("  toc (cells):   ", nrow(toc), "genes x", ncol(toc), "cells\n")
cat("  tod (empties): ", nrow(tod), "genes x", ncol(tod), "empties\n")

# Create SoupChannel
cat("\nCreating SoupChannel and estimating soup profile from empty droplets...\n")
sc <- SoupChannel(tod = tod, toc = toc)

# === SOUP PROFILE DIAGNOSTICS ===
cat("\n", strrep("=", 60), "\n")
cat("SOUP PROFILE (Ambient RNA Composition)\n")
cat(strrep("=", 60), "\n\n")

cat("Top 30 genes in ambient RNA (by expression fraction):\n")
top_soup <- head(sc$soupProfile[order(-sc$soupProfile$est), ], 30)
print(top_soup)

# Save full soup profile
soup_file <- paste0(out_prefix, "_soup_profile.csv")
write.csv(sc$soupProfile, soup_file, row.names = TRUE)
cat("\nFull soup profile saved to:", soup_file, "\n")

# Highlight specific gene classes
cat("\nKey ambient RNA markers:\n")

# Hemoglobin
hb_genes <- c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2")
hb_present <- hb_genes[hb_genes %in% rownames(sc$soupProfile)]
if (length(hb_present) > 0) {
  cat("  Hemoglobin genes:\n")
  for (g in hb_present) {
    pct <- 100 * sc$soupProfile[g, "est"]
    cat(sprintf("    %s: %.2f%% of soup\n", g, pct))
  }
}

# Mitochondrial
mt_genes <- rownames(sc$soupProfile)[grepl("^MT-", rownames(sc$soupProfile))]
if (length(mt_genes) > 0) {
  total_mt <- sum(sc$soupProfile[mt_genes, "est"])
  cat(sprintf("  Total mitochondrial: %.2f%% of soup (%d genes)\n", 
              100 * total_mt, length(mt_genes)))
}

# Ribosomal
ribo_genes <- rownames(sc$soupProfile)[grepl("^RP[SL]", rownames(sc$soupProfile))]
if (length(ribo_genes) > 0) {
  total_ribo <- sum(sc$soupProfile[ribo_genes, "est"])
  cat(sprintf("  Total ribosomal: %.2f%% of soup (%d genes)\n", 
              100 * total_ribo, length(ribo_genes)))
}

# === CLUSTER INFORMATION ===
cat("\n", strrep("=", 60), "\n")
cat("CLUSTER ASSIGNMENT\n")
cat(strrep("=", 60), "\n\n")

clusters <- read.csv(file.path(in_dir, "clusters.csv"), stringsAsFactors = FALSE)
cl <- setNames(clusters$cluster, clusters$barcode)
cl <- cl[colnames(toc)]

if (any(is.na(cl))) {
  stop("ERROR: Some cells missing cluster assignments!")
}

sc <- setClusters(sc, cl)

cat("Cell type distribution:\n")
print(table(cl))

# === CONTAMINATION ESTIMATION ===
cat("\n", strrep("=", 60), "\n")
cat("ESTIMATING CONTAMINATION\n")
cat(strrep("=", 60), "\n\n")

cat("Running autoEstCont()...\n")
sc <- autoEstCont(sc)

# Extract contamination estimates
rho_vec <- sc$metaData$rho

cat("\nOverall contamination summary:\n")
print(summary(rho_vec))

cat("\nContamination by cell type:\n")
rho_by_cluster <- tapply(rho_vec, sc$metaData$clusters, function(x) {
  c(Mean = mean(x), Median = median(x), SD = sd(x), Min = min(x), Max = max(x))
})
print(do.call(rbind, rho_by_cluster))

# Save contamination estimates
rho_file <- paste0(out_prefix, "_rho_per_cell.csv")
write.csv(rho_df <- data.frame(
  barcode = colnames(toc),
  rho = rho_vec,
  cluster = sc$metaData$clusters
), rho_file, row.names = FALSE)
cat("\nContamination estimates saved to:", rho_file, "\n")

# === ADJUST COUNTS ===
cat("\n", strrep("=", 60), "\n")
cat("ADJUSTING COUNTS\n")
cat(strrep("=", 60), "\n\n")

cat("Applying correction and rounding to integers...\n")
out <- adjustCounts(sc, roundToInt = TRUE)

# Calculate total correction
total_before <- sum(toc)
total_after <- sum(out)
pct_removed <- 100 * (1 - total_after / total_before)

cat(sprintf("\nTotal reads removed: %.2f%%\n", pct_removed))
cat(sprintf("  Before: %s reads\n", format(total_before, big.mark = ",")))
cat(sprintf("  After:  %s reads\n", format(total_after, big.mark = ",")))

# === SAVE OUTPUTS ===
cat("\n", strrep("=", 60), "\n")
cat("SAVING OUTPUTS\n")
cat(strrep("=", 60), "\n\n")

# Save SoupChannel object
sc_file <- paste0(out_prefix, "_SoupX_channel.rds")
saveRDS(sc, sc_file)
cat(sprintf("SoupChannel saved to: %s\n", sc_file))

# Save corrected counts
mtx_file <- paste0(out_prefix, "_corrected_counts.mtx")
writeMM(out, mtx_file)
cat(sprintf("Corrected counts saved to: %s\n", mtx_file))

# Save gene and barcode lists
genes_file <- paste0(out_prefix, "_genes.tsv")
barcodes_file <- paste0(out_prefix, "_barcodes.tsv")
writeLines(rownames(out), genes_file)
writeLines(colnames(out), barcodes_file)
cat(sprintf("Genes saved to: %s\n", genes_file))
cat(sprintf("Barcodes saved to: %s\n", barcodes_file))

# === CLOSE LOG ===
cat("\n", strrep("=", 60), "\n")
cat("COMPLETE!\n")
cat(strrep("=", 60), "\n")
cat(sprintf("Log saved to: %s\n", log_file))

sink(type = "message")
sink(type = "output")
close(log_con)

message("SoupX correction completed successfully: ", out_prefix)