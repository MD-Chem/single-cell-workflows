# Bulk RNA-seq exploratory analysis and differential expression
# Source-derived, edited excerpt; private input data are not supplied or rerun.
# Upstream tutorial cited by the supplied script:
# Love, Anders, Kim and Huber, "RNA-seq workflow: gene-level exploratory analysis
# and differential expression" (2019).
# https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# Upstream source-script credit retained: Daniel Muliaditan; Version 16/06/2022.
# Supporting bulk transcriptomics code used/adapted in the analysis;
# this file makes no sole-authorship claim.
# Editorial changes replace sample/contrast/reference lists with private inputs;
# no experimental group assignments, gene results or pathway results are bundled.
# Count filtering, thresholds and model design require review for each study.

if (!identical(Sys.getenv("RUN_PRIVATE_ANALYSIS"), "yes")) {
  stop("Code-only portfolio excerpt. Configure private inputs before setting RUN_PRIVATE_ANALYSIS=yes.")
}

library(tximeta)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(apeglm)
library(AnnotationDbi)
library(fgsea)
library(tidyverse)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(MatrixGenerics::rowVars)

input_root <- "data/private"
out_dir <- "outputs/private/bulk-rna"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Manifest columns: names, experiment, condition. No manifest is distributed.
coldata <- read.csv(file.path(input_root, "bulk_sample_manifest.csv"),
                     stringsAsFactors = FALSE)
stopifnot(all(c("names", "experiment", "condition") %in% names(coldata)))
stopifnot(!anyDuplicated(coldata$names))
coldata$files <- file.path(input_root, "quantification", coldata$names, "quant.sf")
stopifnot(all(file.exists(coldata$files)))
coldata$experiment <- factor(coldata$experiment)
coldata$condition <- factor(coldata$condition)

# The Salmon reference must be compatible with tximeta and gene summarization.
se <- tximeta(coldata)
gse <- summarizeToGene(se)
saveRDS(gse, file.path(out_dir, "gene_counts.rds"))
saveRDS(coldata, file.path(out_dir, "sample_metadata.rds"))

# Supply an AnnotationDb matching the gene identifiers used for quantification.
# The species-specific database and experiment/contrast assignments are omitted.
annotation_db <- AnnotationDbi::loadDb(file.path(input_root, "annotation.sqlite"))
analysis_plan <- readRDS(file.path(input_root, "bulk_analysis_plan.rds"))
# Plan is a list of records with exp_name, conditions_order, ref_condition.
# Optional local GMT must use the same identifiers as the ranked gene list.
gmt_file <- file.path(input_root, "gene_sets.gmt")
if (!file.exists(gmt_file)) gmt_file <- NULL

run_experiment <- function(gse, exp_name, conditions_order,
                           ref_condition,
                           gsea_gmt = NULL, annotation_db) {
  
  cat("\n============================================\n")
  cat("EXPERIMENT:", exp_name, "\n")
  cat("============================================\n")
  
  # ── Subset to this experiment ───────────────────────────
  keep_idx <- gse$experiment == exp_name
  sub_gse  <- gse[, keep_idx]
  sub_gse$condition <- factor(sub_gse$condition,
                              levels = conditions_order)
  sub_gse$condition <- relevel(sub_gse$condition, ref = ref_condition)
  
  # ── Build DESeq2 object ──────────────────────────────────
  dds <- DESeqDataSet(sub_gse, design = ~ condition)
  
  # Filter low counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  cat("Genes after filtering:", nrow(dds), "\n")
  
  # ── Transformations ──────────────────────────────────────
  vsd <- vst(dds,  blind = FALSE)
  rld <- rlog(dds, blind = FALSE)
  
  # ── Sample distances heatmap ─────────────────────────────
  sampleDists      <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$names, sep = " - ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  pdf(file.path(out_dir, paste0(exp_name, "_sample_distances.pdf")),
      width = 8, height = 7)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col   = colors,
           main  = paste(exp_name, "— Sample Euclidean Distances"))
  dev.off()
  
  # ── Poisson distances ────────────────────────────────────
  poisd            <- PoissonDistance(t(counts(dds)))
  samplePoisMatrix <- as.matrix(poisd$dd)
  rownames(samplePoisMatrix) <- paste(dds$condition, dds$names, sep = " - ")
  colnames(samplePoisMatrix) <- NULL
  
  pdf(file.path(out_dir, paste0(exp_name, "_poisson_distances.pdf")),
      width = 8, height = 7)
  pheatmap(samplePoisMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col  = colors,
           main = paste(exp_name, "— Poisson Distances"))
  dev.off()
  
  # ── PCA plot ─────────────────────────────────────────────
  pcaData    <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pca_plot <- ggplot(pcaData,
                     aes(x = PC1, y = PC2,
                         color = condition, label = name)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_bw() +
    ggtitle(paste("PCA —", exp_name))
  
  ggsave(file.path(out_dir, paste0(exp_name, "_PCA.pdf")),
         pca_plot, width = 8, height = 6)
  
  # ── GLM-PCA (fix: increased maxIter) ─────────────────────
  gpca     <- glmpca(counts(dds), L = 2,
                     ctl = list(maxIter = 5000, tol = 1e-4))
  gpca.dat <- as.data.frame(gpca$factors)
  colnames(gpca.dat) <- c("dim1", "dim2")
  gpca.dat$condition <- dds$condition
  
  gpca_plot <- ggplot(gpca.dat,
                      aes(x = dim1, y = dim2, color = condition)) +
    geom_point(size = 4) +
    coord_fixed() +
    ggtitle(paste("GLM-PCA —", exp_name)) +
    theme_bw() +
    theme(axis.text       = element_text(size = 14),
          axis.title      = element_text(size = 16),
          legend.text     = element_text(size = 12),
          legend.title    = element_text(size = 14),
          plot.title      = element_text(size = 16))
  
  ggsave(file.path(out_dir, paste0(exp_name, "_GLMPCA.pdf")),
         gpca_plot, width = 8, height = 6)
  
  # ── Z-score calculation ──────────────────────────────────
  dds  <- estimateSizeFactors(dds)
  vsd2 <- assay(vst(dds))
  Z    <- t(scale(t(vsd2)))
  Z    <- Z[complete.cases(Z), ]
  
  # Annotate rows with gene symbols
  ens.str     <- sub("\\.[0-9]+$", "", rownames(Z))
  sym         <- mapIds(annotation_db,
                        keys      = ens.str,
                        column    = "SYMBOL",
                        keytype   = "ENSEMBL",
                        multiVals = "first")
  rownames(Z) <- ifelse(is.na(sym), ens.str, sym)
  
  # Average Z-score per condition
  conditions_present <- levels(droplevels(sub_gse$condition))
  Z_mean <- sapply(conditions_present, function(cond) {
    cols <- coldata$names[coldata$experiment == exp_name &
                            coldata$condition  == cond]
    rowMeans(Z[, cols, drop = FALSE])
  })
  
  write.table(Z_mean,
              file.path(out_dir, paste0(exp_name, "_Zscore_mean.txt")),
              quote     = FALSE,
              sep       = "\t",
              row.names = TRUE,
              col.names = TRUE)
  
  # ── Run DESeq2 ───────────────────────────────────────────
  dds <- DESeq(dds)
  cat("Contrasts available:\n"); print(resultsNames(dds))
  
  # ── All contrasts vs the configured reference ─────────────────────────────
  non_ref <- conditions_order[conditions_order != ref_condition]
  non_ref <- non_ref[non_ref %in% levels(droplevels(sub_gse$condition))]
  
  all_results <- list()
  
  for (cond in non_ref) {
    
    contrast_name <- paste0(exp_name, "_", cond, "_vs_", ref_condition)
    cat("\nRunning contrast:", contrast_name, "\n")
    
    res <- results(dds,
                   contrast = c("condition", cond, ref_condition),
                   alpha    = 0.05)
    
    # ── LFC shrinkage for MA plot ──────────────────────────
    coef_name  <- paste0("condition_", cond, "_vs_", ref_condition)
    res_shrunk <- tryCatch(
      lfcShrink(dds, coef = coef_name, type = "apeglm"),
      error = function(e) {
        cat("  lfcShrink skipped:", e$message, "\n")
        res
      }
    )
    
    # MA plot
    pdf(file.path(out_dir, paste0(contrast_name, "_MAplot.pdf")),
        width = 7, height = 5)
    plotMA(res_shrunk, ylim = c(-10, 10),
           main = paste("MA —", contrast_name))
    dev.off()
    
    # ── Annotate with gene symbols ─────────────────────────
    ens.str    <- sub("\\.[0-9]+$", "", rownames(res))
    res$symbol <- mapIds(annotation_db,
                         keys      = ens.str,
                         column    = "SYMBOL",
                         keytype   = "ENSEMBL",
                         multiVals = "first")
    res$entrez <- mapIds(annotation_db,
                         keys      = ens.str,
                         column    = "ENTREZID",
                         keytype   = "ENSEMBL",
                         multiVals = "first")
    
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      arrange(padj)
    
    # Save all genes
    write.csv(res_df,
              file.path(out_dir, paste0(contrast_name, "_all_genes.csv")),
              row.names = FALSE)
    
    # Save significant genes only
    sig <- res_df %>% filter(padj <= 0.05, abs(log2FoldChange) > 1)
    write.csv(sig,
              file.path(out_dir, paste0(contrast_name, "_significant.csv")),
              row.names = FALSE)
    
    cat("  Total genes tested:                ", nrow(res_df), "\n")
    cat("  Significant (padj<=0.05, |LFC|>1): ", nrow(sig), "\n")
    cat("  Upregulated:                       ", sum(sig$log2FoldChange > 0, na.rm = TRUE), "\n")
    cat("  Downregulated:                     ", sum(sig$log2FoldChange < 0, na.rm = TRUE), "\n")
    
    # ── Volcano plot ───────────────────────────────────────
    top_labels <- res_df %>%
      filter(!is.na(padj), !is.na(symbol)) %>%
      slice_min(padj, n = 15)
    
    vol_plot <- res_df %>%
      filter(!is.na(padj)) %>%
      mutate(sig_group = case_when(
        padj <= 0.05 & log2FoldChange >  1 ~ "Up",
        padj <= 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "NS"
      )) %>%
      ggplot(aes(log2FoldChange, -log10(pvalue), color = sig_group)) +
      geom_point(alpha = 0.4, size = 0.8) +
      scale_color_manual(
        values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")
      ) +
      ggrepel::geom_text_repel(
        data        = top_labels,
        aes(label   = symbol),
        size        = 2.5,
        color       = "black",
        max.overlaps = 20
      ) +
      geom_vline(xintercept = c(-1, 1),
                 linetype = "dashed", linewidth = 0.4) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", linewidth = 0.4) +
      theme_bw() +
      labs(title    = contrast_name,
           subtitle = paste0("Up: ", sum(sig$log2FoldChange > 0, na.rm=TRUE),
                             "  |  Down: ", sum(sig$log2FoldChange < 0, na.rm=TRUE),
                             "  |  Total sig: ", nrow(sig)),
           x = "log2 Fold Change",
           y = "-log10(p-value)",
           color = NULL)
    
    ggsave(file.path(out_dir, paste0(contrast_name, "_volcano.pdf")),
           vol_plot, width = 8, height = 6)
    
    all_results[[cond]] <- res_df
    
    # ── GSEA ──────────────────────────────────────────────
    if (!is.null(gsea_gmt) && file.exists(gsea_gmt)) {
      
      pathways <- gmtPathways(gsea_gmt)
      
      res_gsea <- res_df %>%
        dplyr::select(symbol, stat) %>%
        filter(!is.na(symbol), !is.na(stat)) %>%
        distinct() %>%
        group_by(symbol) %>%
        summarize(stat = mean(stat), .groups = "drop")
      
      ranks <- deframe(res_gsea)
      
      # Fix: use fgseaMultilevel (no nperm)
      fgseaRes <- fgsea(pathways = pathways,
                        stats    = ranks,
                        minSize  = 15,
                        maxSize  = 500)
      
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))
      
      # Save GSEA table
      write.csv(
        fgseaResTidy %>% dplyr::select(-leadingEdge),
        file.path(out_dir, paste0(contrast_name, "_GSEA.csv")),
        row.names = FALSE
      )
      
      cat("  GSEA significant pathways (padj<=0.05):",
          sum(fgseaResTidy$padj <= 0.05, na.rm = TRUE), "\n")
      
      # ── GSEA plot (fixed) ────────────────────────────────
      # ── Top 10 up + top 10 down pathways ─────────────────────────
      top_up <- fgseaResTidy %>%
        filter(!is.na(padj), !is.na(NES), NES > 0) %>%
        arrange(padj) %>%
        head(10) %>%
        mutate(direction = "Upregulated")
      
      top_down <- fgseaResTidy %>%
        filter(!is.na(padj), !is.na(NES), NES < 0) %>%
        arrange(padj) %>%
        head(10) %>%
        mutate(direction = "Downregulated")
      
      gsea_plot_data <- bind_rows(top_up, top_down) %>%
        mutate(
          pathway_short = pathway,
          pathway_short = gsub("_", " ", pathway_short),
          pathway_short = str_wrap(pathway_short, width = 35),
          direction     = factor(direction,
                                 levels = c("Upregulated", "Downregulated")),
          sig_label = paste0("p=", formatC(padj, format = "e", digits = 1))
        ) %>%
        arrange(NES)   # order bars from most down to most up
      
      if (nrow(gsea_plot_data) > 0) {
        
        # Count sig up and down for subtitle
        n_sig_up   <- sum(top_up$padj   <= 0.05, na.rm = TRUE)
        n_sig_down <- sum(top_down$padj <= 0.05, na.rm = TRUE)
        
        gsea_plot <- ggplot(
          gsea_plot_data,
          aes(x    = reorder(pathway_short, NES),
              y    = NES,
              fill = direction)) +
          geom_col(width = 0.75) +
          geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30") +
          scale_fill_manual(
            values = c("Upregulated"   = "#d73027",
                       "Downregulated" = "#4575b4")
          ) +
          # # Significance labels on bars
          geom_text(
            aes(label = sig_label,
                hjust = ifelse(NES > 0, -0.1, 1.1)),
            size  = 2.8,
            color = "black"
          ) +
          # Dividing line between up and down sections
          geom_hline(yintercept = 0, linewidth = 1,
                     color = "grey20", linetype = "solid") +
          coord_flip() +
          labs(
            x        = NULL,
            y        = "Normalized Enrichment Score (NES)",
            fill     = NULL,
            title    = paste("GSEA —", contrast_name),
            subtitle = paste0(
              "Top 10 upregulated + top 10 downregulated pathways",
              "  |  Sig up: ", n_sig_up,
              "  |  Sig down: ", n_sig_down
            )
          ) +
          # Facet to visually separate up and down
          facet_grid(direction ~ .,
                     scales = "free_y",
                     space  = "free_y") +
          theme_bw() +
          theme(
            axis.text.y     = element_text(size = 9),
            axis.text.x     = element_text(size = 9),
            plot.title      = element_text(size = 12, face = "bold"),
            plot.subtitle   = element_text(size = 9),
            legend.position = "none",        # direction shown by facet label
            strip.text      = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "grey92", color = "grey70"),
            panel.spacing   = unit(0.8, "lines"),
            plot.margin     = margin(10, 55, 10, 10)
          ) +
          expand_limits(
            y = c(min(gsea_plot_data$NES, na.rm = TRUE) * 1.3,
                  max(gsea_plot_data$NES, na.rm = TRUE) * 1.3)
          )
        
        ggsave(
          file.path(out_dir, paste0(contrast_name, "_GSEA.pdf")),
          gsea_plot,
          width  = 11,
          height = max(7, nrow(gsea_plot_data) * 0.42),
          units  = "in"
        )
        
        cat("  GSEA plot saved — Up:", nrow(top_up),
            "| Down:", nrow(top_down), "pathways\n")
        
      } else {
        cat("  No GSEA pathways to plot for", contrast_name, "\n")
      }
    }
  }
  
  # ── Top variable genes heatmap (fix: MatrixGenerics::rowVars) ──
  topVarGenes <- head(
    order(MatrixGenerics::rowVars(assay(rld)), decreasing = TRUE),
    50
  )
  mat <- assay(rld)[topVarGenes, ]
  mat <- mat - rowMeans(mat)
  
  ens.str       <- sub("\\.[0-9]+$", "", rownames(mat))
  sym           <- mapIds(annotation_db,
                          keys      = ens.str,
                          column    = "SYMBOL",
                          keytype   = "ENSEMBL",
                          multiVals = "first")
  rownames(mat) <- ifelse(is.na(sym), ens.str, sym)
  
  anno <- as.data.frame(colData(rld)[, "condition", drop = FALSE])
  
  pdf(file.path(out_dir, paste0(exp_name, "_topvar_heatmap.pdf")),
      width = 9, height = 11)
  pheatmap(mat,
           annotation_col = anno,
           show_rownames  = TRUE,
           fontsize_row   = 7,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           main  = paste("Top 50 Variable Genes —", exp_name))
  dev.off()
  
  cat("\nAll outputs saved for", exp_name, "\n")
  return(list(dds = dds, rld = rld, vsd = vsd, results = all_results))
}

# Run the same experiment-level function for each privately configured analysis.
analysis_results <- lapply(analysis_plan, function(plan) {
  run_experiment(gse, exp_name = plan$exp_name,
                 conditions_order = plan$conditions_order,
                 ref_condition = plan$ref_condition,
                 gsea_gmt = gmt_file, annotation_db = annotation_db)
})
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
