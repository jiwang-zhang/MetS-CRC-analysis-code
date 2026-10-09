library(enrichR)
library(ggplot2)
library(dplyr)
options(bitmapType = "cairo")

enrich_dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                "GO_Molecular_Function_2023", "KEGG_2021_Human")

hubGenes <- c("GPAT3", "PLAU")
baseDir  <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/31.Enrich"
oeDir    <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/30.scTenifoldOE"

for (target_gene in hubGenes) {
  cat("\n========================================\n")
  cat("Enrichment for", target_gene, "OE\n")
  cat("========================================\n")
  
  geneDir <- file.path(baseDir, target_gene)
  dir.create(geneDir, showWarnings = FALSE, recursive = TRUE)
  setwd(geneDir)
  
  sigFile <- file.path(oeDir, target_gene, "sigDiff.txt")
  if (!file.exists(sigFile)) { cat("  sigDiff.txt not found, skipping\n"); next }
  rt <- read.table(sigFile, header = TRUE, sep = "\t", check.names = FALSE)
  genes <- unique(as.vector(rt[, 1]))
  cat("  Genes:", length(genes), "\n")
  # 2026-08-24: 下游聚焦 GPAT3/PLAU 两 hub; <1 阈值保留以便小扰动集也能出最小富集图
  if (length(genes) < 1) { cat("  <1 genes, skipping enrichment\n"); next }
  
  enriched <- enrichr(genes, enrich_dbs)
  all_results <- list()
  for (db in enrich_dbs) {
    df <- enriched[[db]]
    if (nrow(df) > 0) {
      df$Database <- db
      all_results[[length(all_results) + 1]] <- df
    }
  }
  if (length(all_results) == 0) { cat("  No results\n"); next }
  
  GK <- do.call(rbind, all_results)
  colnames(GK)[colnames(GK) == "Term"] <- "Description"
  colnames(GK)[colnames(GK) == "P.value"] <- "pvalue"
  write.table(GK, file = "enrichment.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  data_selected <- GK %>% arrange(Database, pvalue) %>%
    group_by(Database) %>% slice_head(n = 5) %>% ungroup()
  data_selected$log10p <- -log10(data_selected$pvalue)
  data_selected <- data_selected %>%
    arrange(Database, desc(log10p)) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  db_colors <- c("GO_Biological_Process_2023" = "#E178A8",
                 "GO_Cellular_Component_2023" = "#A0E070",
                 "GO_Molecular_Function_2023" = "#B5E1FF",
                 "KEGG_2021_Human" = "#FFB143")
  # 2026-08-25: 图例缩短为 BP/CC/MF/KEGG, 防止底部单行图例超出 PDF 宽度被截断
  db_labels <- c("GO_Biological_Process_2023" = "GO Biological Process",
                 "GO_Cellular_Component_2023" = "GO Cellular Component",
                 "GO_Molecular_Function_2023" = "GO Molecular Function",
                 "KEGG_2021_Human" = "KEGG")

  p <- ggplot(data_selected, aes(x = log10p, y = Description, fill = Database)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = Description, x = 0.05), hjust = 0, size = 3, color = "black") +
    scale_fill_manual(values = db_colors, labels = db_labels, name = "Database") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = paste0(target_gene, " OE: Enrichment"), x = "-log10(pvalue)", y = "") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), legend.position = "bottom")
  
  pdf(file = "barplot.pdf", width = 9, height = 6)
  print(p)
  dev.off()
  cat("  Plot saved: barplot.pdf\n")
}
cat("\n>>> ALL DONE: 2 genes (GPAT3/PLAU) OE enrichment (enrichR) completed.\n")
