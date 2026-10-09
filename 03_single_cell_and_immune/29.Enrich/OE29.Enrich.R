# install.packages("enrichR")
# install.packages("ggplot2")

library(enrichR)
library(ggplot2)
library(dplyr)
options(bitmapType = "cairo")

# enrichR 需要设置物种数据库
listEnrichrDbs()  # 可选：查看可用数据库
enrich_dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                "GO_Molecular_Function_2023", "KEGG_2021_Human")

hubGenes <- c("GPAT3", "PLAU")
baseDir  <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/29.Enrich"
koDir    <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/28.scTenifoldKnk"

for (target_gene in hubGenes) {
  cat("\n========================================\n")
  cat("Enrichment for", target_gene, "KO\n")
  cat("========================================\n")

  geneDir <- file.path(baseDir, target_gene)
  dir.create(geneDir, showWarnings = FALSE, recursive = TRUE)
  setwd(geneDir)

  sigFile <- file.path(koDir, target_gene, "sigDiff.txt")
  if (!file.exists(sigFile)) { cat("  sigDiff.txt not found, skipping\n"); next }
  rt <- read.table(sigFile, header = TRUE, sep = "\t", check.names = FALSE)
  genes <- unique(as.vector(rt[, 1]))
  cat("  Genes:", length(genes), "\n")
  # 2026-08-24: 下游聚焦 GPAT3/PLAU 两 hub; <1 阈值保留以便小扰动集也能出最小富集图
  if (length(genes) < 1) { cat("  <1 genes, skipping enrichment\n"); next }

  # enrichR 直接接受 gene symbols
  enriched <- enrichr(genes, enrich_dbs)

  all_results <- list()
  for (db in enrich_dbs) {
    df <- enriched[[db]]
    if (nrow(df) > 0) {
      df$Database <- db
      all_results[[length(all_results) + 1]] <- df
      cat("  ", db, ":", sum(df$Adjusted.P.value < 0.05), "sig /", nrow(df), "total\n")
    }
  }

  if (length(all_results) == 0) { cat("  No enrichment results\n"); next }

  # 合并所有结果
  GK <- do.call(rbind, all_results)
  colnames(GK)[colnames(GK) == "Term"] <- "Description"
  colnames(GK)[colnames(GK) == "P.value"] <- "pvalue"
  colnames(GK)[colnames(GK) == "Adjusted.P.value"] <- "p.adjust"
  colnames(GK)[colnames(GK) == "Overlap"] <- "GeneRatio"

  write.table(GK, file = "enrichment.txt", sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- 柱状图 (top 5 per database) ----
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

  p <- ggplot(data_selected, aes(x = log10p, y = Description, fill = Database)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = Description, x = 0.05), hjust = 0, size = 3,
              color = "black") +
    scale_fill_manual(values = db_colors) +
    labs(title = paste0(target_gene, " KO: Enrichment"),
         x = "-log10(pvalue)", y = "") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "bottom")

  pdf(file = "barplot.pdf", width = 9, height = 6)
  print(p)
  dev.off()

  cat("  Plots saved: barplot.pdf\n")
}

cat("\n>>> ALL DONE: 2 genes (GPAT3/PLAU) KO enrichment (enrichR) completed.\n")
