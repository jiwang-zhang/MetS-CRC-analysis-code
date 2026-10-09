# Replot OE enrichment barplots from existing enrichment.txt (no enrichR re-query)
# 2026-08-25: fix bottom legend clipping (short labels + 2-row legend)
library(ggplot2)
library(dplyr)
options(bitmapType = "cairo")

baseDir <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/31.Enrich"
hubGenes <- c("GPAT3", "PLAU")

db_colors <- c("GO_Biological_Process_2023" = "#E178A8",
               "GO_Cellular_Component_2023" = "#A0E070",
               "GO_Molecular_Function_2023" = "#B5E1FF",
               "KEGG_2021_Human" = "#FFB143")
db_labels <- c("GO_Biological_Process_2023" = "GO Biological Process",
               "GO_Cellular_Component_2023" = "GO Cellular Component",
               "GO_Molecular_Function_2023" = "GO Molecular Function",
               "KEGG_2021_Human" = "KEGG")

for (target_gene in hubGenes) {
  geneDir <- file.path(baseDir, target_gene)
  enrFile <- file.path(geneDir, "enrichment.txt")
  if (!file.exists(enrFile)) { cat(target_gene, ": enrichment.txt not found, skipping\n"); next }
  GK <- read.table(enrFile, header = TRUE, sep = "\t", check.names = FALSE, quote = "")

  data_selected <- GK %>% arrange(Database, pvalue) %>%
    group_by(Database) %>% slice_head(n = 5) %>% ungroup()
  data_selected$log10p <- -log10(data_selected$pvalue)
  data_selected <- data_selected %>%
    arrange(Database, desc(log10p)) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))

  p <- ggplot(data_selected, aes(x = log10p, y = Description, fill = Database)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = Description, x = 0.05), hjust = 0, size = 3, color = "black") +
    scale_fill_manual(values = db_colors, labels = db_labels, name = "Database") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = paste0(target_gene, " OE: Enrichment"), x = "-log10(pvalue)", y = "") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), legend.position = "bottom")

  setwd(geneDir)
  pdf(file = "barplot.pdf", width = 9, height = 6)
  print(p)
  dev.off()
  cat(target_gene, ": barplot.pdf re-saved\n")
}
cat("\n>>> Replot done.\n")
