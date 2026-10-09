suppressMessages(library(ggVennDiagram))
suppressMessages(library(ggplot2))

crc  <- readLines("D:/geotmp/crc_genes.txt")
mets <- readLines("D:/geotmp/mets_genes.txt")
rese <- readLines("D:/geotmp/rese_genes.txt")

x <- list("CRC DEGs" = crc, "MetS DEGs" = mets, "TNF/NF-kB signature" = rese)

p <- ggVennDiagram(x,
                   category.names = c("CRC DEGs", "MetS DEGs", "TNF/NF-kB signature"),
                   label = "count",
                   label.color = "black",
                   label.alpha = 0,
                   edge.lty = "dashed",
                   edge.size = 1) +
  scale_fill_gradient(low = "white", high = "#1F7A8C", name = "gene count") +
  theme(legend.position = "right")

ggsave("D:/geotmp/venn_gg.png", plot = p, width = 8, height = 8, dpi = 600)
ggsave("D:/geotmp/venn_gg.pdf", plot = p, width = 8, height = 8)
cat("OK: ggVennDiagram figure saved -> D:/geotmp/venn_gg.png | .pdf\n")
