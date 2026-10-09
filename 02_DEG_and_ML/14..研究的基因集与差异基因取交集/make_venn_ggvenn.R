suppressMessages(library(ggvenn))
suppressMessages(library(ggplot2))

# ===== 读取三个基因集（ASCII 临时目录，绕开中文路径坑）=====
crc_genes  <- readLines("D:/geotmp/crc_genes.txt")
mets_genes <- readLines("D:/geotmp/mets_genes.txt")
rese_genes <- readLines("D:/geotmp/rese_genes.txt")

# 去空 + 转大写，保证交集准确
clean <- function(x) sort(unique(toupper(trimws(x[x != ""]))))
crc  <- clean(crc_genes)
mets <- clean(mets_genes)
rese <- clean(rese_genes)
cat(sprintf("CRC=%d  MetS=%d  TNF/NF-kB=%d\n",
            length(crc), length(mets), length(rese)))

# ===== 构造命名列表（顺序 = 配色顺序）=====
gene_list <- list(
  "CRC"         = crc,
  "MetS"        = mets,
  "TNF/NF-kB"   = rese
)

# ===== 照 差异分析.R 的 ggvenn 风格（白描边 + 蓝橙 + 百分比）=====
p <- ggvenn(gene_list,
            show_percentage = TRUE,
            stroke_color    = "white",
            stroke_size     = 1.0,
            fill_color      = c("#1E90FF", "#FF8C00", "#984EA3"),
            set_name_color  = c("#1E90FF", "#FF8C00", "#984EA3"),
            set_name_size   = 7,
            text_size       = 5)

# ===== 输出 PNG(600dpi) + PDF =====
ggsave("D:/geotmp/venn_ggvenn.png", p, width = 8, height = 8, dpi = 600)
ggsave("D:/geotmp/venn_ggvenn.pdf", p, width = 8, height = 8)
cat("Venn saved -> D:/geotmp/venn_ggvenn.png / .pdf\n")
