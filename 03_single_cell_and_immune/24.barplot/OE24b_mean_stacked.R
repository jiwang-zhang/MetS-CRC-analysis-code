# OE24b_mean_stacked.R
# Purpose: replace the dense per-sample CIBERSORT stacked bar (barplot.pdf, 276 tumour bars)
#          with a clean two-bar MEAN composition plot (Normal vs CRC).
# Input  : ../23.CIBERSORT/CIBERSORT-Results.txt  (CIBERSORT LM22 deconvolution, 282 samples)
# Output : mean_stacked_bar.pdf
#
# NOTE on sample imbalance: CIBERSORT was run on a tumour-enriched subset with
# sparse normal colon samples (n = 5 normal vs n = 276 CRC). Per-group means are
# still comparable, but report this imbalance in the caption / methods.

library(reshape2)
library(ggplot2)

# ---- read CIBERSORT results (first column = sample ID used as row name) ----
rt <- read.table("../23.CIBERSORT/CIBERSORT-Results.txt",
                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# ---- 22 immune-cell fractions are the first 22 columns ----
cells <- colnames(rt)[1:22]

# ---- group samples by tumour/normal (tumour IDs end with "_tre") ----
is_tre <- grepl("_tre$", rownames(rt))
df <- data.frame(
  Type = ifelse(is_tre, "CRC", "Normal"),
  rt[, cells, drop = FALSE]
)

# ---- per-group mean fraction ----
mean_df   <- aggregate(. ~ Type, data = df, FUN = mean)
mean_long <- melt(mean_df, id.vars = "Type",
                  variable.name = "Cell", value.name = "Fraction")
mean_long$Type <- factor(mean_long$Type, levels = c("Normal", "CRC"))

# ---- 22 distinct, perceptually separated colours (HCL) ----
# 旧版 rainbow(s=0.7, v=0.7) 太浅 + 小图会把 22 行图例截断成 3 行，
# 这里改用高饱和 HCL 调色板并把图例放底部 + 增大画布。
cols <- hcl(h = seq(0, 360, length.out = 23)[1:22], c = 90, l = 55)
names(cols) <- cells

# ---- plot ----
p <- ggplot(mean_long, aes(x = Type, y = Fraction, fill = Cell)) +
  geom_bar(stat = "identity", width = 0.6, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = NULL, y = "Mean CIBERSORT fraction",
       title = "Mean immune composition: Normal vs CRC",
       caption = "CIBERSORT LM22; Normal n = 5, CRC n = 276.") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        plot.caption = element_text(size = 9, hjust = 0.5)) +
  guides(fill = guide_legend(ncol = 4, title.position = "top"))

# 增大 height 给底部图例留够空间，避免再被截断
ggsave("mean_stacked_bar.pdf", p, width = 8, height = 7)
message("Saved mean_stacked_bar.pdf")
