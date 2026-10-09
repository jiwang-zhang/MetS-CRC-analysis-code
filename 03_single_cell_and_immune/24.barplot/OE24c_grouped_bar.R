# OE24c_grouped_bar.R
# Purpose: alternative to OE24b_mean_stacked.R -- grouped horizontal bar chart
#          showing Normal vs CRC side-by-side for each of the 22 CIBERSORT cell types.
# Input  : ../23.CIBERSORT/CIBERSORT-Results.txt  (CIBERSORT LM22, 282 samples)
# Output : mean_grouped_bar.pdf
#
# Advantage over stacked bar: every cell type is directly comparable between groups,
# making increases/decreases obvious (e.g., B naive down, CD8 up in CRC).
# Disadvantage: loses the "fractions sum to 1" composition intuition.

library(reshape2)
library(ggplot2)

# ---- read CIBERSORT results ----
rt <- read.table("../23.CIBERSORT/CIBERSORT-Results.txt",
                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
cells <- colnames(rt)[1:22]

# ---- group and compute per-group means ----
is_tre <- grepl("_tre$", rownames(rt))
df <- data.frame(
  Type = ifelse(is_tre, "CRC", "Normal"),
  rt[, cells, drop = FALSE]
)
mean_df   <- aggregate(. ~ Type, data = df, FUN = mean)
mean_long <- melt(mean_df, id.vars = "Type",
                  variable.name = "Cell", value.name = "Fraction")
mean_long$Type <- factor(mean_long$Type, levels = c("Normal", "CRC"))

# ---- plot: horizontal grouped bars ----
p <- ggplot(mean_long, aes(x = Fraction, y = reorder(Cell, Fraction), fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7,
           colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("Normal" = "#6699FF", "CRC" = "#E6550D")) +
  labs(x = "Mean CIBERSORT fraction", y = NULL,
       title = "Mean CIBERSORT fraction by cell type: Normal vs CRC",
       caption = "CIBERSORT LM22; Normal n = 5, CRC n = 276.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom right",
        plot.caption = element_text(size = 9, hjust = 0.5))

ggsave("mean_grouped_bar.pdf", p, width = 7, height = 9)
message("Saved mean_grouped_bar.pdf")
