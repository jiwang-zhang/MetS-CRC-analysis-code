# ============================================================
# 三集合 Venn 图 —— RStudio 复现版
# 集合: CRC差异基因 ∩ MetS差异基因 ∩ TNF/NF-κB(freq2)研究基因集
# 运行前首次安装: install.packages("VennDiagram")
# ============================================================
library(VennDiagram)

# 1) 读取三个集合(取第一列基因名, 自动转大写去重)
crc  <- unique(toupper(read.csv("D:/daixie-colorectal/第二部分DEG/09.两个疾病分别DEG/DE_significant_CRC.csv")$Gene))
mets <- unique(toupper(read.csv("D:/daixie-colorectal/第二部分DEG/09.两个疾病分别DEG/DE_significant_MetS.csv")$Gene))
rese <- unique(toupper(read.csv("D:/daixie-colorectal/第二部分DEG/11..研究的基因集与差异基因取交集/TNF_NFKB_freq2_with_counts.csv")$Gene))

cat("CRC:", length(crc), "| MetS:", length(mets), "| TNF/NF-kB:", length(rese), "\n")

# 2) 画 Venn(填充色与 Python 版一致: 青绿/橙/紫)
venn.plot <- venn.diagram(
  x = list(CRC = crc, MetS = mets, TNF_NFKB = rese),
  filename = "D:/daixie-colorectal/第二部分DEG/11..研究的基因集与差异基因取交集/three_way_venn_R.png",
  fill = c("#2ca19f", "#e8893b", "#7b5ea7"),
  alpha = 0.55,
  cat.cex = 1.1, cex = 1.3, cat.fontface = "bold",
  main = "CRC DEGs ∩ MetS DEGs ∩ TNF/NF-κB signature",
  main.cex = 1.1
)

# 3) 计算并导出三交集(核心结果)
triple <- intersect(intersect(crc, mets), rese)
cat("★ 三交集基因数 =", length(triple), "\n")
write.csv(data.frame(Gene = triple),
          "D:/daixie-colorectal/第二部分DEG/11..研究的基因集与差异基因取交集/region_triple_R.csv",
          row.names = FALSE)
