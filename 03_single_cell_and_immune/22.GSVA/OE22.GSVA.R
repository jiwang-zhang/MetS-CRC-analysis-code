# ============================================================
# Step 7d: GSVA ssGSEA — 7 hub 基因高/低表达组间 KEGG 通路差异
# 输入: 合并表达矩阵 + hubGenes.List.txt + KEGG基因集
# 输出: 每个基因一张 barplot.pdf (通路 t 值柱状图)
#
# 依赖: limma, GSEABase, GSVA, reshape2, ggpubr
# ============================================================

library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile   <- "../../第二部分DEG/25.diagnostic_model/merged_after_batch_removal.csv"  # 表达数据
geneFile  <- "hubGenes.List.txt"          # 7 hub 基因列表
gmtFile   <- "c2.cp.kegg.Hs.symbols.gmt" # KEGG 基因集
setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/22.GSVA")

# ============================================================
# 1. 读取表达矩阵
# ============================================================
rt <- read.table(expFile, header = TRUE, sep = ",", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dn  <- list(rownames(exp), colnames(exp))
exp <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dn)
exp <- avereps(exp)

# 仅肿瘤样本 (列名含 _tre)
treCols <- grepl("_tre$", colnames(exp))
data <- exp[, treCols, drop = FALSE]
cat("表达矩阵 (仅肿瘤):", nrow(data), "基因 x", ncol(data), "样本\n")

# ============================================================
# 2. 读取基因集, 运行 ssGSEA
# ============================================================
geneSets <- getGmt(gmtFile, geneIdType = SymbolIdentifier())
ssgseaPar <- ssgseaParam(data, geneSets, minSize = 5, maxSize = 500)
gsvaScore <- gsva(ssgseaPar)

# 归一化到 [0, 1]
normalize <- function(x) { return((x - min(x)) / (max(x) - min(x))) }
gsvaScore <- normalize(gsvaScore)
gsvaScore <- gsvaScore[apply(gsvaScore, 1, sd) > 0.01, ]
cat("ssGSEA 完成:", nrow(gsvaScore), "通路\n")

# 保存
write.csv(gsvaScore, "gsva_score.csv", quote = FALSE)
cat("gsva_score.csv 已保存\n")

# ============================================================
# 3. 按 hub 基因中位数分高低组, 通路差异 t-test → barplot
# ============================================================
geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
hubGenes <- as.vector(geneRT[, 1])

for (gene in hubGenes) {
  if (!gene %in% rownames(data)) {
    cat("跳过", gene, ": 不在表达矩阵中\n")
    next
  }
  # 分组: 中位数
  geneExp <- data[gene, ]
  lowName  <- colnames(data)[geneExp < median(geneExp)]
  highName <- colnames(data)[geneExp >= median(geneExp)]
  
  lowScore  <- gsvaScore[, lowName,  drop = FALSE]
  highScore <- gsvaScore[, highName, drop = FALSE]
  scoreData <- cbind(lowScore, highScore)
  conNum    <- ncol(lowScore)
  treatNum  <- ncol(highScore)
  Type <- c(rep("Low", conNum), rep("High", treatNum))
  
  # 逐通路 t-test
  outTab <- data.frame()
  for (i in rownames(scoreData)) {
    test <- t.test(scoreData[i, ] ~ Type)
    t <- test$statistic
    pvalue <- test$p.value
    if (test$estimate[2] > test$estimate[1]) t <- abs(t) else t <- -abs(t)
    Sig <- ifelse(pvalue > 0.05, "Not", ifelse(t > 0, "Up", "Down"))
    outTab <- rbind(outTab, cbind(Pathway = i, t, pvalue, Sig))
  }
  
  notSigTab <- outTab[outTab$Sig == "Not", ]
  notSigTab <- notSigTab[order(as.numeric(notSigTab$t)), ]
  sigTab <- outTab[outTab$Sig != "Not", ]
  sigTab <- sigTab[order(as.numeric(sigTab$t)), ]
  
  if (nrow(sigTab) > 20) {
    outTab <- rbind(sigTab[c(1:10, ((nrow(sigTab) - 9):nrow(sigTab))), ],
                    notSigTab[c(1:5, ((nrow(notSigTab) - 4):nrow(notSigTab))), ])
  } else {
    outTab <- rbind(sigTab,
                    notSigTab[c(1:5, ((nrow(notSigTab) - 4):nrow(notSigTab))), ])
  }
  
  # 柱状图
  pdf(file = paste0(gene, ".barplot.pdf"), width = 10, height = 7)
  outTab$t <- as.numeric(outTab$t)
  outTab$Sig <- factor(outTab$Sig, levels = c("Down", "Not", "Up"))
  gg1 <- ggbarplot(outTab, x = "Pathway", y = "t", fill = "Sig", color = "white",
                   palette = c("royalblue1", "gray80", "tomato"),
                   sort.val = "asc", sort.by.groups = TRUE,
                   rotate = TRUE, title = gene, legend.title = "Group",
                   legend = "right", xlab = "", ylab = "t value of GSVA score",
                   x.text.angle = 60)
  print(gg1)
  dev.off()
  cat(gene, "barplot 已保存\n")
}

cat("\n=== GSVA 全部完成\n")
