#install.packages("colorspace")
#install.packages("stringi")
#install.packages("tidytable")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("tidytable")
library("ggplot2")

pvalueFilter=0.05      #p值过滤条件
adjPvalFilter=0.05     #矫正后p值的过滤条件

#定义图形的颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
  colorSel="pvalue"
}

drugFile="DSigDB_All_detailed.txt"      #药物和基因的关系文件
hubFile="interGenes.txt"                  #下游聚焦两交叉复现 hub: GPAT3(DSigDB中作AGPAT9)/PLAU (课程原示例5基因已备份至 interGenes_course_example.txt)
setwd("D:/BaiduSyncdisk/daixie-colorectal/5.drug_AND_dock/27.drugEnrich")      #设置工作目录(已改为本项目路径)

#读取交集特征基因的列表文件
rt=read.table(hubFile, header=F, sep="\t", check.names=F)
#提取基因的名称
genes=unique(as.vector(rt[,1]))

#读取药物和基因的关系文件
drugRT=read.table(drugFile, header=T, sep="\t", check.names=F, quote="", comment.char="")
drugRT=drugRT[,1:2]

#药物富集分析
# 注意: DSigDB 以 AGPAT9 收录 GPAT3,故 interGenes.txt 须用 AGPAT9 才能捕获 GPAT3 的药物靶点
kk=suppressWarnings(enricher(genes,
            pvalueCutoff=1, qvalueCutoff=1,
            minGSSize = 1, maxGSSize = 500,   # 输入仅2个hub基因(AGPAT9/PLAU)，放宽至1以捕获靶向单基因的药物(默认10会把小靶点集药物漏掉)
            TERM2GENE=drugRT))

# qvalue(Storey)在仅2个hub、p值分布极端退化时会计算失败->用 BH p.adjust 兜底填充 qvalue 列,避免整列 NA
if (any(is.na(kk@result$qvalue))) {
  kk@result$qvalue[is.na(kk@result$qvalue)] <- kk@result$p.adjust[is.na(kk@result$qvalue)]
  message("提示: 部分/全部 qvalue 不可用(小样本 p 值分布退化),已用 BH p.adjust 兜底。过滤与出图均以 p.adjust 为准。")
}

#输出显著富集的结果
DRUG=as.data.frame(kk)
DRUG=DRUG[(DRUG$pvalue<pvalueFilter & DRUG$p.adjust<adjPvalFilter),]
write.table(DRUG[,-c(3,4)], file="DRUG.enrich.xls", sep="\t", quote=F, row.names = F)

#设置展示药物的数目
showNum=30
if(nrow(DRUG)<showNum){
  showNum=nrow(DRUG)
}

#柱状图 / 气泡图(ggplot2 直接绘制, 绕开 enrichplot 在 qvalue 退化时输出空/损坏 PDF 的坑)
# 注: 每个药物仅靶向 1 个 hub(Count=1), 故用 -log10(BH p) 表征显著性, fill/分面按 hub 基因
plotDF <- DRUG[order(DRUG$p.adjust), ]
plotDF$Drug <- factor(plotDF$ID, levels=plotDF$ID[order(-log10(plotDF$p.adjust))])
pdat <- plotDF[1:showNum, ]

pdf(file="barplot.pdf", width=11, height=9)
print(ggplot(pdat, aes(x=-log10(p.adjust), y=Drug, fill=geneID)) +
  geom_col() +
  scale_fill_manual(values=c("PLAU"="#2E9BD6","AGPAT9"="#E69F00"), name="Hub gene") +
  labs(x=expression(-log[10](BH~adjusted~p)), y="Drug", title="Drug enrichment: GPAT3(AGPAT9) & PLAU") +
  theme_bw(base_size=10) + theme(plot.title=element_text(hjust=0.5), axis.text.y=element_text(size=7)))
dev.off()

pdf(file="bubble.pdf", width=11, height=9)
print(ggplot(pdat, aes(x=geneID, y=Drug, size=-log10(p.adjust), color=-log10(p.adjust))) +
  geom_point(alpha=0.8) +
  scale_color_gradient(low="#FDDBC7", high="#2166AC", name=expression(-log[10](BH~adj~p))) +
  scale_size(range=c(2,9), name=expression(-log[10](BH~adj~p))) +
  labs(x="Hub gene", y="Drug", title="Drug enrichment: GPAT3(AGPAT9) & PLAU") +
  theme_bw(base_size=10) + theme(plot.title=element_text(hjust=0.5), axis.text.y=element_text(size=7)))
dev.off()

#输出网络关系文件
DRUG=DRUG[,c("ID", "geneID")]
DRUG2=separate_rows(DRUG, geneID, sep="/")
DRUG2=DRUG2 %>% group_by(geneID) %>% slice_head(n=20)     #设置展示药物的数目
networkTab=cbind(DRUG2, "Drug")
colnames(networkTab)=c("Node1", "Node2", "Type")
write.table(file="net.network.txt", networkTab, sep="\t", quote=F, row.names=F)
#输出节点属性文件
drugNode=cbind(unique(networkTab[,1]),"Drug")
colnames(drugNode)=c("Node", "Type")
geneNode=cbind(unique(networkTab[,2]),"Gene")
colnames(geneNode)=c("Node", "Type")
nodeTab=rbind(drugNode, geneNode)
write.table(file="net.node.txt", nodeTab, sep="\t", quote=F, row.names=F)

