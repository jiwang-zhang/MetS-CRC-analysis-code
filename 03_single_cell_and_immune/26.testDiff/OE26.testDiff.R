#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="GSE36980.normalize.txt"      #表达数据文件
geneFile="interGenes.List.txt"        #基因列表文件
diseaseName="Alzheimer"                 #设置图形中展示疾病的名称
setwd("C:\\Users\\Administrator\\Desktop\\OE\\26.testDiff")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取交集基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=t(data[as.vector(geneRT[,1]),])

#获取样品的分组信息(对照组和疾病组)
Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
Type=ifelse(Type=="Control", "Control", diseaseName)
Type=factor(Type, levels=c("Control", diseaseName))
rt=cbind(as.data.frame(data), Type)

#将数据转换为箱线图的输入文件
data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", fill = "Type",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Type",
	     palette = c("#0088FF", "#FF5555"), width=0.6)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出图形
pdf(file="boxplot.pdf", width=5, height=4.5)
print(p1)
dev.off()
		

######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


