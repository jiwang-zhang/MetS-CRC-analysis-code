#install.packages("e1071")


#引用包
set.seed(12345)
library(limma)
library(e1071)

expFile="merged_after_batch_removal.csv"   #表达数据文件（项目CRC合并矩阵,逗号分隔）
geneFile="region_triple_R.csv"             #36个三交集候选基因
# 设置工作目录：RStudio 中先把目录设为脚本所在文件夹再 Source；
# 若用 Rscript 命令行请加编码参数：Rscript --encoding=UTF-8 AI21.SVM.R
try(setwd("D:/daixie-colorectal/第二部分DEG/21.SVM"), silent=TRUE)
source("AI21.msvmRFE.R")

#读取表达数据文件
rt=read.table(expFile, header=T, sep=",", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件,提取核心基因的表达量
geneRT=read.table(geneFile, header=T, sep=",", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]
data=as.data.frame(t(data))

#获取样品的分组信息(对照组和疾病组)
group=sub(".*_", "", row.names(data))   #取样本名最后一段(tre/con)
data=cbind(group, data)
data$group=factor(data$group, levels=c("tre","con"))

#构建机器学习-支持向量机递归特征消除算法(SVM-RFE)
svmRFE(data, k=10, halve.above=50)
nfold=10
geneNum=nrow(data)
folds=rep(1:nfold, len=geneNum)[sample(geneNum)]
folds=lapply(1:nfold, function(x) which(folds == x))
results=lapply(folds, svmRFE.wrap, data, k=10, halve.above=50)

#根据基因的重要性对基因进行排序
top.features=WriteFeatures(results, data, save=F)
#输出基因重要性排序的结果
write.table(top.features, file="feature_svm.txt", sep="\t", quote=F,row.names=F)

#交叉验证
num=ifelse((ncol(data)-1)>30, 30, (ncol(data)-1))
featsweep=lapply(1:num, FeatSweep.wrap, results, data)

#获取交叉验证的误差
no.info=min(prop.table(table(data[,1])))
errors=sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#绘制交叉验证误差的图形
pdf(file="errors.pdf", width=5, height=5)
PlotErrors(errors, no.info=no.info)
dev.off()

#绘制交叉验证准确性的图形
pdf(file="accuracy.pdf", width=5, height=5)
Plotaccuracy(1-errors, no.info=no.info)
dev.off()

#输出SVM机器学习方法得到的疾病特征基因
featureGenes=top.features[1:which.min(errors),1,drop=F]
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信