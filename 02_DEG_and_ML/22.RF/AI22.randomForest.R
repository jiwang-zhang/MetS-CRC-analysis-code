#install.packages("randomForest")
#install.packages("ggplot2")


#引用包
set.seed(12345)
library(limma)
library(randomForest)
library(ggplot2)

expFile="merged_after_batch_removal.csv"   #表达数据文件（项目CRC合并矩阵,逗号分隔）
geneFile="region_triple_R.csv"             #36个三交集候选基因
# 设置工作目录：RStudio 中先把目录设为脚本所在文件夹再 Source；
# 若用 Rscript 命令行请加编码参数：Rscript --encoding=UTF-8 AI22.randomForest.R
try(setwd("D:/daixie-colorectal/第二部分DEG/22.RF"), silent=TRUE)

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
data=t(data)
group=sub(".*_", "", row.names(data))   #取样本名最后一段作为分组(tre/con)

#随机森林树
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=5, height=5)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#获取基因重要性的评分
importance=importance(x=rf2)

#挑选疾病的特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>2])       #挑选重要性评分大于2的基因（放宽阈值，对齐课程13.RF版）
write.table(rfGenes, file="RF.gene.txt", sep="\t", quote=F, col.names=F, row.names=F)

##绘制基因重要性的气泡图
importance=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),,drop=F]
rt=cbind(row.names(importance), as.data.frame(importance))
colnames(rt)=c("ID", "Importance")
if(nrow(rt)>30){rt=rt[1:30,]}
rt$ID=factor(rt$ID, levels=rev(rt$ID))
rt$Importance=as.numeric(rt$Importance)
p=ggplot(rt, aes(Importance, ID)) + 
    geom_segment(aes(x=0, xend=Importance, y=ID, yend=ID), color="grey", cex=2) + 
    geom_point(aes(color=Importance), cex=5)+
    scale_colour_gradient(low="#5050FFFF", high="#BB0021FF") + 
    labs(x="Importance", y="")+
    theme(axis.text.x=element_text(color="black", size=10),
          axis.text.y=element_text(color="black", size=10))+
    theme_bw()

#输出图形
pdf(file="geneImportance.pdf", width=5, height=4.5)
print(p)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


