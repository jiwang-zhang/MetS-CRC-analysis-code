#install.packages("glmnet")


#引用包
library(limma)
library(glmnet)
set.seed(12345)

expFile="merged_after_batch_removal.csv"      #表达数据文件（项目CRC合并矩阵,逗号分隔）
geneFile="region_triple_R.csv"                 #36个三交集候选基因
setwd("D:\\daixie-colorectal\\第二部分DEG\\20.lasso")      #设置工作目录（改成你本机路径）

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
rt=t(data)

#构建lasso回归模型
x=as.matrix(rt)
y=sub(".*_", "", row.names(rt))   #取样本名最后一段作为分组(tre/con)
fit=glmnet(x, y, family = "binomial", alpha=1)
#绘制Lasso回归的图形
pdf(file="lasso.pdf", width=6, height=5.5)
plot(fit)
dev.off()
#绘制交叉验证的图形
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf", width=6, height=5.5)
plot(cvfit)
dev.off()

#输出筛选疾病的特征基因
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信