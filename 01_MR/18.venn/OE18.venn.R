#install.packages("ggvenn")
#install.packages("ggplot2")


library(ggvenn)
library(ggplot2)

setwd("C:\\Users\\Administrator\\Desktop\\OE\\18.venn")     #设置工作目录
upList=list()        #上调基因的列表
downList=list()      #下调基因的列表

#读取差异分析的结果文件
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
upRT=rt[rt$logFC>0,]          #提取上调基因的数据
downRT=rt[rt$logFC<0,]        #提取下调基因的数据
upGenes=as.vector(upRT[,1])         #提取上调基因名称
downGenes=as.vector(downRT[,1])     #提取下调基因名称
upList[["DEG_up"]]=upGenes
downList[["DEG_down"]]=downGenes

#读取孟德尔随机化分析结果文件
rt=read.csv("IVW.filter.csv", header=T, sep=",", check.names=F)
rt=rt[rt$method=="Inverse variance weighted",]
upRT=rt[rt$or>1,]          #提取高风险基因的数据
downRT=rt[rt$or<1,]        #提取低风险基因的数据
upGenes=unique(upRT[,"exposure"])         #提取高风险基因名称
downGenes=unique(downRT[,"exposure"])     #提取低风险基因名称
upList[["MR_or>1"]]=upGenes
downList[["MR_or<1"]]=downGenes

#==================== 绘制上调venn图 ====================
pdf(file="up.pdf", width=5, height=5)
ggvenn(upList,
       show_percentage = T,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#E41A1C","#1E90FF"),
       set_name_color =c("#E41A1C","#1E90FF"),
       set_name_size=6, text_size=4.5)
dev.off()

#==================== 绘制下调venn图 ====================
pdf(file="down.pdf", width=5, height=5)
ggvenn(downList,
       show_percentage = T,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#E41A1C","#1E90FF"),
       set_name_color =c("#E41A1C","#1E90FF"),
       set_name_size=6, text_size=4.5)
dev.off()

#输出交集基因的文件
# 上调交集基因
upInterGenes = Reduce(intersect, upList)
if(length(upInterGenes) > 0){
  upTab = cbind(Gene = upInterGenes, Type = "up")
} else {
  upTab = data.frame(Gene = character(), Type = character(), stringsAsFactors = FALSE)
}
# 下调交集基因
downInterGenes = Reduce(intersect, downList)
if(length(downInterGenes) > 0){
  downTab = cbind(Gene = downInterGenes, Type = "down")
} else {
  downTab = data.frame(Gene = character(), Type = character(), stringsAsFactors = FALSE)
}
interTab=rbind(upTab, downTab)
colnames(interTab)=c("Gene", "Type")
#输出交集基因的属性文件
write.table(file="interGenes.Type.txt", interTab, sep="\t", quote=F, col.names=T, row.names=F)
#输出交集基因的列表文件
write.table(file="interGenes.List.txt", interTab[,1], sep="\t", quote=F, col.names=F, row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: seqBio


