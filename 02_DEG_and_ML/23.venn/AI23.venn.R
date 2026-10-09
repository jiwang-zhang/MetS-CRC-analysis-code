#install.packages("ggvenn")

#引用包
library(ggvenn)

setwd("D:/daixie-colorectal/第二部分DEG/23.venn")    #设置工作目录（已改为本项目路径）

#只匹配三个 ML 方法产出的基因文件，避免误读其它 txt
files=list.files(pattern="\\.gene\\.txt$")
files=sort(files)                                    #固定顺序：LASSO / RF / SVM-RFE
files=setdiff(files, "interFeatureGenes.txt")
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(inputFile in files){
  rt=read.table(inputFile, header=F, sep="\t", check.names=F, blank.lines.skip=T)  #读取输入文件
  geneNames=as.vector(rt[,1])
  geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
  geneNames=geneNames[geneNames!=""]       #过滤空行（避免末尾空行污染交集）
  uniqGene=unique(geneNames)               #基因取unique
  setName=sub("\\.txt$","",inputFile)      #集合名=文件名去后缀（LASSO / RF / SVM-RFE）
  geneList[[setName]]=uniqGene
}

#绘制venn图（颜色统一为蓝/橙/紫，与项目其它 Venn 同款）
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList, show_percentage = T,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#1E90FF","#FF8C00","#984EA3"),
       set_name_color = c("#1E90FF","#FF8C00","#984EA3"),
       set_name_size = 6, text_size = 4.5)
dev.off()

#输出三种方法的交集特征基因
interGenes=Reduce(intersect, geneList)
write.table(file="interFeatureGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信
