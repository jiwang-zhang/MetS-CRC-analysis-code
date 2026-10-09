#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#引用包
# 注:若用VCF格式结局(如gwasglue的readVcf)才需下面两个包;
# 当前CRC结局为制表符GWAS,仅需 TwoSampleMR
# library(VariantAnnotation)
# library(gwasglue)
library(TwoSampleMR)

exposureFile="(p5e-08)exposure_data.csv"       #暴露:eQTL基因表达(全eQTLGen基因;可筛选为7 hub)
outcomeFile="finngen_R13_C3_COLORECTAL.gz"     #结局:FinnGen CRC GWAS(制表符分隔,非VCF)
# 备选结局(基因→MetS方向): "vanWalree_2022/CF_MetS_results_adapt.txt.gz" (需另写读取映射)
outcomeName="Colorectal cancer"        #设置图形中展示疾病的名称
setwd("D:/daixie-colorectal/3.MR/16.MR")     #设置工作目录

#读取暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

#读取结局数据:制表符GWAS汇总(直接读,无需VCF)
#FinnGen列: chrom pos ref alt rsids ... pval mlogp beta sebeta af_alt ...
#beta以alt为效应allele -> effect_allele=alt, other_allele=ref, eaf=af_alt
outcomeData=read_outcome_data(snps=exposure_dat$SNP,
                              filename=outcomeFile, sep="\t",
                              snp_col="rsids",
                              beta_col="beta",
                              se_col="sebeta",
                              effect_allele_col="alt",
                              other_allele_col="ref",
                              pval_col="pval",
                              eaf_col="af_alt",
                              chr_col="chrom", pos_col="pos")
# 注: 老版本 TwoSampleMR 的 read_outcome_data 不认 outcome/clump 参数,
#     故此处只传数据映射; 结局名由下方 outcomeData$outcome=outcomeName 补设,
#     clump 在老版本默认不执行(如需发表级clump需另跑 clump_data())

#将暴露数据和结局数据合并
outcomeData$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcomeData)
dat=dat[dat$pval.outcome>0.05,]

#输出用于孟德尔随机化的工具变量
dat=dat[dat$mr_keep=="TRUE",]
write.csv(dat, file="table.SNP.csv", row.names=F)

#孟德尔随机化分析
mrResult=mr(dat)

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


