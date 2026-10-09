# ============================================================
# CRC x MetS 双疾病分别DEG → 取交集co-DEGs（按模板忠实复刻）
# CRC: merged_after_batch_removal.csv (GPL570 三合一去批次, 748样本)
# MetS: GSE98895.csv (GPL570芯片, 40样本)
# ============================================================
setwd("D:/daixie-colorectal/第二部分DEG/09.两个疾病分别DEG")

threshold_logFC    <- 0.5
threshold_adjP     <- 0.05
max_display_genes  <- 50

# ==================== Step 1: 加载R包 ====================
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(grid)
suppressMessages({library(ComplexHeatmap)})
library(circlize)

# ==================== 通用函数：单疾病完整DEG流程 ====================
run_one_disease <- function(file_path, disease_label) {
  cat("\n", paste(rep("=",60),collapse=""), "\n")
  cat("  ", disease_label, " DEG 分析\n")
  cat(paste(rep("=",60),collapse=""), "\n")
  
  # --- 读取表达矩阵 ---
  tmp_head <- readLines(file_path, 1)
  sep <- ifelse(grepl(",", tmp_head), ",", "\t")
  expr_raw <- read.table(file_path, header=TRUE, sep=sep, check.names=FALSE, stringsAsFactors=FALSE)
  rownames(expr_raw) <- expr_raw[,1]
  expr_mat <- expr_raw[,-1, drop=FALSE]
  expr_mat <- as.matrix(expr_mat)
  expr_mat <- apply(expr_mat, 2, as.numeric)
  rownames(expr_mat) <- rownames(expr_raw)
  combined_expr <- expr_mat
  cat("维度:", nrow(expr_mat), "基因 x", ncol(expr_mat), "样本\n")
  
  # --- 自动分组 ---
  sample_names <- colnames(expr_mat)
  group_info <- ifelse(grepl("_con$", sample_names, ignore.case=TRUE), "Control",
                       ifelse(grepl("_tre$", sample_names, ignore.case=TRUE), "Disease", "Unknown"))
  if(any(group_info=="Unknown")) stop("分组未知！")
  num_ctrl  <- sum(group_info=="Control")
  num_treat <- sum(group_info=="Disease")
  cat("Control:", num_ctrl, "  Disease:", num_treat, "\n")
  
  # --- limma DEG ---
  group_labels <- factor(group_info, levels=c("Control","Disease"))
  design_mat <- model.matrix(~0+group_labels)
  colnames(design_mat) <- c("Control","Disease")
  fit <- lmFit(expr_mat, design_mat)
  fit2 <- contrasts.fit(fit, makeContrasts(Disease-Control, levels=design_mat))
  fit2 <- eBayes(fit2)
  all_diff_results <- topTable(fit2, adjust.method="fdr", number=Inf)
  cat("adj.P<0.05 基因数:", sum(all_diff_results$adj.P.Val<0.05), "\n")
  
  # 保存全基因结果
  fname <- paste0("DE_results_", disease_label, ".csv")
  write.table(cbind(Gene=rownames(all_diff_results), all_diff_results),
              file=fname, sep=",", quote=FALSE, row.names=FALSE)
  cat("保存:", fname, "\n")
  
  # --- 筛选显著 DEG ---
  significant_DEGs <- all_diff_results[abs(all_diff_results$logFC) > threshold_logFC &
                                         all_diff_results$adj.P.Val < threshold_adjP, ]
  cat("显著 DEG (|logFC|>", threshold_logFC, " & adjP<", threshold_adjP, "): ",
      nrow(significant_DEGs),
      " (上调", sum(significant_DEGs$logFC>0),
      ", 下调", sum(significant_DEGs$logFC<0), ")\n")
  
  output_DEGs <- cbind(Gene=rownames(significant_DEGs), significant_DEGs)
  SE <- ifelse(as.numeric(output_DEGs[,"t"]) != 0,
               abs(as.numeric(output_DEGs[,"logFC"]) / as.numeric(output_DEGs[,"t"])), NA)
  output_DEGs <- cbind(output_DEGs, SE=SE)
  output_DEGs <- as.data.frame(output_DEGs, stringsAsFactors=FALSE)
  output_DEGs$Regulation <- ifelse(as.numeric(output_DEGs$logFC)>0, "Up", "Down")
  desired_order <- c("Gene","logFC","Regulation","SE","AveExpr","t","P.Value","adj.P.Val","B")
  existing_cols <- intersect(desired_order, colnames(output_DEGs))
  output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]
  
  fname2 <- paste0("DE_significant_", disease_label, ".csv")
  write.csv(output_DEGs, file=fname2, row.names=FALSE)
  cat("保存:", fname2, "\n")
  
  return(list(
    expr = combined_expr,
    deg = all_diff_results,
    sig = significant_DEGs,
    sig_output = output_DEGs,
    group_info = group_info,
    num_ctrl = num_ctrl,
    num_treat = num_treat,
    label = disease_label
  ))
}

# ==================== CRC 完整流程 ====================
CRC <- run_one_disease("merged_after_batch_removal.csv", "CRC")

# ---- CRC 渐变色火山图（Step 8 模板）----
cat("\n[Step 8] CRC 渐变色火山图...\n")
volcano_data <- CRC$deg
volcano_data$Gene <- rownames(volcano_data)
intersect_deg_data <- volcano_data[volcano_data$Gene %in% rownames(CRC$sig), ]
if (nrow(intersect_deg_data) <= 40) {
  label_genes_gradient <- intersect_deg_data$Gene
} else {
  up_int   <- intersect_deg_data[intersect_deg_data$logFC > 0, ]
  up_int   <- up_int[order(-(-log10(up_int$adj.P.Val))), ]
  down_int <- intersect_deg_data[intersect_deg_data$logFC < 0, ]
  down_int <- down_int[order(-(-log10(down_int$adj.P.Val))), ]
  up_label   <- head(up_int$Gene, min(20, nrow(up_int)))
  down_label <- head(down_int$Gene, min(20, nrow(down_int)))
  label_genes_gradient <- unique(c(up_label, down_label))
}
volcano_data$label_gene <- ifelse(volcano_data$Gene %in% label_genes_gradient, volcano_data$Gene, "")
x_limit <- max(abs(volcano_data$logFC), na.rm=TRUE) * 1.1
y_max   <- max(-log10(volcano_data$adj.P.Val), na.rm=TRUE) * 1.05
y_vals  <- -log10(volcano_data$adj.P.Val); y_vals <- y_vals[is.finite(y_vals)]
size_breaks <- pretty(range(y_vals, na.rm=TRUE), n=5)
size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max(y_vals)]
up_count   <- sum(volcano_data$logFC > threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)
down_count <- sum(volcano_data$logFC < -threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)

gradient_volcano <- ggplot(volcano_data, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=logFC, size=-log10(adj.P.Val)), alpha=0.85) +
  geom_point(data=subset(volcano_data, label_gene!=""),
             aes(size=-log10(adj.P.Val)), shape=21, fill=NA, color="black", stroke=1.2) +
  scale_color_gradientn(
    colors=c("#2B83BA","#5AAE61","#A6D96A","#FFFFBF","#FEE08B","#FDAE61","#D7191C"),
    values=scales::rescale(c(-6,-4,-2,0,2,4,6)), limits=c(-x_limit,x_limit), name="log2FC") +
  scale_size_continuous(range=c(0.8,7), name="-log10(p_val)",
                        breaks=size_breaks, labels=as.character(round(size_breaks,0))) +
  geom_vline(xintercept=c(-threshold_logFC, threshold_logFC), linetype="dashed", color="grey50", linewidth=0.6) +
  geom_hline(yintercept=-log10(threshold_adjP), linetype="dashed", color="grey50", linewidth=0.6) +
  geom_text_repel(data=subset(volcano_data, label_gene!=""), aes(label=label_gene),
                  size=3.2, max.overlaps=Inf, box.padding=0.5, point.padding=0.3,
                  segment.color="grey40", segment.size=0.3, fontface="italic", color="black", show.legend=FALSE) +
  annotate("segment", x=-x_limit*0.55, xend=-x_limit*0.9, y=y_max*0.97, yend=y_max*0.97,
           arrow=arrow(length=unit(0.25,"cm"), type="closed"), color="#2B83BA", linewidth=1.2) +
  annotate("text", x=-x_limit*0.72, y=y_max*0.97, label=paste0("Down (",down_count,")"),
           color="#2B83BA", fontface="bold", size=4.5, vjust=-0.8) +
  annotate("segment", x=x_limit*0.55, xend=x_limit*0.9, y=y_max*0.97, yend=y_max*0.97,
           arrow=arrow(length=unit(0.25,"cm"), type="closed"), color="#D7191C", linewidth=1.2) +
  annotate("text", x=x_limit*0.72, y=y_max*0.97, label=paste0("Up (",up_count,")"),
           color="#D7191C", fontface="bold", size=4.5, vjust=-0.8) +
  annotate("text", x=x_limit*0.98, y=-log10(threshold_adjP),
           label=paste0("p=",threshold_adjP), hjust=1, vjust=-0.5, size=3.5, color="grey30") +
  scale_x_continuous(limits=c(-x_limit,x_limit), expand=c(0.02,0)) +
  scale_y_continuous(limits=c(0,y_max), expand=c(0.02,0)) +
  labs(title=paste0("CRC Volcano Plot (", nrow(CRC$sig), " DEGs)"),
       x="log2FC", y="-log10(adj.P)") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold",hjust=0.5,size=16),
        axis.title=element_text(face="bold",size=12),
        axis.text=element_text(color="black",size=10),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",linewidth=1.2),
        legend.position="right",
        legend.background=element_rect(fill="white",color="grey80",linewidth=0.5),
        legend.key=element_rect(fill="white"),
        legend.title=element_text(face="bold",size=10),
        legend.text=element_text(size=9),
        legend.spacing.y=unit(0.3,"cm"),
        plot.margin=margin(15,15,10,10)) +
  guides(color=guide_colorbar(title="log2FC",title.position="top",title.hjust=0.5,
         barwidth=1.2,barheight=10,frame.colour="black",frame.linewidth=0.5,
         ticks.colour="black",ticks.linewidth=0.5,order=1),
         size=guide_legend(title="-log10(p_val)",title.position="top",title.hjust=0.5,
         override.aes=list(alpha=1,color="grey30"),order=2))
ggsave(paste0("CRC_volcano_gradient.pdf"), gradient_volcano, width=9, height=7.5, dpi=300)

# ---- CRC 基因排序点图 ----
gene_data_real <- CRC$deg[order(CRC$deg$logFC, decreasing=TRUE),]
plot_rank <- ggplot(gene_data_real, aes(x=seq_along(logFC), y=logFC)) +
  geom_point(aes(color=adj.P.Val, size=abs(logFC)), alpha=0.8) +
  scale_color_gradientn(colors=rev(brewer.pal(9,"RdYlBu")),
                        values=scales::rescale(c(0,0.25,0.5,0.75,1))) +
  scale_size_continuous(range=c(1,3)) +
  labs(title="CRC Gene Ranking Dotplot", x="Ranking of DEGs", y="Log2 Fold Change",
       color="Adjusted P-value", size="Log2 FC") +
  theme_minimal(base_size=16) +
  theme(plot.title=element_text(size=18,face="bold",color="#2F4F4F",hjust=0.5),
        axis.title=element_text(size=16,face="bold",color="#2F4F4F"),
        axis.text=element_text(size=14,color="#555555"),
        legend.position="right", legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        panel.grid=element_blank(), axis.line=element_line(color="black",linewidth=0.5))
ggsave("CRC_gene_ranking_dotplot.pdf", plot_rank, width=10, height=8, dpi=300)

# ---- CRC PCA ----
pca_res <- prcomp(t(CRC$expr), scale.=TRUE)
pca_var_pct <- round(100*pca_res$sdev^2/sum(pca_res$sdev^2), 2)
pca_df <- data.frame(Sample=colnames(CRC$expr),
                     PC1=pca_res$x[,1], PC2=pca_res$x[,2],
                     Group=factor(CRC$group_info, levels=c("Control","Disease")))
pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Group, fill=Group)) +
  stat_ellipse(geom="polygon", level=0.95, alpha=0.2, linewidth=1) +
  geom_point(size=3, shape=15) +
  scale_color_manual(values=c(Control="#E64B35", Disease="#4DBBD5")) +
  scale_fill_manual(values=c(Control="#E64B35", Disease="#4DBBD5")) +
  labs(x=paste0("PC1(", pca_var_pct[1], "%)"), y=paste0("PC2(", pca_var_pct[2], "%)")) +
  theme_bw(base_size=14) +
  theme(axis.title=element_text(face="bold",size=14,color="black"),
        axis.text=element_text(size=12,color="black"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=11),
        legend.position="right",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",linewidth=1))
ggsave("CRC_PCA.pdf", pca_plot, width=7, height=6, dpi=300)

# ---- CRC ComplexHeatmap（top 50上调+50下调显著DEG）----
if (nrow(CRC$sig) > 0) {
  cat("\n[Step 11] CRC ComplexHeatmap...\n")
  ordered_DEGs <- CRC$sig[order(CRC$sig$logFC), ]
  ordered_names <- rownames(ordered_DEGs)
  n_total <- length(ordered_names)
  if (n_total > (max_display_genes*2)) {
    sel_genes <- ordered_names[c(1:max_display_genes, (n_total-max_display_genes+1):n_total)]
  } else {
    sel_genes <- ordered_names
  }
  heatmap_expr <- CRC$expr[sel_genes, ]
  heatmap_expr_scaled <- t(scale(t(heatmap_expr)))
  
  sample_annotation <- data.frame(Group=factor(c(rep("Control",CRC$num_ctrl), rep("Disease",CRC$num_treat))))
  rownames(sample_annotation) <- colnames(CRC$expr)
  annotation_colors <- list(Group=c(Control="#66C2A5", Disease="#FC8D62"))
  
  ha_top <- HeatmapAnnotation(
    Distribution=anno_density(heatmap_expr_scaled, type="heatmap", which="column", height=unit(2,"cm")),
    Group=sample_annotation$Group, col=annotation_colors,
    annotation_name_side="left", annotation_name_gp=gpar(fontsize=12))
  
  color_pal <- colorRamp2(c(-2,0,2), c("#313695","white","#A50026"))
  gene_logFC <- CRC$sig[sel_genes, "logFC"]
  gene_labels <- paste0(sel_genes, ifelse(gene_logFC>0," (Up)"," (Down)"))
  names(gene_labels) <- sel_genes
  
  ht <- Heatmap(heatmap_expr_scaled, name="Expression\n(Z-score)", col=color_pal,
                top_annotation=ha_top, cluster_columns=FALSE, show_column_names=FALSE,
                row_names_gp=gpar(fontsize=8), row_labels=gene_labels[rownames(heatmap_expr_scaled)],
                column_title=paste("CRC DEG Heatmap\nControl:", CRC$num_ctrl, " | Disease:", CRC$num_treat),
                column_title_gp=gpar(fontsize=14, fontface="bold"),
                heatmap_legend_param=list(title_gp=gpar(fontsize=12,fontface="bold"), labels_gp=gpar(fontsize=10)))
  
  pdf("CRC_heatmap_with_distribution.pdf", width=12, height=10)
  tryCatch({draw(ht)}, error=function(e){cat("ComplexHeatmap error:", e$message, "\n"); print(ht)})
  dev.off()
  cat("保存: CRC_heatmap_with_distribution.pdf\n")
}

# ==================== MetS DEG ====================
MetS <- run_one_disease("GSE98895.csv", "MetS")
# MetS 出DEG表格 + 火山图（风格与CRC一致，供补充材料使用）

# ---- MetS 渐变色火山图 ----
cat("\n[Step 8b] MetS 渐变色火山图...\n")
volcano_data_mets <- MetS$deg
volcano_data_mets$Gene <- rownames(volcano_data_mets)
intersect_deg_data_mets <- volcano_data_mets[volcano_data_mets$Gene %in% rownames(MetS$sig), ]
if (nrow(intersect_deg_data_mets) <= 40) {
  label_genes_mets <- intersect_deg_data_mets$Gene
} else {
  up_int_m   <- intersect_deg_data_mets[intersect_deg_data_mets$logFC > 0, ]
  up_int_m   <- up_int_m[order(-(-log10(up_int_m$adj.P.Val))), ]
  down_int_m <- intersect_deg_data_mets[intersect_deg_data_mets$logFC < 0, ]
  down_int_m <- down_int_m[order(-(-log10(down_int_m$adj.P.Val))), ]
  up_label_m   <- head(up_int_m$Gene, min(20, nrow(up_int_m)))
  down_label_m <- head(down_int_m$Gene, min(20, nrow(down_int_m)))
  label_genes_mets <- unique(c(up_label_m, down_label_m))
}
volcano_data_mets$label_gene <- ifelse(volcano_data_mets$Gene %in% label_genes_mets,
                                       volcano_data_mets$Gene, "")
x_limit_m <- max(abs(volcano_data_mets$logFC), na.rm=TRUE) * 1.1
y_max_m   <- max(-log10(volcano_data_mets$adj.P.Val), na.rm=TRUE) * 1.05
y_vals_m  <- -log10(volcano_data_mets$adj.P.Val); y_vals_m <- y_vals_m[is.finite(y_vals_m)]
size_breaks_m <- pretty(range(y_vals_m, na.rm=TRUE), n=5)
size_breaks_m <- size_breaks_m[size_breaks_m > 0 & size_breaks_m <= max(y_vals_m)]
up_count_m   <- sum(volcano_data_mets$logFC > threshold_logFC & volcano_data_mets$adj.P.Val < threshold_adjP)
down_count_m <- sum(volcano_data_mets$logFC < -threshold_logFC & volcano_data_mets$adj.P.Val < threshold_adjP)

gradient_volcano_mets <- ggplot(volcano_data_mets, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=logFC, size=-log10(adj.P.Val)), alpha=0.85) +
  geom_point(data=subset(volcano_data_mets, label_gene!=""),
             aes(size=-log10(adj.P.Val)), shape=21, fill=NA, color="black", stroke=1.2) +
  scale_color_gradientn(
    colors=c("#2B83BA","#5AAE61","#A6D96A","#FFFFBF","#FEE08B","#FDAE61","#D7191C"),
    values=scales::rescale(c(-6,-4,-2,0,2,4,6)), limits=c(-x_limit_m,x_limit_m), name="log2FC") +
  scale_size_continuous(range=c(0.8,7), name="-log10(p_val)",
                        breaks=size_breaks_m, labels=as.character(round(size_breaks_m,0))) +
  geom_vline(xintercept=c(-threshold_logFC, threshold_logFC), linetype="dashed", color="grey50", linewidth=0.6) +
  geom_hline(yintercept=-log10(threshold_adjP), linetype="dashed", color="grey50", linewidth=0.6) +
  geom_text_repel(data=subset(volcano_data_mets, label_gene!=""), aes(label=label_gene),
                  size=3.2, max.overlaps=Inf, box.padding=0.5, point.padding=0.3,
                  segment.color="grey40", segment.size=0.3, fontface="italic", color="black", show.legend=FALSE) +
  annotate("segment", x=-x_limit_m*0.55, xend=-x_limit_m*0.9, y=y_max_m*0.97, yend=y_max_m*0.97,
           arrow=arrow(length=unit(0.25,"cm"), type="closed"), color="#2B83BA", linewidth=1.2) +
  annotate("text", x=-x_limit_m*0.72, y=y_max_m*0.97, label=paste0("Down (",down_count_m,")"),
           color="#2B83BA", fontface="bold", size=4.5, vjust=-0.8) +
  annotate("segment", x=x_limit_m*0.55, xend=x_limit_m*0.9, y=y_max_m*0.97, yend=y_max_m*0.97,
           arrow=arrow(length=unit(0.25,"cm"), type="closed"), color="#D7191C", linewidth=1.2) +
  annotate("text", x=x_limit_m*0.72, y=y_max_m*0.97, label=paste0("Up (",up_count_m,")"),
           color="#D7191C", fontface="bold", size=4.5, vjust=-0.8) +
  annotate("text", x=x_limit_m*0.98, y=-log10(threshold_adjP),
           label=paste0("p=",threshold_adjP), hjust=1, vjust=-0.5, size=3.5, color="grey30") +
  scale_x_continuous(limits=c(-x_limit_m,x_limit_m), expand=c(0.02,0)) +
  scale_y_continuous(limits=c(0,y_max_m), expand=c(0.02,0)) +
  labs(title=paste0("MetS Volcano Plot (", nrow(MetS$sig), " DEGs)"),
       x="log2FC", y="-log10(adj.P)") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold",hjust=0.5,size=16),
        axis.title=element_text(face="bold",size=12),
        axis.text=element_text(color="black",size=10),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",linewidth=1.2),
        legend.position="right",
        legend.background=element_rect(fill="white",color="grey80",linewidth=0.5),
        legend.key=element_rect(fill="white"),
        legend.title=element_text(face="bold",size=10),
        legend.text=element_text(size=9),
        legend.spacing.y=unit(0.3,"cm"),
        plot.margin=margin(15,15,10,10)) +
  guides(color=guide_colorbar(title="log2FC",title.position="top",title.hjust=0.5,
         barwidth=1.2,barheight=10,frame.colour="black",frame.linewidth=0.5,
         ticks.colour="black",ticks.linewidth=0.5,order=1),
         size=guide_legend(title="-log10(p_val)",title.position="top",title.hjust=0.5,
         override.aes=list(alpha=1,color="grey30"),order=2))
ggsave("MetS_volcano_gradient.pdf", gradient_volcano_mets, width=9, height=7.5, dpi=300)
cat("保存: MetS_volcano_gradient.pdf\n")

# ==================== co-DEGs = CRC ∩ MetS ====================
cat("\n", paste(rep("=",60),collapse=""), "\n")
cat("  co-DEGs: CRC DEG ∩ MetS DEG\n")
cat(paste(rep("=",60),collapse=""), "\n")

crc_deg_genes  <- rownames(CRC$sig)
mets_deg_genes <- rownames(MetS$sig)
co_genes <- intersect(crc_deg_genes, mets_deg_genes)
cat("CRC 显著:", length(crc_deg_genes), "  MetS 显著:", length(mets_deg_genes), "\n")
cat("co-DEGs (交集):", length(co_genes), "\n")

if (length(co_genes) > 0) {
  # 保存 co-DEGs 表
  co_df <- data.frame(
    Gene = co_genes,
    CRC_log2FC  = CRC$deg[co_genes, "logFC"],
    CRC_adjP    = CRC$deg[co_genes, "adj.P.Val"],
    MetS_log2FC = MetS$deg[co_genes, "logFC"],
    MetS_adjP   = MetS$deg[co_genes, "adj.P.Val"],
    stringsAsFactors = FALSE
  )
  write.csv(co_df, "co_DEGs_CRC_MetS.csv", row.names=FALSE)
  cat("保存: co_DEGs_CRC_MetS.csv\n")
  
  # TXT 基因清单(后续GO/KEGG富集用)
  write.table(data.frame(gene=co_genes), "co_DEGs_geneList.txt",
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  cat("保存: co_DEGs_geneList.txt\n")
  
  # --- Venn图 CRC DEG ∩ MetS DEG ---
  if (require("ggvenn", quietly=TRUE)) {
    venn_list <- list(CRC=length(crc_deg_genes), MetS=length(mets_deg_genes),
                      intersection=length(co_genes))
    # 简单文字版（ggvenn需要完整基因列表）
    gene_list_venn <- list(CRC=crc_deg_genes, MetS=mets_deg_genes)
    pdf("co_DEGs_CRC_MetS_Venn.pdf", width=8, height=8)
    print(ggvenn(gene_list_venn, show_percentage=TRUE, stroke_color="white", stroke_size=0.5,
                 fill_color=c("#E64B35","#4DBBD5"), set_name_color=c("#E64B35","#4DBBD5"),
                 set_name_size=6, text_size=4.5))
    dev.off()
    cat("保存: co_DEGs_CRC_MetS_Venn.pdf\n")
  }
  
  # --- 用co-DEGs在CRC火山图上重新标注 ---
  cat("\nCRC火山图(co-DEGs标注版)...\n")
  v2 <- CRC$deg; v2$Gene <- rownames(v2)
  v2$Group <- "Not Significant"
  v2$Group[v2$logFC > threshold_logFC & v2$adj.P.Val < threshold_adjP] <- "Up regulated"
  v2$Group[v2$logFC < -threshold_logFC & v2$adj.P.Val < threshold_adjP] <- "Down regulated"
  v2$label <- ifelse(v2$Gene %in% co_genes, v2$Gene, "")
  
  gc <- table(v2$Group)
  up_c   <- ifelse("Up regulated"   %in% names(gc), gc["Up regulated"], 0)
  down_c <- ifelse("Down regulated" %in% names(gc), gc["Down regulated"], 0)
  not_c  <- ifelse("Not Significant"%in% names(gc), gc["Not Significant"], 0)
  leg_labels <- c(paste0("Up regulated (",up_c,")"),
                  paste0("Down regulated (",down_c,")"),
                  paste0("Not Significant (",not_c,")"))
  
  vp <- ggplot(v2, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
    geom_point(alpha=0.7, size=2) +
    geom_text_repel(data=subset(v2, label!=""), aes(label=label),
                    size=3, max.overlaps=Inf, box.padding=0.3, point.padding=0.2,
                    segment.color="grey50", color="black", show.legend=FALSE) +
    scale_color_manual(values=c("Up regulated"="#FF4500","Down regulated"="#1E90FF","Not Significant"="#808080"),
                       breaks=c("Up regulated","Down regulated","Not Significant"), labels=leg_labels) +
    geom_vline(xintercept=c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black") +
    geom_hline(yintercept=-log10(threshold_adjP), linetype="dashed", color="black") +
    labs(title=paste0("CRC Volcano (co-DEGs with MetS: ", length(co_genes), " genes)"),
         x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
    theme_minimal(base_size=14) +
    theme(plot.title=element_text(face="bold",hjust=0.5,color="#2F4F4F"),
          panel.grid=element_blank(), axis.line=element_line(color="black",linewidth=0.5))
  ggsave("CRC_volcano_coDEGs_labeled.pdf", vp, width=8, height=8, dpi=300)
  cat("保存: CRC_volcano_coDEGs_labeled.pdf\n")
  
  # --- co-DEGs 热图（CRC侧） ---
  cat("co-DEGs 热图(CRC侧)...\n")
  plot_genes <- co_genes[co_genes %in% rownames(CRC$expr)]
  if (length(plot_genes) > 0) {
    heatmap_mat <- CRC$expr[plot_genes, , drop=FALSE]
    heatmap_scaled <- t(scale(t(heatmap_mat)))
    
    ha_col <- HeatmapAnnotation(
      Group=factor(c(rep("Control",CRC$num_ctrl), rep("Disease",CRC$num_treat))),
      col=list(Group=c(Control="#2166AC", Disease="#B2182B")),
      annotation_name_side="left", annotation_name_gp=gpar(fontsize=11))
    col_fun <- colorRamp2(c(-2,0,2), c("#2166AC","white","#B2182B"))
    
    ht_int <- Heatmap(heatmap_scaled, name="Z-score", col=col_fun,
                      top_annotation=ha_col, cluster_columns=FALSE, cluster_rows=TRUE,
                      show_column_names=FALSE, row_names_gp=gpar(fontsize=10,fontface="italic"),
                      column_title=paste0("co-DEGs Heatmap (n=", length(plot_genes), ")"),
                      column_title_gp=gpar(fontsize=14,fontface="bold"),
                      heatmap_legend_param=list(title_gp=gpar(fontsize=11,fontface="bold"),labels_gp=gpar(fontsize=9)))
    
    pdf("co_DEGs_heatmap.pdf", width=10, height=max(4, length(plot_genes)*0.45))
    draw(ht_int)
    dev.off()
    cat("保存: co_DEGs_heatmap.pdf\n")
    
    # --- co-DEGs 箱线图(CRC侧, top30) ---
    if (require("ggpubr", quietly=TRUE)) {
      top_plot <- head(plot_genes, 30)
      group_vec <- factor(c(rep("Control",CRC$num_ctrl), rep("Disease",CRC$num_treat)), levels=c("Control","Disease"))
      box_list <- list()
      for (g in top_plot) {
        box_list[[g]] <- data.frame(Gene=g, Expression=as.numeric(CRC$expr[g,]), Group=group_vec)
      }
      box_df <- do.call(rbind, box_list)
      box_df$Gene <- factor(box_df$Gene, levels=top_plot)
      
      pvals <- CRC$deg[top_plot, "adj.P.Val"]
      fmt_p <- function(p){if(p<0.001)"***" else if(p<0.01)"**" else if(p<0.05)"*" else "ns"}
      p_labs <- sapply(pvals, fmt_p)
      ypos <- sapply(top_plot, function(g) max(as.numeric(CRC$expr[g,]))*1.08)
      anno_df <- data.frame(Gene=factor(top_plot,levels=top_plot), label=p_labs, y=ypos)
      
      bp <- ggplot(box_df, aes(x=Gene, y=Expression, fill=Group)) +
        geom_boxplot(outlier.size=0.8, width=0.65, alpha=0.85, position=position_dodge(0.75)) +
        scale_fill_manual(values=c(Control="#2166AC", Disease="#B2182B")) +
        geom_text(data=anno_df, aes(x=Gene, y=y, label=label), inherit.aes=FALSE, size=4, fontface="bold") +
        labs(x="", y="Expression", fill="Group") +
        theme_bw(base_size=13) +
        theme(axis.text.x=element_text(size=11,color="black",face="bold.italic",angle=45,hjust=1),
              axis.text.y=element_text(size=10,color="black"),
              axis.title.y=element_text(face="bold",size=13),
              legend.position="top", legend.title=element_text(face="bold",size=11),
              legend.text=element_text(size=10),
              panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.border=element_rect(color="black",linewidth=1))
      
      pw <- max(6, length(top_plot)*1.2+2)
      ggsave("co_DEGs_boxplot.pdf", bp, width=pw-5, height=5.5, dpi=300)
      cat("保存: co_DEGs_boxplot.pdf\n")
    }
  }
} else {
  cat("\n⚠️ 警告: CRC 和 MetS 没有共同显著差异基因！\n")
}

cat("\n", paste(rep("=",60),collapse=""), "\n")
cat("全部完成!\n")
cat("CRC 显著 DEG:", nrow(CRC$sig), "\n")
cat("MetS 显著 DEG:", nrow(MetS$sig), "\n")
cat("co-DEGs:", length(co_genes), "\n")
cat(paste(rep("=",60),collapse=""), "\n")
