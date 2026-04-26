# ============================================================
# Script 15: 同一患者(p4)不同血管节段比较
#   AA（升主动脉）vs LSA（左锁骨下）vs LCA（左颈总）
# ============================================================
source("Script00_Setup.R")

library(UpSetR)
library(AUCell)
library(GSEABase)
library(ggrepel)

setwd("/data/home/caoy/Xinjiang_AAD/fig5/")
load("../ad_Anchors_sctgroup2.RData")

# ==================================================
# 第一步：提取p4样本，添加分组标签
# ==================================================
p4 <- subset(ad, sample %in% c("s5p4","s6p4","s16p4","s17p4","s18p4","s19p4"))

p4$sample2 <- NA
p4$sample2[p4$sample %in% c("s5p4","s6p4")]   <- "AA"   # 升主动脉
p4$sample2[p4$sample %in% c("s16p4","s17p4")] <- "LSA"  # 左锁骨下动脉
p4$sample2[p4$sample %in% c("s18p4","s19p4")] <- "LCA"  # 左颈总动脉

# ==================================================
# 第二步：SCTransform + Harmony 批次整合
# ==================================================
p4 <- SCTransform(p4, assay="Spatial",
                  vars.to.regress="sample", verbose=FALSE)
p4 <- RunPCA(p4, assay="SCT", verbose=FALSE)
p4 <- RunHarmony(p4, group.by.vars="sample")
p4 <- FindNeighbors(p4, reduction="harmony", dims=1:30)
p4 <- FindClusters(p4, resolution=c(0.1,0.2,0.3))
p4 <- FindClusters(p4, resolution=0.1)
p4 <- RunTSNE(p4, reduction="harmony", dims=1:30)
p4 <- RunUMAP(p4, reduction="harmony", dims=1:30)

save(p4, file="./p4.RData")

# ==================================================
# 第三步：三组两两DEG（6次比较）
# ==================================================
Idents(p4) <- p4$sample2

comparisons_p4 <- list(
  c("LCA","AA"), c("LSA","AA"), c("LCA","LSA"),
  c("AA","LCA"), c("AA","LSA"), c("LSA","LCA")
)

deg_results_p4 <- list()
for(comp in comparisons_p4){
  comp_name <- paste0(comp[1], "_vs_", comp[2])
  message("比较:", comp_name)
  markers <- FindMarkers(
    object          = p4,
    ident.1         = comp[1],
    ident.2         = comp[2],
    assay           = "SCT",
    logfc.threshold = 0.25,
    min.pct         = 0.1,
    test.use        = "wilcox",
    only.pos        = FALSE,
    p.adjust.method = "BH"
  )
  deg_results_p4[[comp_name]] <- markers
  write.csv(markers, file=paste0("./markers_", comp_name, ".csv"))
}

# ==================================================
# 第四步：火山图（6张合并）
# ==================================================
plot_volcano_p4 <- function(markers, title){
  markers$gene <- rownames(markers)
  markers <- markers %>%
    filter(!is.na(avg_log2FC) & !is.na(p_val_adj)) %>%
    mutate(
      significant = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.25  ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      gene_label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 5, gene, "")
    )

  ggplot(markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=significant)) +
    geom_point(alpha=0.6, size=1.5) +
    scale_color_manual(values=c("Upregulated"="#FF9999",
                                "Downregulated"="#99CC00",
                                "Not Significant"="gray80")) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
    geom_vline(xintercept=c(-0.25,0.25), linetype="dashed", color="blue") +
    geom_text_repel(aes(label=gene_label), size=3, max.overlaps=20) +
    labs(title=title, x="Log2 Fold Change", y="-Log10 Adj P-value") +
    theme_bw() +
    theme(plot.title=element_text(size=8, hjust=0.5),
          legend.position="top")
}

plot_list <- mapply(
  function(res, nm) plot_volcano_p4(res, nm),
  deg_results_p4, names(deg_results_p4),
  SIMPLIFY=FALSE
)

pdf("./volcano_6comparisons.pdf", width=18, height=5)
plot_grid(plotlist=plot_list, nrow=2)
dev.off()

# ==================================================
# 第五步：UpSet图展示基因交集
# ==================================================
sig_gene_sets <- lapply(deg_results_p4, function(df){
  rownames(df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.25, ])
})

pdf("./upset_6comparisons.pdf", width=8, height=5)
upset(fromList(sig_gene_sets),
      main.bar.color="steelblue",
      sets.bar.color="darkred",
      order.by="freq",
      set_size.show=TRUE)
dev.off()

# ==================================================
# 第六步：AUCell 评分（用wd特异性基因集）
# ==================================================
# 读取wd1/wd2/wd3特异性基因
wd1_spec <- read.csv("../fig3/gene_avg_expression_wd1.csv", header=TRUE)
wd2_spec <- read.csv("../fig3/gene_avg_expression_wd2.csv", header=TRUE)
wd3_spec <- read.csv("../fig3/gene_avg_expression_wd3.csv", header=TRUE)

gene_sets_auc <- list(
  wd1_specific = GeneSet(wd1_spec$gene, setName="wd1_specific"),
  wd2_specific = GeneSet(wd2_spec$gene, setName="wd2_specific"),
  wd3_specific = GeneSet(wd3_spec$gene, setName="wd3_specific")
)

expr_matrix <- GetAssayData(p4, assay="SCT", layer="data")
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores=1, plotStats=FALSE)

auc_results <- list()
for(gs_name in names(gene_sets_auc)){
  cells_AUC <- AUCell_calcAUC(gene_sets_auc[[gs_name]], cells_rankings,
                               aucMaxRank=nrow(cells_rankings)*1)
  auc_results[[gs_name]] <- as.numeric(getAUC(cells_AUC))
}

# 添加到metadata
for(gs_name in names(auc_results)){
  p4[[paste0("AUCell_", gs_name)]] <- auc_results[[gs_name]]
}

# 可视化AUCell评分
score_features <- paste0("AUCell_", names(gene_sets_auc))

pdf("./AUCell_VlnPlot_p4.pdf", width=8, height=4)
VlnPlot(p4, features=score_features, group.by="sample2",
        pt.size=0.1, ncol=3)
dev.off()

pdf("./AUCell_SpatialFeaturePlot_p4.pdf", width=18, height=15)
SpatialFeaturePlot(p4, features=score_features, ncol=3)
dev.off()

# DotPlot: 三个节段 × 三个wd期特征评分
pdf("./AUCell_DotPlot_p4.pdf", width=5, height=3)
DotPlot(p4, features=score_features, group.by="sample2", assay="SCT") +
  scale_color_gradientn(colors=brewer.pal(9,"Blues")) +
  theme_bw() + coord_flip()
dev.off()

# 保存
save(p4, file="./p4_final.RData")
message("Script 15 完成：p4样本分析完毕")