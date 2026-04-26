# ============================================================
# Script 10: 差异基因分析（多层次）
# ============================================================
source("Script00_Setup.R")

library(ggvenn)
library(fgsea)

setwd("/data/home/caoy/Xinjiang_AAD/gene/")
load("../ad_Anchors_sctgroup2.RData")

Idents(ad) <- "new_seurat_clusters"
DefaultAssay(ad) <- "SCT"

# ==================================================
# 第一部分：全样本各cluster的marker基因
# ==================================================
ad_markers <- FindAllMarkers(ad,
                             test.use       = "wilcox",
                             only.pos       = FALSE,
                             min.pct        = 0.25,
                             logfc.threshold = 0.25)
write.csv(ad_markers, file="./ad_markers.csv", row.names=FALSE)
message("全样本markers：", nrow(ad_markers), "条")

# ==================================================
# 第二部分：按疾病分组的DEG
# ==================================================
ad@meta.data$group3 <- NA
group3_map <- c(
  s1p2="wd", s2p2="wd", s3p1="wd", s4p3="wd",
  s5p4="wd", s6p4="wd", s9p5="wd", s10p6="wd", s11p6="wd",
  s12p7="tbg", s13p7="tbg", s14p8="tbg", s15p8="tbg",
  s16p4="zsgx", s17p4="zsgx",
  s18p4="zjz", s19p4="zjz"
)
for(s in names(group3_map)){
  ad@meta.data$group3[ad@meta.data$sample==s] <- group3_map[s]
}

wd   <- subset(x=ad, subset=(group3=="wd"))
tbg  <- subset(x=ad, subset=(group3=="tbg"))
zjz  <- subset(x=ad, subset=(group3=="zjz"))
zsgx <- subset(x=ad, subset=(group3=="zsgx"))

# 批量FindAllMarkers
group_list <- list(wd=wd, tbg=tbg, zjz=zjz, zsgx=zsgx)
deg_results <- lapply(group_list, function(obj){
  FindAllMarkers(obj, test.use="wilcox", only.pos=FALSE,
                 min.pct=0.25, logfc.threshold=0.25)
})

write.csv(deg_results[["wd"]],   file="./aa_markers.csv",   row.names=FALSE)
write.csv(deg_results[["tbg"]],  file="./ba_markers.csv",   row.names=FALSE)
write.csv(deg_results[["zjz"]],  file="./lcca_markers.csv", row.names=FALSE)
write.csv(deg_results[["zsgx"]], file="./lsa_markers.csv",  row.names=FALSE)

# ==================================================
# 第三部分：韦恩图展示各组DEG交集
# ==================================================
gene_list_venn <- list(
  "AA"   = deg_results[["wd"]]$gene,
  "BA"   = deg_results[["tbg"]]$gene,
  "LCCA" = deg_results[["zjz"]]$gene,
  "LSA"  = deg_results[["zsgx"]]$gene
)

pdf("./cluster_total_venn.pdf", width=2.6, height=3.9)
ggvenn(gene_list_venn, stroke_color="white",
       set_name_size=3, stroke_size=0.25, text_size=3)
dev.off()

# ==================================================
# 第四部分：中膜各细胞类型在不同分组中的DEG
# ==================================================
load("../fig2/20231204/wd.RData")

# 提取中膜各细胞类型的子集
Intima    <- subset(x=wd, loc=="Intima")
Media     <- subset(x=wd, loc=="Media")
Adventitia <- subset(x=wd, loc=="Adventitia")

M_Fb  <- subset(Media, max_cell=="Fb")
M_MSC <- subset(Media, max_cell=="MSC")
M_SMC <- subset(Media, max_cell=="SMC")
M_Mp  <- subset(Media, max_cell=="Mp")

# 只保留比例>1的spot
M_SMC1 <- subset(M_SMC, SMC>=1)
M_Fb1  <- subset(M_Fb,  Fb>=1)
M_Mp1  <- subset(M_Mp,  Mp>=1)
M_MSC1 <- subset(M_MSC, MSC>=1)

cell_subset_list <- list(
  M_Fb=M_Fb1, M_MSC=M_MSC1, M_SMC=M_SMC1, M_Mp=M_Mp1
)

M_deg_list <- lapply(cell_subset_list, function(obj){
  FindAllMarkers(obj, only.pos=FALSE,
                 min.pct=0.25, logfc.threshold=0.25)
})

write.csv(M_deg_list[["M_Fb"]],  file="./M.Fb.csv",  row.names=FALSE)
write.csv(M_deg_list[["M_MSC"]], file="./M.MSC.csv", row.names=FALSE)
write.csv(M_deg_list[["M_SMC"]], file="./M.SMC.csv", row.names=FALSE)
write.csv(M_deg_list[["M_Mp"]],  file="./M.Mp.csv",  row.names=FALSE)

# ==================================================
# 第五部分：top10基因热图
# ==================================================
library(ComplexHeatmap)

top10_per_celltype <- lapply(M_deg_list, function(deg){
  sig_deg <- deg[deg$p_val < 0.05, ]
  as.data.frame(sig_deg %>%
    group_by(cluster) %>%
    slice_max(n=5, order_by=avg_log2FC))
})

# 以M_SMC为例绘制热图
gene_cell_exp <- AverageExpression(
  M_SMC, features=top10_per_celltype[["M_SMC"]]$gene,
  group.by='group'
)$SCT

library(vegan)
df_scaled <- decostand(gene_cell_exp, "standardize", MARGIN=1)

pdf("./Media_SMC_deg_Heatmap.pdf", width=4, height=5)
Heatmap(df_scaled,
        cluster_rows=TRUE, cluster_columns=FALSE,
        show_column_names=TRUE, show_row_names=TRUE,
        col=colorRampPalette(c("#0033FF","#D6E7A3","#FF0000"))(80),
        row_names_gp=gpar(fontsize=5),
        column_names_gp=gpar(fontsize=5))
dev.off()

message("Script 10 完成：差异基因分析完毕")