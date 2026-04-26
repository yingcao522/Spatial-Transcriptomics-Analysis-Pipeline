# ============================================================
# Script 04: Seurat Anchor Label Transfer（空间 ← scRNA标签）
# ============================================================
source("Script00_Setup.R")
setwd("/data/home/caoy/Xinjiang_AAD/fig3/")

load("./ad_clustered.RData")
load("../taa_combined.RData")

# ---- FindTransferAnchors ----
message("计算Transfer Anchors...")
anchors <- FindTransferAnchors(
  reference = combined,
  query     = ad,
  normalization.method = "SCT"
)

# ---- TransferData ----
# 注意：refdata使用粗分类celltype
predictions <- TransferData(
  anchorset        = anchors,
  refdata          = combined$celltype,
  prediction.assay = TRUE,
  weight.reduction = ad[["harmony"]],
  dims             = 1:30
)

# 挂载预测结果
ad[["predictions"]] <- predictions

# ---- 保存预测结果 ----
write.csv(ad@assays$predictions@data,
          file="./tsne_umap/ad_Anchors_predictions_sctgroup2.csv")
save(ad, file="./ad_Anchors_sctgroup2.RData")

# ---- 可视化预测分数 ----
ctype <- c("Tcell","MonoMaphDC","NK","Plasma","Mastcell","Bcell",
           "MSC","SMC2","Fibroblast","SMC1","EC")

DefaultAssay(ad) <- "predictions"
pdf("./tsne_umap/ad_Anchors_SpatialFeaturePlot.pdf",width=60,height=30)
SpatialFeaturePlot(object=ad, features=ctype, pt.size.factor=1.6, crop=TRUE)
dev.off()

# ---- 整合预测分数到metadata ----
cell_an <- t(as.data.frame(ad@assays$predictions@data))
cell_an <- as.data.frame(cell_an)

# 合并SMC1和SMC2
cell_an$SMC <- cell_an$SMC1 + cell_an$SMC2
# 去除不需要的列（根据实际列名调整）
cell_an <- cell_an[, !colnames(cell_an) %in% c("SMC1","SMC2","Mastcell","max")]

cell_an$cells <- rownames(cell_an)

meta_data <- ad@meta.data
metadata  <- merge(meta_data, cell_an, by="cells", all=TRUE)
rownames(metadata) <- metadata[,1]
ad@meta.data <- metadata

DefaultAssay(ad) <- "SCT"

# ---- 定义 max_cell（占比最高的细胞类型）----
celltype_cols <- c("EC","Fb","MSC","SMC","MonoMaphDC","Plasma","Tcell","NK","Bcell")
avail_cols    <- intersect(celltype_cols, colnames(ad@meta.data))

ad@meta.data$max_cell <- apply(
  ad@meta.data[, avail_cols], 1,
  function(x) avail_cols[which.max(x)]
)
ad@meta.data$sec_cell <- apply(
  ad@meta.data[, avail_cols], 1,
  function(x) avail_cols[order(x, decreasing=TRUE)[2]]
)

# ---- 保存 ----
save(ad, file="./ad_Anchors_sctgroup2.RData")
message("Script 04 完成：Label Transfer 完毕")