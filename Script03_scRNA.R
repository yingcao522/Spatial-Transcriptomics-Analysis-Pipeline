# ============================================================
# Script 03: 单细胞RNA数据处理、整合、细胞类型注释
# ============================================================
source("Script00_Setup.R")
setwd("/data/home/caoy/Xinjiang_AAD/")

# ---- 加载单细胞数据 ----
load("./taa.Rdata")  # 包含TAA1~TAA8

# ---- 分样本标准化 ----
taa_list <- SplitObject(taa, split.by="sample")

taa_list <- lapply(taa_list, function(x){
  NormalizeData(x, verbose=FALSE)
})

# ---- Seurat整合 ----
message("计算整合锚点（耗时较长）...")
anchors <- FindIntegrationAnchors(
  object.list = taa_list,
  dims = 1:20
)

combined <- IntegrateData(anchorset=anchors, dims=1:20)

# ---- 修复metadata列名 ----
metadata <- combined@meta.data
names(metadata)[names(metadata)=="seurat_clusters"] <- "seurat_clusters_ori"
names(metadata)[names(metadata)=="integrated_snn_res.0.6"] <- "integrated_snn_res.0.6_ori"
combined@meta.data <- metadata

# ---- 质量过滤 ----
combined <- subset(x=combined,
                   subset = (nUMI  >= 200) &
                            (nGene >= 200) &
                            (nGene <= 5000) &
                            (mitoRatio < 0.10))

# ---- 标准化 ----
DefaultAssay(combined) <- "integrated"
combined <- SCTransform(combined, assay="RNA", verbose=FALSE)
combined <- ScaleData(combined, verbose=FALSE)

# ---- 降维与聚类 ----
combined <- RunPCA(combined, npcs=30, verbose=FALSE)
combined <- RunTSNE(combined, reduction="pca", dims=1:20)
combined <- FindNeighbors(combined, reduction="pca", dims=1:20)
combined <- FindClusters(combined, resolution=0.6)

# ---- 可视化 ----
pdf("./scRNA/combined_DimPlot_sample.pdf", width=5, height=4)
DimPlot(combined, reduction="tsne", group.by="sample")
dev.off()

pdf("./scRNA/combined_DimPlot_celltype.pdf", width=6, height=4)
DimPlot(combined, group.by="celltype2", label=TRUE)
dev.off()

# ---- 细胞类型统计（注释完成后）----
meta <- combined@meta.data
message("细胞类型分布：")
print(table(meta$celltype))

# 添加粗分类 (celltype3)
meta$celltype3 <- NA
meta$celltype3[str_detect(meta$celltype,"^Tcell")]       <- "Tcell"
meta$celltype3[str_detect(meta$celltype,"^MonoMaphDC")] <- "Mp"
meta$celltype3[str_detect(meta$celltype,"^Plasma")]      <- "Plasma"
meta$celltype3[str_detect(meta$celltype,"^MSC")]         <- "MSC"
meta$celltype3[str_detect(meta$celltype,"^SMC1")]        <- "SMC"
meta$celltype3[str_detect(meta$celltype,"^SMC2")]        <- "SMC"
meta$celltype3[str_detect(meta$celltype,"^Fibroblast")]  <- "Fb"
meta$celltype3[str_detect(meta$celltype,"^EC")]          <- "EC"
meta$celltype3[str_detect(meta$celltype,"^NK")]          <- "NK"
meta$celltype3[str_detect(meta$celltype,"^Bcell")]       <- "Bcell"
combined@meta.data <- meta

# ---- 保存 ----
save(combined, file="./taa_combined.RData")
message("Script 03 完成：单细胞整合完毕，细胞数 =", ncol(combined))