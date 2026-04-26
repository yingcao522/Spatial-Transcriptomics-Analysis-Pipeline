# ============================================================
# Script 09: 以RCTD细胞类型比例为特征，重建Seurat对象
#            目的：基于细胞组成相似性对spot进行新的聚类
# ============================================================
source("Script00_Setup.R")
setwd("/data/home/caoy/Xinjiang_AAD/fig1/")

load("../fig3/ad_0722.RData")

# ==================================================
# 第一步：提取RCTD比例矩阵，转置后作为"基因表达"
# ==================================================
meta.data  <- ad@meta.data
RCTD_cols  <- c("EC","Fb","MSC","SMC","Mp","NK","Tcell","Bcell","Plasma")

# 提取比例并转换为数值
meta.data2 <- as.data.frame(lapply(meta.data[, RCTD_cols], as.numeric))
rownames(meta.data2) <- rownames(meta.data)

# 转置：基因(细胞类型) × spots
expr_matrix  <- as.matrix(meta.data2)
expr_matrix  <- t(expr_matrix)  # 9 × n_spots

# ==================================================
# 第二步：创建新的Seurat对象
# ==================================================
st3 <- CreateSeuratObject(
  counts    = expr_matrix,
  meta.data = ad@meta.data,
  assay     = "RNA"
)
st3@images <- ad@images

# ==================================================
# 第三步：标准化与聚类
# ==================================================
st3 <- NormalizeData(st3, normalization.method="LogNormalize",
                     scale.factor=10000)
st3 <- FindVariableFeatures(st3, selection.method="vst", nfeatures=9)
all.genes <- rownames(st3)
st3 <- ScaleData(st3, features=all.genes)
st3 <- RunPCA(st3, assay="RNA", verbose=FALSE, approx=FALSE)

# 分辨率选择
st3 <- FindNeighbors(object=st3, reduction="pca", dims=1:8)
st3 <- FindClusters(object=st3,
                    resolution=c(0.01,0.02,0.03,0.05,0.1,0.2,0.3))

pdf("./st3_clustree.pdf", width=10, height=8)
clustree(st3)
dev.off()

# 使用分辨率0.1
st3 <- FindClusters(object=st3, resolution=0.1)

# 可视化嵌入
set.seed(12598326)
st3 <- RunTSNE(st3, reduction="pca",
               check_duplicates=FALSE, dims=1:8)
st3 <- RunUMAP(st3, reduction="pca", dims=1:8)

# ==================================================
# 第四步：可视化
# ==================================================
pdf("./ad_cell_dimplot_tsne.pdf", width=4.3, height=3.4)
DimPlot(object=st3, reduction='tsne',
        label=TRUE, group.by="seurat_clusters",
        raster=FALSE, repel=TRUE,
        cols=c("#FFCCCC99","#CCFF0099","#FFCC0099","#00FF0099","#6699FF99","#CC33FF99"),
        pt.size=0.2, label.size=4)
dev.off()

pdf("./ad_cell_dimplot_maxcell.pdf", width=4.7, height=3.4)
DimPlot(st3, reduction="tsne", group.by="max_cell",
        label=TRUE, raster=FALSE, repel=TRUE,
        pt.size=0.2, label.size=4) +
  scale_color_brewer(palette="Set3")
dev.off()

pdf("./ad_cell_dimplot_loc.pdf", width=4.7, height=3.4)
DimPlot(st3, reduction="tsne", group.by="loc",
        label=TRUE, raster=FALSE, repel=TRUE,
        pt.size=0.2, label.size=4)
dev.off()

# 空间分布
pdf("./ad_cells_SpatialDimPlot.pdf", width=60, height=3)
SpatialDimPlot(st3, image.alpha=0.3, alpha=c(0.59,1),
               group.by="seurat_clusters",
               label=TRUE, label.size=2,
               ncol=19, pt.size.factor=2)
dev.off()

# cluster与位置的统计
ad_cells <- FetchData(st3, vars=c("seurat_clusters","loc")) %>%
  dplyr::count(seurat_clusters, loc) %>%
  tidyr::spread(seurat_clusters, n)

ad_cells <- FetchData(st3, vars=c("seurat_clusters","max_cell")) %>%
  dplyr::count(seurat_clusters, max_cell) %>%
  tidyr::spread(seurat_clusters, n)

# ==================================================
# 第五步：合并wd对象的st3聚类结果回ad
# ==================================================
save(st3, file="./ad_cell.rdata")
message("Script 09 完成：st3对象聚类完毕，cluster数 =",
        length(unique(st3@meta.data$seurat_clusters)))