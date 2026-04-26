# ============================================================
# Script 02: 去线粒体基因 → SCTransform → PCA → Harmony → 聚类
# ============================================================
source("Script00_Setup.R")
load("./rawdata/ad_qc.RData")
setwd("/data/home/caoy/Xinjiang_AAD/")

# ---- 去除线粒体基因和血红蛋白基因 ----
all_counts <- as.matrix(GetAssayData(object=ad, slot="counts"))

mt_hb_genes <- c(
  "MT-ATP6","MT-ATP8","MTCH1","MT-CO1","MT-CO2","MT-CO3","MT-CYB",
  "MTHFD2","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5",
  "HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"
)

all_filtered_counts <- all_counts[!rownames(all_counts) %in% mt_hb_genes, ]
message("去除前基因数:", nrow(all_counts), "  去除后:", nrow(all_filtered_counts))

ad2 <- CreateSeuratObject(all_filtered_counts,
                          meta.data = ad@meta.data,
                          assay = "Spatial")
ad2@active.ident <- ad@active.ident
ad2@images       <- ad@images
ad2@version      <- ad@version
ad <- ad2
remove(ad2, all_counts, all_filtered_counts)

# ---- SCTransform 标准化 ----
# 用sample回归批次效应
ad <- SCTransform(ad, assay="Spatial",
                  vars.to.regress="sample",
                  verbose=FALSE)

# ---- 降维 ----
dir.create("pca", showWarnings=FALSE)

ad <- RunPCA(ad, assay="SCT", verbose=FALSE)

# Elbow plot 选择PC数
pdf("./pca/ElbowPlot.pdf", width=6, height=4)
ElbowPlot(ad, ndims=50)
dev.off()

# ---- Harmony 批次整合 ----
ad <- RunHarmony(ad, group.by.vars="sample", reduction="pca")

# ---- 聚类 ----
# 用Harmony降维结果进行聚类
ad <- FindNeighbors(object=ad, reduction="harmony", dims=1:30)

# 尝试多个分辨率
ad <- FindClusters(object=ad,
                   resolution=c(0.1, 0.2, 0.3, 0.4, 0.5))

pdf("./pca/all_harmony_clustree.pdf", width=12, height=8)
clustree(ad)
dev.off()

# 选用分辨率0.1（根据clustree结果）
set.seed(123456)
ad <- FindClusters(object=ad, resolution=0.1)

# ---- 可视化嵌入 ----
ad <- RunTSNE(ad, reduction="harmony", dims=1:30)
ad <- RunUMAP(ad, reduction="harmony", dims=1:30)

# ---- 可视化 ----
pdf("./pca/DimPlot_tsne_cluster.pdf", width=5.5, height=4)
DimPlot(ad, reduction="tsne", label=TRUE, raster=FALSE,
        repel=TRUE, pt.size=0.2, label.size=4) +
  scale_color_brewer(palette="Set2")
dev.off()

pdf("./pca/DimPlot_tsne_sample.pdf", width=6, height=4)
DimPlot(ad, reduction="tsne", group.by="sample",
        raster=FALSE, pt.size=0.2) +
  scale_color_manual(values=colorRampPalette(brewer.pal(7,"YlOrBr"))(17))
dev.off()

pdf("./pca/SpatialDimPlot_all.pdf", width=60, height=4)
SpatialDimPlot(ad, image.alpha=0.8, alpha=c(0.59,1),
               group.by="seurat_clusters",
               label=TRUE, label.size=2,
               ncol=19, pt.size.factor=1.5)
dev.off()

# ---- 重新注释：合并cluster为new_seurat_clusters ----
meta_data <- ad@meta.data
meta_data$new_seurat_clusters <- as.character(meta_data$seurat_clusters)

# 合并规则（根据生物学意义）
merge_rules <- list(
  "c0"=0, "c1"=1, "c2"=2, "c3"=3,
  "c4"=c(4,5), "c5"=c(6,8,14), "c6"=7,
  "c7"=9, "c8"=c(10,12), "c9"=11,
  "c10"=13, "c11"=15
)
for(new_name in names(merge_rules)){
  old_ids <- as.character(merge_rules[[new_name]])
  meta_data$new_seurat_clusters[meta_data$seurat_clusters %in% old_ids] <- new_name
}
meta_data$new_seurat_clusters <- factor(
  meta_data$new_seurat_clusters,
  levels=c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11")
)

ad@meta.data <- meta_data

# ---- 保存 ----
save(ad, file="./fig3/ad_clustered.RData")
message("Script 02 完成：聚类数 =", length(unique(ad@meta.data$seurat_clusters)))