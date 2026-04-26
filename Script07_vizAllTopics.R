# ============================================================
# Script 07: vizAllTopics 可视化 RCTD 解卷积结果
# 每个spot显示为饼图，颜色按细胞类型比例
# ============================================================
source("Script00_Setup.R")

library(STdeconvolve)  # 提供 vizAllTopics 函数

setwd("/data/home/caoy/Xinjiang_AAD/fig3/")
load("./ad_0722.RData")

# ---- 定义可视化函数 ----
# 参数说明：
#   - deconProp: RCTD比例矩阵 (spots × cell_types)
#   - coords:    空间坐标 (spots × 2)，列名需改为 x, y
#   - annot:     分群标签 (spots，用于外圈颜色)
#   - r:         饼图半径
vizAllTopics_sample <- function(seurat_sub, loc_data, output_pdf,
                                 meta_col_range = 19:27,
                                 width=5, height=6.4, r=75){

  # 提取坐标
  coords <- dplyr::select(loc_data, -c(1:3))
  names(coords)[names(coords)=="imagerow"] <- "x"
  names(coords)[names(coords)=="imagecol"] <- "y"

  # 提取元数据
  metadata <- seurat_sub@meta.data

  # 分群标签（外圈颜色）
  annot <- metadata$new_seurat_clusters
  names(annot) <- metadata$cells

  # RCTD 比例（第19-27列 = EC/Fb/MSC/SMC/Mp/NK/Tcell/Bcell/Plasma）
  deconProp  <- metadata[, meta_col_range]
  deconProp2 <- as.data.frame(lapply(deconProp, as.numeric))
  rownames(deconProp2) <- rownames(deconProp)

  # 绘图
  pdf(output_pdf, width=width, height=height)
  vizAllTopics(
    deconProp2, coords,
    groups     = annot,
    lwd        = 0.00006,
    group_cols = rainbow(length(levels(annot))),
    r          = r
  )
  dev.off()
  message("已保存：", output_pdf)
}

# ---- 提取每个样本的位置信息 ----
sample_locs <- list(
  s1p2  = ad@images[["s1p2"]]@coordinates,
  s2p2  = ad@images[["s2p2"]]@coordinates,
  s3p1  = ad@images[["s3p1"]]@coordinates,
  s4p3  = ad@images[["s4p3"]]@coordinates,
  s5p4  = ad@images[["s5p4"]]@coordinates,
  s6p4  = ad@images[["s6p4"]]@coordinates,
  s9p5  = ad@images[["s9p5"]]@coordinates,
  s10p6 = ad@images[["s10p6"]]@coordinates,
  s11p6 = ad@images[["s11p6"]]@coordinates,
  s12p7 = ad@images[["s12p7"]]@coordinates,
  s13p7 = ad@images[["s13p7"]]@coordinates,
  s14p8 = ad@images[["s14p8"]]@coordinates,
  s15p8 = ad@images[["s15p8"]]@coordinates,
  s16p4 = ad@images[["s16p4"]]@coordinates,
  s17p4 = ad@images[["s17p4"]]@coordinates,
  s18p4 = ad@images[["s18p4"]]@coordinates,
  s19p4 = ad@images[["s19p4"]]@coordinates
)

# ---- 批量绘图 ----
sample_names <- names(sample_locs)

for(sname in sample_names){
  # 提取子集
  ad_sub  <- subset(x=ad, sample %in% sname)

  # 修正images
  img_idx <- which(names(ad@images) == sname)
  ad_sub@images <- ad@images[img_idx]

  # 运行
  vizAllTopics_sample(
    seurat_sub   = ad_sub,
    loc_data     = sample_locs[[sname]],
    output_pdf   = paste0("./", sname, "_vizAllTopics_sctgroup2.pdf"),
    meta_col_range = 19:27,   # EC~Plasma在meta.data中的列范围
    r = 75
  )
}

# ---- 空间特征图：各细胞类型比例 ----
DefaultAssay(ad) <- "SCT"

pdf("./cell_SpatialDimPlot.pdf", width=52, height=36)
SpatialFeaturePlot(ad,
                   features = c("EC","Fb","MSC","SMC",
                                "Mp","NK","Tcell","Bcell","Plasma"),
                   ncol=17)
dev.off()

# ---- Marker基因空间分布 ----
marker <- c(
  "PECAM1","VWF",      # EC
  "ACTA2","DCN","LOX", # Fb
  "THY1","MYH11",      # MSC
  "CALD1",             # SMC
  "CD163","CD68",      # Mp
  "NKG7","KLRD1",      # NK
  "CD3D","CD3G",       # T
  "HLA-DRA","CD79A","MZB1" # B/Plasma
)

pdf("./cell_marker_SpatialDimPlot.pdf", width=52, height=60)
SpatialFeaturePlot(ad, features=marker, ncol=17)
dev.off()

message("Script 07 完成：vizAllTopics 绘图完毕")