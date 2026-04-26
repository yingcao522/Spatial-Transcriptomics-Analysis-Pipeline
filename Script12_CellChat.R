# ============================================================
# Script 12: CellChat 细胞通信分析
#   - 模式1：按分组（wd1/wd2/wd3），非空间模式
#   - 模式2：按样本，空间模式（考虑距离）
# ============================================================
source("Script00_Setup.R")

library(CellChat)
library(NMF)
library(ggalluvial)

options(stringsAsFactors=FALSE)
setwd("/data/home/caoy/Xinjiang_AAD/fig4/20231204/")
load("/data/home/caoy/Xinjiang_AAD/fig3/wd.RData")

# 只保留主要非免疫细胞类型
wd <- subset(wd, max_cell %in% c("EC","Fb","SMC","MSC","Mp"))

# 按wd分组拆分
wd1 <- subset(wd, group=="wd1")
wd2 <- subset(wd, group=="wd2")
wd3 <- subset(wd, group=="wd3")

# 修正images
for(grp in c("wd1","wd2","wd3")){
  img_names <- switch(grp,
    "wd1" = c("s1p2","s2p2","s3p1"),
    "wd2" = c("s4p3","s5p4","s6p4"),
    "wd3" = c("s9p5","s10p6","s11p6")
  )
  obj_name <- get(grp)
  obj_name@images <- wd@images[img_names]
  assign(grp, obj_name)
}

# ==================================================
# 通用CellChat运行函数（非空间模式）
# ==================================================
run_cellchat_group <- function(seurat_obj, group_name, output_dir){
  message("处理:", group_name)

  # SCTransform
  seurat_obj <- SCTransform(seurat_obj, assay="Spatial", verbose=FALSE)

  # 创建CellChat对象
  cellchat_obj <- createCellChat(object=seurat_obj, group.by="max_cell")

  # 设置数据库
  CellChatDB <- CellChatDB.human
  cellchat_obj@DB <- CellChatDB

  # 预处理
  cellchat_obj <- subsetData(cellchat_obj)
  future::plan("multisession", workers=10)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- projectData(cellchat_obj, PPI.human)

  # 计算通信概率
  cellchat_obj <- computeCommunProb(cellchat_obj,
                                    raw.use=FALSE, population.size=TRUE)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells=10)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name="netP")

  # 保存
  save(cellchat_obj, file=paste0("./", group_name, ".cellchat.rdata"))

  # 可视化
  df.net <- subsetCommunication(cellchat_obj, slot.name="netP")
  write.csv(df.net, paste0("./", group_name, ".cellchat.df.net.csv"))

  groupSize <- as.numeric(table(cellchat_obj@idents))

  pdf(paste0("./", group_name, ".cellchat_Number.pdf"), width=5, height=4)
  netVisual_circle(cellchat_obj@net$count, vertex.weight=groupSize,
                   weight.scale=TRUE, label.edge=FALSE,
                   title.name="Number of interactions")
  dev.off()

  pdf(paste0("./", group_name, ".cellchat_heatmap.pdf"), width=4.5, height=4)
  netVisual_heatmap(cellchat_obj, measure="weight")
  dev.off()

  return(cellchat_obj)
}

# 运行三个分组
wd1.cellchat <- run_cellchat_group(wd1, "wd1", ".")
wd2.cellchat <- run_cellchat_group(wd2, "wd2", ".")
wd3.cellchat <- run_cellchat_group(wd3, "wd3", ".")

# ==================================================
# 空间模式CellChat（单样本，考虑距离）
# ==================================================
# 单独样本：s1p2, s4p3, s10p6

metadata <- wd@meta.data
metadata2 <- unite(metadata, "cluster2",
                   c("max_cell","loc"), sep="_", remove=FALSE)
wd@meta.data <- metadata2

run_spatial_cellchat <- function(sample_name, seurat_obj_full, wd_obj_full,
                                 scalefactors_dir){
  message("空间CellChat:", sample_name)

  sub_obj <- subset(wd_obj_full, sample==sample_name)
  img     <- wd_obj_full@images[[sample_name]]
  sub_obj@images <- wd_obj_full@images[sample_name]

  data.input <- GetAssayData(sub_obj, slot="data", assay="SCT")
  Idents(sub_obj) <- "cluster2"
  meta_sub <- data.frame(labels=Idents(sub_obj),
                         row.names=names(Idents(sub_obj)))

  # 空间坐标
  spatial.locs <- dplyr::select(img@coordinates, -c(1:3))

  # 比例因子
  sf_json <- jsonlite::fromJSON(
    file.path(scalefactors_dir, sample_name, "spatial/scalefactors_json.json")
  )
  scale.factors <- list(
    spot.diameter = sf_json[["spot_diameter_fullres"]],
    spot          = sf_json[["spot_diameter_fullres"]],
    fiducial      = sf_json[["fiducial_diameter_fullres"]],
    hires         = sf_json[["tissue_hires_scalef"]],
    lowres        = sf_json[["tissue_lowres_scalef"]]
  )

  # 创建空间CellChat对象
  cc <- createCellChat(object=data.input, meta=meta_sub,
                       group.by="labels",
                       datatype="spatial",
                       coordinates=spatial.locs,
                       scale.factors=scale.factors)
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  future::plan("multisession", workers=10)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- projectData(cc, PPI.human)

  # 空间距离约束
  cc <- computeCommunProb(cc,
                          type="truncatedMean", trim=0.1,
                          interaction.length=200,
                          scale.distance=0.01, k.min=10,
                          distance.use=TRUE)
  cc <- filterCommunication(cc, min.cells=3)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  save(cc, file=paste0("./", sample_name, ".cellchat.rdata"))

  # 通路可视化 - SPP1
  pdf(paste0("./", sample_name, "_SPP1_circle.pdf"), width=4, height=4.3)
  netVisual_aggregate(cc, signaling="SPP1", layout="circle")
  dev.off()

  pdf(paste0("./", sample_name, "_SPP1_LR.pdf"), width=6, height=4.3)
  plotGeneExpression(cc, signaling="SPP1")
  dev.off()

  return(cc)
}

rawdata_dir <- "/data/home/caoy/Xinjiang_AAD/rawdata"
s1p2.cellchat  <- run_spatial_cellchat("s1p2",  ad, wd, rawdata_dir)
s4p3.cellchat  <- run_spatial_cellchat("s4p3",  ad, wd, rawdata_dir)
s10p6.cellchat <- run_spatial_cellchat("s10p6", ad, wd, rawdata_dir)

message("Script 12 完成：CellChat分析完毕")