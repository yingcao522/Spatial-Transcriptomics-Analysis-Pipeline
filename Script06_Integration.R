# ============================================================
# Script 06: RCTD结果 → 整合到Seurat对象 → 定义解剖位置分层
# ============================================================
source("Script00_Setup.R")
setwd("/data/home/caoy/Xinjiang_AAD/RCTD2/")

load("../ad_Anchors_sctgroup2.RData")

# ---- 读取RCTD归一化权重 ----
norm_weights <- read.csv("./all_result_normalize.csv",
                         header=TRUE, row.names=1, check.names=FALSE)

# 确保数值类型
norm_weights2 <- as.data.frame(lapply(norm_weights, as.numeric))
rownames(norm_weights2) <- rownames(norm_weights)

# 保留辅助列
norm_weights2$cells    <- rownames(norm_weights2)
norm_weights2$max_cell <- apply(
  norm_weights2[, c("EC","Fb","MSC","SMC","Mp","NK","Tcell","Bcell","Plasma")],
  1, function(x) names(x)[which.max(x)]
)
norm_weights2$sec_cell <- apply(
  norm_weights2[, c("EC","Fb","MSC","SMC","Mp","NK","Tcell","Bcell","Plasma")],
  1, function(x){ o<-order(x,decreasing=TRUE); names(x)[o[2]] }
)

# ---- 合并到ad@meta.data ----
meta_data <- ad@meta.data
meta_data <- dplyr::select(meta_data, -any_of(c("EC","Fb","MSC","SMC","Mp",
                                                  "NK","Tcell","Bcell","Plasma",
                                                  "max_cell","sec_cell")))

meta <- merge(meta_data, norm_weights2, by="cells", all=TRUE)
rownames(meta) <- meta[, 1]
ad@meta.data   <- meta

# ---- 添加 new_seurat_clusters ----
meta_data <- ad@meta.data
meta_data$new_seurat_clusters <- as.character(meta_data$seurat_clusters)

cluster_map <- list(
  "c0"="0","c1"="1","c2"="2","c3"="3",
  "c4"=c("4","5"),"c5"=c("6","8","14"),"c6"="7",
  "c7"="9","c8"=c("10","12"),"c9"="11","c10"="13","c11"="15"
)
for(new_c in names(cluster_map)){
  old_cs <- as.character(cluster_map[[new_c]])
  meta_data$new_seurat_clusters[meta_data$seurat_clusters %in% old_cs] <- new_c
}
meta_data$new_seurat_clusters <- factor(
  meta_data$new_seurat_clusters,
  levels = c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11")
)

# ---- 解剖位置标注（需要手动注释CSV）----
# loc 列：Intima / Media / Adventitia
# 从外部CSV读入
# loc_info <- read.csv("./loc/all_loc.csv", header=TRUE, row.names=1)
# loc_info$cells <- rownames(loc_info)
# meta_data <- merge(meta_data, loc_info, by="cells", all.x=TRUE)

# ---- Factor 设置 ----
meta_data$max_cell <- factor(
  meta_data$max_cell,
  levels = c("EC","Fb","MSC","SMC","Mp","NK","Tcell","Bcell","Plasma")
)
meta_data$sample <- factor(
  meta_data$sample,
  levels = c("s1p2","s2p2","s3p1","s4p3","s5p4","s6p4",
             "s9p5","s10p6","s11p6","s12p7","s13p7",
             "s14p8","s15p8","s16p4","s17p4","s18p4","s19p4")
)

ad@meta.data <- meta_data

# ---- 保存 ----
write.csv(meta_data, file="./ad_meta_data.csv")
save(ad, file="../fig3/ad_0722.RData")
message("Script 06 完成：元数据整合完毕，列数 =", ncol(ad@meta.data))