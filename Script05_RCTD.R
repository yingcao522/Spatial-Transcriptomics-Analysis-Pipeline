# ============================================================
# Script 05: RCTD 解卷积（对19个样本逐一进行）
# ============================================================
source("Script00_Setup.R")

library(spacexr)
library(Matrix)

setwd("/data/home/caoy/Xinjiang_AAD/RCTD2/")
load("../taa_combined.RData")
load("../ad_Anchors_sctgroup2.RData")

# ==================================================
# 第一步：构建scRNA参考
# ==================================================
# 过滤，只保留9种主要细胞类型
combined_sub <- subset(
  x = combined,
  celltype3 %in% c("EC","Fb","SMC","MSC","Mp","Bcell","NK","Plasma","Tcell")
)

meta_ref  <- combined_sub@meta.data
cell_types <- as.factor(meta_ref$celltype3)
names(cell_types) <- meta_ref$cells

nUMI_ref  <- meta_ref$nUMI
names(nUMI_ref) <- meta_ref$cells

# 去除线粒体基因
counts_ref <- as.matrix(combined_sub@assays$RNA@counts)
mt_genes   <- grep("^MT-", rownames(counts_ref), value=TRUE)
counts_ref <- counts_ref[!rownames(counts_ref) %in% mt_genes, ]

reference  <- Reference(counts_ref, cell_types, nUMI_ref)
message("参考数据集：", print(table(reference@cell_types)))

# ==================================================
# 第二步：提取各样本坐标位置
# ==================================================
sample_list <- c("s1p2","s2p2","s3p1","s4p3","s5p4","s6p4",
                 "s9p5","s10p6","s11p6","s12p7","s13p7",
                 "s14p8","s15p8","s16p4","s17p4","s18p4","s19p4")

# ==================================================
# 第三步：定义RCTD运行函数
# ==================================================
run_RCTD_sample <- function(sample_name, ad_obj, reference){
  message("\n====== 处理样本：", sample_name, " ======")

  # 提取子集
  ad_sub <- ad_obj[, ad_obj@meta.data$sample == sample_name]

  # 坐标（去掉前3列：tissue/row/col）
  coords        <- dplyr::select(ad_obj@images[[sample_name]]@coordinates, -c(1:3))
  spatialcounts <- as.matrix(ad_sub@assays$Spatial@counts)
  nUMI_sp       <- ad_sub@meta.data$nCount_Spatial
  names(nUMI_sp) <- rownames(ad_sub@meta.data)

  # 创建 SpatialRNA 对象
  puck <- SpatialRNA(coords, spatialcounts, nUMI_sp)
  message("  基因数:", nrow(puck@counts), "  Spot数:", ncol(puck@counts))

  # 运行 RCTD
  myRCTD <- create.RCTD(
    puck, reference,
    max_cores      = 10,
    UMI_min        = 100,
    counts_MIN     = 100,
    class_df       = NULL,
    MAX_MULTI_TYPES = 20,
    keep_reference  = TRUE
  )
  myRCTD <- run.RCTD(myRCTD, doublet_mode='full')

  # 保存RData
  save(myRCTD, file=paste0("./myRCTD_", sample_name, ".RData"))

  # 提取并归一化权重
  results      <- myRCTD@results
  norm_weights <- normalize_weights(results$weights)

  # 保存CSV
  write.csv(results$weights,
            file=paste0("./", sample_name, "_result.csv"))
  write.csv(norm_weights,
            file=paste0("./", sample_name, "_result_normalize.csv"))

  # 绘制每种细胞类型的空间分布
  cell_type_names <- myRCTD@cell_type_info$info[[2]]
  spatialRNA_obj  <- myRCTD@spatialRNA
  resultsdir      <- sample_name
  dir.create(resultsdir, showWarnings=FALSE)

  plot_weights(cell_type_names, spatialRNA_obj, resultsdir, norm_weights)
  plot_weights_unthreshold(cell_type_names, spatialRNA_obj, resultsdir, norm_weights)
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA_obj)

  message("  样本", sample_name, "完成！")
  return(norm_weights)
}

# ==================================================
# 第四步：循环处理所有样本
# ==================================================
all_results <- list()
for(s in sample_list){
  result <- tryCatch(
    run_RCTD_sample(s, ad, reference),
    error = function(e){
      message("样本", s, "出错:", e$message)
      return(NULL)
    }
  )
  if(!is.null(result)) all_results[[s]] <- result
}

# ==================================================
# 第五步：合并所有样本结果
# ==================================================
norm_weights_all <- do.call(rbind, all_results)
write.csv(norm_weights_all, file="./all_result_normalize.csv")

message("Script 05 完成：RCTD解卷积，共", nrow(norm_weights_all), "个spots")