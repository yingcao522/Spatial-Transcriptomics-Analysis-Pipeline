# ============================================================
# Script 13: MistyR 空间相互作用分析
#   - Level 1: 细胞类型共定位 (使用st3对象的RCTD比例)
#   - Level 2: 基因共表达   (使用wd对象的SCT数据)
# ============================================================
source("Script00_Setup.R")

library(mistyR)
library(tidyverse)

future::plan(future::multisession)
setwd("/data/home/caoy/Xinjiang_AAD/fig6/")

source("/data/home/caoy/flu_st_result_exon/st_all_result_0313/fig2/analysis/utils/misty_utilities.R")

# ==================================================
# 通用分析函数
# ==================================================
run_colocalization <- function(slide, assay, useful_features,
                               out_label, misty_out_alias="result"){
  view_assays  <- list("main"="intra", "juxta"="juxta", "para"="para")
  view_types   <- list("main"="intra", "juxta"="juxta", "para"="para")
  view_features <- list("main"=useful_features,
                        "juxta"=useful_features, "para"=useful_features)
  view_params  <- list("main"=NULL, "juxta"=2, "para"=5)

  misty_out <- paste0(misty_out_alias, out_label, "_", assay)

  run_misty_seurat(
    visium.slide  = slide,
    view.assays   = list("main"=assay,"juxta"=assay,"para"=assay),
    view.features = view_features,
    view.types    = view_types,
    view.params   = view_params,
    spot.ids      = NULL,
    out.alias     = misty_out
  )
  return(misty_out)
}

# ==================================================
# Level 1: 细胞类型共定位（st3对象）
# ==================================================
load("../fig1/ad_cell.rdata")  # st3对象

wd1_st3 <- subset(x=st3, sample %in% c("s1p2","s2p2","s3p1"))
wd2_st3 <- subset(x=st3, sample %in% c("s4p3","s5p4","s6p4"))
wd3_st3 <- subset(x=st3, sample %in% c("s9p5","s10p6","s11p6"))

# 修正images
wd1_st3@images <- st3@images[1:3]
wd2_st3@images <- st3@images[4:6]
wd3_st3@images <- st3@images[7:9]

# 运行
mout_ct_wd1 <- run_colocalization(wd1_st3, "RNA",
                                   rownames(wd1_st3), NULL, "celltype_wd1")
mout_ct_wd2 <- run_colocalization(wd2_st3, "RNA",
                                   rownames(wd2_st3), NULL, "celltype_wd2")
mout_ct_wd3 <- run_colocalization(wd3_st3, "RNA",
                                   rownames(wd3_st3), NULL, "celltype_wd3")

save(mout_ct_wd1, file="./mout_celltype_wd1.RData")
save(mout_ct_wd2, file="./mout_celltype_wd2.RData")
save(mout_ct_wd3, file="./mout_celltype_wd3.RData")

# ==================================================
# Level 2: 基因共表达（wd对象，使用空间基因集）
# ==================================================
load("../fig3/wd.RData")
gene_sets <- read.csv("/data/home/caoy/Xinjiang_AAD/fig3/20231225/wd123.1sample.marker.deg.top50.csv",
                      header=TRUE, sep=",")

wd1 <- subset(wd, sample %in% c("s1p2","s2p2","s3p1"))
wd2 <- subset(wd, sample %in% c("s4p3","s5p4","s6p4"))
wd3 <- subset(wd, sample %in% c("s9p5","s10p6","s11p6"))

wd1@images <- wd@images[1:3]
wd2@images <- wd@images[4:6]
wd3@images <- wd@images[7:9]

mout_gene_wd1 <- run_colocalization(wd1, "Spatial",
                                     gene_sets$wd1, NULL, "gene_wd1")
mout_gene_wd2 <- run_colocalization(wd2, "Spatial",
                                     gene_sets$wd2, NULL, "gene_wd2")
mout_gene_wd3 <- run_colocalization(wd3, "Spatial",
                                     gene_sets$wd3, NULL, "gene_wd3")

save(mout_gene_wd1, file="./mout_gene_wd1.RData")
save(mout_gene_wd2, file="./mout_gene_wd2.RData")
save(mout_gene_wd3, file="./mout_gene_wd3.RData")

# ==================================================
# 可视化函数
# ==================================================
plot_misty_results <- function(mout_path, output_prefix, view="intra"){
  misty_res <- collect_results(mout_path)
  a <- misty_res[["importances.aggregated"]]

  df_importance <- subset(a, view==view) %>%
    dplyr::select(c(Predictor, Target, Importance)) %>%
    pivot_wider(names_from=Target, values_from=Importance, values_fill=0) %>%
    column_to_rownames('Predictor')

  col_fun <- colorRamp2(
    seq(0, 2.5, length=9),
    RColorBrewer::brewer.pal(name="BuGn", n=9)
  )

  ht <- Heatmap(as.matrix(df_importance),
                name="Importance",
                cluster_columns=TRUE, cluster_rows=TRUE,
                rect_gp=gpar(col="black", lwd=0.2),
                col=col_fun, na_col="white",
                row_names_side="left")

  pdf(paste0("./", output_prefix, "_misty_heatmap.pdf"),
      height=3.84, width=4.5)
  draw(ht)
  dev.off()
}

# 批量绘图
plot_misty_results(mout_ct_wd1, "celltype_wd1")
plot_misty_results(mout_ct_wd2, "celltype_wd2")
plot_misty_results(mout_ct_wd3, "celltype_wd3")
plot_misty_results(mout_gene_wd1, "gene_wd1")
plot_misty_results(mout_gene_wd2, "gene_wd2")
plot_misty_results(mout_gene_wd3, "gene_wd3")

message("Script 13 完成：MistyR空间分析完毕")