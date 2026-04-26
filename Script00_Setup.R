# ============================================================
# Script 00: 环境设置、包加载、全局参数
# ============================================================

# ---- 包加载 ----
library(ggplot2)
library(Seurat)
library(dplyr)
library(hdf5r)
library(cowplot)
library(ggsci)
library(gridExtra)
library(tidyverse)
library(patchwork)
library(harmony)
library(clustree)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)

# ---- 工作目录 ----
setwd("/data/home/caoy/Xinjiang_AAD/")

# ---- 全局配色方案 ----
cb_palette <- c(
  "#8F7700FF","#D198D6","#D43F3AFF","#EEA236FF","#5CB85CFF",
  "#46B8DAFF","#9632B8FF",
  "#FFDC91FF","#ADB6B6FF","#20854EFF","#7876B1FF","#6F99ADFF",
  "#E18727FF","#EE4C97FF",
  "#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF",
  "#0073C2FF","#EFC000FF","#59B7D2","#CD534CFF","#FF7F0EFF"
)

my36colors <- c(
  '#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87',
  '#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658',
  '#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4C755','#F7F398',
  '#AA9A59','#E63863','#E39A35','#C1E6F3','#6778AE','#91D0BE','#B53E2B',
  '#712820','#DCC1DD','#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963','#968175'
)

# ---- 样本信息表 ----
# 样本名 → 分组 → 患者 的映射关系
sample_info <- data.frame(
  sample  = c("s1p2","s2p2","s3p1","s4p3","s5p4","s6p4",
              "s9p5","s10p6","s11p6",
              "s12p7","s13p7","s14p8","s15p8",
              "s16p4","s17p4","s18p4","s19p4"),
  group   = c("wd1","wd1","wd1","wd2","wd2","wd2",
              "wd3","wd3","wd3",
              "tbg","tbg","tbg","tbg",
              "zsgx","zsgx","zjz","zjz"),
  patient = c("p2","p2","p1","p3","p4","p4",
              "p5","p6","p6",
              "p7","p7","p8","p8",
              "p4","p4","p4","p4"),
  stringsAsFactors = FALSE
)

# ---- 细胞类型颜色 ----
celltype_colors <- c(
  "EC"     = "#FF00CC",
  "Fb"     = "#FFFF66",
  "MSC"    = "#ADD26D",
  "SMC"    = "#336633",
  "Mp"     = "#FF7F0EFF",
  "NK"     = "#951E23",
  "Tcell"  = "#A25B2B",
  "Bcell"  = "#741E72",
  "Plasma" = "#7876B1FF"
)

message("Script 00 完成：环境设置完毕")