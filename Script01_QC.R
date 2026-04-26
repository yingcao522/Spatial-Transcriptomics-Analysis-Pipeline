# ============================================================
# Script 01: 原始数据加载、元数据添加、QC过滤
# ============================================================
source("Script00_Setup.R")
setwd("/data/home/caoy/Xinjiang_AAD/")

# ---- 加载原始Seurat对象 ----
load("./rawdata/ad.RData")

# ---- 添加样本元数据 ----
metadata <- ad@meta.data
metadata$cells <- rownames(metadata)

# 添加 sample 列
for(i in 1:nrow(sample_info)){
  s <- sample_info$sample[i]
  idx <- which(str_detect(metadata$cells, paste0("^", s)))
  metadata$sample[idx]  <- s
  metadata$group[idx]   <- sample_info$group[i]
  metadata$patient[idx] <- sample_info$patient[i]
}

# 添加 group2 列（S1~S11编号）
group2_map <- list(
  "s1p2"="S1","s2p2"="S1","s3p1"="S2","s4p3"="S3",
  "s5p4"="S4","s6p4"="S4","s7p2"="S5","s8p2"="S5",
  "s9p5"="S6","s10p6"="S7","s11p6"="S7","s12p7"="S8",
  "s13p7"="S8","s14p8"="S9","s15p8"="S9","s16p4"="S10",
  "s17p4"="S10","s18p4"="S11","s19p4"="S11"
)
for(s in names(group2_map)){
  metadata$group2[which(str_detect(metadata$cells, paste0("^",s)))] <- group2_map[[s]]
}

# 重命名列
metadata <- metadata %>%
  dplyr::rename(nCount_Spatial = nUMI,
                nFeature_Spatial = nGene)

ad@meta.data <- metadata

# ---- 计算QC指标 ----
ad$log10GenesPerUMI <- log10(ad$nFeature_Spatial) / log10(ad$nCount_Spatial)
ad$mitoPercent       <- PercentageFeatureSet(object = ad, pattern = "^MT-")
ad$mitoRatio         <- ad@meta.data$mitoPercent / 100

# 血红蛋白基因比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m     <- match(HB.genes, rownames(ad@assays$Spatial))
HB.genes <- rownames(ad@assays$Spatial)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
ad[["percent.HB"]] <- PercentageFeatureSet(ad, features = HB.genes)

# ---- QC可视化（过滤前）----
dir.create("raw", showWarnings = FALSE)
metadata <- ad@meta.data

pdf("./raw/nCount_per_sample.pdf", width=20, height=7)
p1 <- VlnPlot(ad, features="nUMI", group.by="sample",
              split.by="sample", cols=cb_palette, pt.size=0.1) + NoLegend()
p2 <- SpatialFeaturePlot(ad, features="nUMI") +
      theme(legend.position="right")
plot_grid(p1, p2, ncol=1)
dev.off()

pdf("./raw/nGene_per_sample.pdf", width=20, height=7)
p1 <- VlnPlot(ad, features="nGene", group.by="sample",
              split.by="sample", cols=cb_palette, pt.size=0.1) + NoLegend()
p2 <- SpatialFeaturePlot(ad, features="nGene") +
      theme(legend.position="right")
plot_grid(p1, p2, ncol=1)
dev.off()

pdf("./raw/all_UMIs_counts.pdf", width=6, height=4)
metadata %>%
  ggplot(aes(color=sample, x=nUMI, fill=sample)) +
  geom_density(alpha=0.4) + scale_x_log10() + theme_classic() +
  ylab("Cell density") + geom_vline(xintercept=500) +
  scale_fill_igv() + scale_color_igv() +
  ggtitle("UMIs/transcripts per spot")
dev.off()

pdf("./raw/all_genes_per_spot.pdf", width=6, height=4)
metadata %>%
  ggplot(aes(color=sample, x=nGene, fill=sample)) +
  geom_density(alpha=0.2) + theme_classic() + scale_x_log10() +
  geom_vline(xintercept=300) + scale_fill_igv() + scale_color_igv() +
  ggtitle("genes detected per spot")
dev.off()

pdf("./raw/all_genes_per_UMI.pdf", width=6, height=4)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color=sample, fill=sample)) +
  geom_density(alpha=0.2) + theme_classic() +
  geom_vline(xintercept=0.8) + scale_fill_igv() + scale_color_igv() +
  ggtitle("genes detected per UMI")
dev.off()

# ---- 质量过滤 ----
ad <- subset(x = ad,
             subset = (nCount_Spatial   >= 200) &
                      (nFeature_Spatial >= 200) &
                      (percent.HB        < 3)   &
                      (mitoPercent       <= 20)  &
                      (log10GenesPerUMI  > 0.70))

# ---- QC可视化（过滤后）----
dir.create("filtered", showWarnings = FALSE)
metadata <- ad@meta.data

pdf("./filtered/all_UMIs_counts_clean.pdf", width=6, height=4)
metadata %>%
  ggplot(aes(color=sample, x=nUMI, fill=sample)) +
  geom_density(alpha=0.4) + scale_x_log10() + theme_classic() +
  ylab("Cell density") + geom_vline(xintercept=500) +
  scale_fill_igv() + scale_color_igv()
dev.off()

pdf("./filtered/all_genes_per_spot_boxplot_clean.pdf", width=5, height=4)
metadata %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot(notch=TRUE, notchwidth=0.3, outlier.size=1) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  scale_fill_igv() + ggtitle("genes detected per spot")
dev.off()

# ---- 保存 ----
save(ad, file="./rawdata/ad_qc.RData")
message("Script 01 完成：QC过滤后 spots =", ncol(ad))