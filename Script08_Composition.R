# ============================================================
# Script 08: 细胞类型组成统计、百分比图、显著性检验
# ============================================================
source("Script00_Setup.R")

library(ggsignif)
library(ggpubr)
library(reshape2)
library(plyr)

setwd("/data/home/caoy/Xinjiang_AAD/fig3/")
load("./ad_0722.RData")

# ==================================================
# 第一部分：提取宽格式数据转长格式
# ==================================================
metadata   <- ad@meta.data
deconProp  <- dplyr::select(metadata, c(19:27))
deconProp2 <- as.data.frame(lapply(deconProp, as.numeric))
deconProp2$cells <- rownames(deconProp)
rownames(deconProp2) <- rownames(deconProp)

# 宽→长格式
deconProp3 <- melt(deconProp2,
                   id.vars       = "cells",
                   variable.name = "Cell",
                   value.name    = "Nor_pro")

# 添加分组信息
for(sname in unique(metadata$sample)){
  deconProp3$sample[str_detect(deconProp3$cells, paste0("^",sname))] <- sname
  deconProp3$group[str_detect(deconProp3$cells, paste0("^",sname))]  <-
    metadata$group[metadata$sample==sname][1]
  deconProp3$patient[str_detect(deconProp3$cells, paste0("^",sname))] <-
    metadata$patient[metadata$sample==sname][1]
}

# Factor 设置
deconProp3$Cell <- factor(deconProp3$Cell,
                          levels=c("EC","Fb","MSC","SMC","Mp","NK","Tcell","Bcell","Plasma"))
deconProp3$group <- factor(deconProp3$group,
                           levels=c("wd1","wd2","wd3","tbg","zsgx","zjz"))
deconProp3$sample <- factor(deconProp3$sample,
                             levels=c("s1p2","s2p2","s3p1","s4p3","s5p4","s6p4",
                                      "s9p5","s10p6","s11p6","s12p7","s13p7",
                                      "s14p8","s15p8","s16p4","s17p4","s18p4","s19p4"))

# ==================================================
# 第二部分：各样本细胞比例堆积图
# ==================================================
col1 <- c("#FB8072","#FFFFB3","#BEBADA","#80B1D3","#8DD3C7",
          "#D9D9D9","#B3DE69","#FCCDE5","#FDB462")

pdf("./seurat_cell/cells_Prop_all_samples.pdf", width=18, height=3)
ggplot(deconProp3, aes(x=cells, y=Nor_pro, fill=Cell)) +
  scale_fill_manual(values=col1) +
  geom_bar(stat="identity", position='fill') +
  facet_grid(. ~ sample, scales="free") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

# ==================================================
# 第三部分：各分组间细胞比例箱线图+显著性
# ==================================================
wd_data <- subset(deconProp3, group %in% c("wd1","wd2","wd3"))
wd_data$group <- factor(wd_data$group, levels=c("wd1","wd2","wd3"))

compaired <- list(c("wd1","wd2"), c("wd1","wd3"), c("wd2","wd3"))

pdf("./seurat_cell/cell_signif_wd_groups.pdf", width=8, height=3)
ggplot(wd_data, aes(x=group, y=Nor_pro, color=group, fill=group)) +
  geom_boxplot(outlier.shape=NA, width=0.6) +
  facet_grid(. ~ Cell) +
  ylim(0, 1.3) +
  geom_signif(comparisons=compaired,
              test="wilcox.test",
              step_increase=0.1,
              textsize=2,
              map_signif_level=FALSE) +
  theme_classic()
dev.off()

# ==================================================
# 第四部分：按cluster统计细胞组成（百分比堆积）
# ==================================================
merge_data <- metadata[, c("new_seurat_clusters")]
deconProp2_merge <- cbind(deconProp2, cluster=metadata$new_seurat_clusters)

d02 <- melt(deconProp2_merge,
            id.vars       = c("cells","cluster"),
            variable.name = "Cell",
            value.name    = "Normalised_probabilities")

d02$cluster <- factor(d02$cluster,
                      levels=c("c0","c1","c2","c3","c4","c5",
                               "c6","c7","c8","c9","c10","c11"))
d02$Cell    <- factor(d02$Cell,
                      levels=c("EC","Fb","MSC","SMC","Mp","NK","Tcell","Bcell","Plasma"))

pdf("./seurat_cell/cluster_cell_composition.pdf", width=10, height=2.8)
ggplot(data=d02) +
  geom_boxplot(mapping=aes(x=cluster, y=Normalised_probabilities, colour=Cell),
               alpha=0.5, size=0.3, outlier.shape=NA, width=0.6) +
  theme_bw(base_line_size=0.5) +
  scale_color_igv()
dev.off()

# ==================================================
# 第五部分：按cluster统计各细胞类型比例（堆积柱状图）
# ==================================================
cell_count_list <- list()
for(cl in levels(d02$cluster)){
  sub  <- d02[d02$cluster==cl, ]
  agg  <- aggregate(sub$Normalised_probabilities,
                    by=list(type=sub$Cell), sum)
  agg$group <- cl
  cell_count_list[[cl]] <- agg
}
cell_count <- do.call(rbind, cell_count_list)

data3 <- ddply(cell_count, 'group', transform,
               percent_x = x/sum(x)*100)
data3$group <- factor(data3$group,
                      levels=c("c0","c1","c2","c3","c4","c5",
                               "c6","c7","c8","c9","c10","c11"))

pdf("./cell_cluster_Prop.pdf", width=4, height=3)
ggplot(data3, aes(x=group, y=percent_x, fill=type)) +
  scale_fill_brewer(palette="Set1") +
  geom_bar(stat="identity", position='fill', width=0.5) +
  theme_classic()
dev.off()

write.csv(data3, file="./fig2e.csv")
message("Script 08 完成：细胞组成分析完毕")