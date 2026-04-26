# ============================================================
# Script 14: Mfuzz 模糊聚类 - AAD进展中的基因趋势
# ============================================================
source("Script00_Setup.R")

library(Mfuzz)
library(clusterProfiler)
library(GOplot)
library(GOSemSim)

setwd("/data/home/caoy/Xinjiang_AAD/fig3/")
load("./ad_Anchors_sctgroup2.RData")

# ==================================================
# 第一步：计算各样本平均表达量（三个代表性样本）
# ==================================================
wd_subset <- subset(ad, sample %in% c("s1p2","s4p3","s10p6"))

ave <- AverageExpression(wd_subset, group.by='sample', assays='SCT')
av  <- as.matrix(ave$SCT)

write.csv(av, file="./wd_AverageExpression.csv")
av <- read.csv("./wd_AverageExpression.csv", header=TRUE, row.names=1)
av <- as.matrix(av)

# ==================================================
# 第二步：Mfuzz 模糊聚类
# ==================================================
eset <- new('ExpressionSet', exprs=av)

# 预处理
set  <- filter.NA(eset, thres=0.25)   # 过滤缺失基因
set  <- fill.NA(set, mode='knn')       # 均值填充
set  <- filter.std(set, min.std=0)     # 过滤低方差
set  <- standardise(set)               # Z-score标准化

# 最优模糊参数
m_est <- mestimate(set)
message("最优模糊参数 m =", m_est)

# 聚类（c=5个趋势）
set.seed(123)
cluster_num  <- 5
mfuzz_cluster <- mfuzz(set, c=cluster_num, m=m_est)

# 可视化
pdf("./mfuzz_5clusters.pdf", width=10, height=10)
mfuzz.plot(set, cl=mfuzz_cluster, mfrow=c(2,3),
           colo=brewer.pal(11,"PrGn"),
           time.labels=c('s1p2','s4p3','s10p6'))
dev.off()

# ==================================================
# 第三步：提取各趋势的基因
# ==================================================
# 计算各cluster的平均表达趋势
cluster_assignments <- mfuzz_cluster$cluster
gene_cluster_df <- data.frame(
  gene    = names(cluster_assignments),
  cluster = cluster_assignments
)

# 手动计算趋势（基于平均表达量的变化方向）
avg_exp_long <- read.csv("./wd_AverageExpression.csv",
                          header=TRUE, row.names=1)

trend_data <- data.frame(
  gene = rownames(avg_exp_long),
  wd1  = avg_exp_long$s1p2,
  wd2  = avg_exp_long$s4p3,
  wd3  = avg_exp_long$s10p6
)

trend_data$trend_wd2_vs_wd1 <- ifelse(trend_data$wd2 > trend_data$wd1,
                                       "Up", "Down")
trend_data$trend_wd3_vs_wd1 <- ifelse(trend_data$wd3 > trend_data$wd1,
                                       "Up", "Down")

# 提取各方向的基因
genes_wd1_high <- trend_data$gene[
  trend_data$trend_wd2_vs_wd1=="Down" & trend_data$trend_wd3_vs_wd1=="Down"]
genes_wd3_high <- trend_data$gene[
  trend_data$trend_wd2_vs_wd1=="Up"   & trend_data$trend_wd3_vs_wd1=="Up"]

message("wd1特高基因:", length(genes_wd1_high))
message("wd3特高基因:", length(genes_wd3_high))

# ==================================================
# 第四步：GO富集分析
# ==================================================
deg_list_mfuzz <- list(
  wd1_high = genes_wd1_high,
  wd3_high = genes_wd3_high
)

go_results <- lapply(deg_list_mfuzz, function(genes){
  enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
})

for(gname in names(go_results)){
  write.csv(as.data.frame(go_results[[gname]]),
            file=paste0("./go_bp_", gname, ".csv"),
            row.names=FALSE)
}

# ==================================================
# 第五步：GOBubble 可视化
# ==================================================
# 准备 genelist（带logFC）
for(gname in names(deg_list_mfuzz)){
  genes <- deg_list_mfuzz[[gname]]
  gl <- data.frame(
    ID    = genes,
    logFC = ifelse(gname=="wd3_high", 1, -1)  # 方向标记
  )
  write.csv(gl, file=paste0("./genelist_", gname, ".csv"),
            row.names=FALSE)
}

# GOBubble（以wd3_high为例）
go_df <- as.data.frame(go_results[["wd3_high"]])
if(nrow(go_df) > 0){
  genelist <- read.csv("./genelist_wd3_high.csv", header=TRUE)
  circ     <- circle_dat(go_df[1:30,], genelist)

  pdf("./GOBubble_wd3_high.pdf", width=6, height=4)
  GOBubble(circ, display="single", title="wd3特异性高表达基因通路",
           table.legend=TRUE, labels=5)
  dev.off()
}

message("Script 14 完成：Mfuzz趋势分析完毕")