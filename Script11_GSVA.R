# ============================================================
# Script 11: GSVA/ssGSEA 通路富集分析
# ============================================================
source("Script00_Setup.R")

library(GSVA)
library(msigdbr)
library(pheatmap)
library(clusterProfiler)

setwd("/data/home/caoy/Xinjiang_AAD/fig2/20231204/")
load("../wd.RData")

DefaultAssay(wd) <- "SCT"
Idents(wd)       <- "new_seurat_clusters"

# ==================================================
# 第一步：提取各cluster平均表达矩阵
# ==================================================
expr <- AverageExpression(wd, assays="SCT", slot="data")[[1]]
expr <- expr[rowSums(expr) > 0, ]  # 过滤全0基因
expr <- as.matrix(expr)
message("基因数:", nrow(expr), "  Cluster数:", ncol(expr))

# ==================================================
# 第二步：准备基因集
# ==================================================
# KEGG 基因集
genesets_kegg <- msigdbr(species="Homo sapiens", category="C2")
kegg_list <- subset(genesets_kegg, gs_subcat=="CP:KEGG",
                    select=c("gs_description","human_gene_symbol")) %>%
  as.data.frame()
kegg_list2 <- split(kegg_list$human_gene_symbol,
                    kegg_list$gs_description)

# GO 基因集（BP）
genesets_go <- msigdbr(species="Homo sapiens", category="C5")
go_list <- subset(genesets_go, gs_subcat=="GO:BP",
                  select=c("gs_name","gene_symbol"))
go_list2 <- split(go_list$gene_symbol, go_list$gs_name)

# ==================================================
# 第三步：运行 ssGSEA
# ==================================================
message("运行KEGG ssGSEA...")
gsva_kegg <- gsva(expr, kegg_list2, method="ssgsea")
message("运行GO ssGSEA...")
gsva_go   <- gsva(expr, go_list2, method="ssgsea")

# 保存结果
gsva_kegg_df <- data.frame(Genesets=rownames(gsva_kegg), gsva_kegg, check.names=FALSE)
gsva_go_df   <- data.frame(Genesets=rownames(gsva_go),   gsva_go,   check.names=FALSE)

write.csv(gsva_kegg_df, file="./GSVA/ad_gsva_kegg_res_sct.csv", row.names=FALSE)
write.csv(gsva_go_df,   file="./GSVA/ad_gsva_go_res_sct.csv",   row.names=FALSE)

# ==================================================
# 第四步：热图可视化
# ==================================================
pdf("./GSVA/ad_gsva_kegg_pheatmap.pdf", width=5, height=15)
pheatmap(gsva_kegg,
         show_colnames=TRUE, scale="row",
         color=colorRampPalette(rev(brewer.pal(7,"RdBu")))(60),
         show_rownames=TRUE, fontsize_row=3.5)
dev.off()

# 读取手动筛选的GO BP通路
gsva_go_bp <- read.csv(file="./GSVA/ad_gsva_go_bp_sct.csv",
                       header=TRUE, row.names=1)
pdf("./GSVA/ad_gsva_go_bp_pheatmap.pdf", width=5, height=10)
pheatmap(gsva_go_bp,
         show_colnames=TRUE, scale="row",
         color=colorRampPalette(rev(brewer.pal(7,"RdBu")))(60),
         show_rownames=TRUE, fontsize_row=3.5)
dev.off()

# ==================================================
# 第五步：fGSEA 分析（按分组）
# ==================================================
library(fgsea)
library(org.Hs.eg.db)

Idents(wd) <- "group"
Intima <- subset(wd, loc=="Intima")
Idents(Intima) <- "group"

Intima_markers <- FindAllMarkers(Intima, test.use="wilcox",
                                 only.pos=FALSE, min.pct=0.25,
                                 logfc.threshold=0.25)

# ENTREZ ID 转换
data_name <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = Intima_markers$gene,
  keytype = "SYMBOL",
  column  = "ENTREZID"
)
data_name$gene <- data_name$SYMBOL
inter <- merge(Intima_markers, data_name, all=FALSE)
Intima_markers <- inter

# 按组进行fGSEA
group_markers <- split(Intima_markers, Intima_markers$cluster)

fgsea_results <- lapply(group_markers, function(deg){
  deg <- deg[order(deg$avg_log2FC, decreasing=TRUE), ]
  deg <- deg[!duplicated(deg$gene), ]
  genelist <- structure(deg$avg_log2FC, names=deg$gene)
  fgsea(go_list2, stats=genelist, nperm=1000)
})

for(grp in names(fgsea_results)){
  write.csv(fgsea_results[[grp]],
            file=paste0("./Intima_", grp, "_go.csv"), row.names=FALSE)
}

message("Script 11 完成：GSVA/fGSEA分析完毕")