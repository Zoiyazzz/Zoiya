if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
install.packages("ggrepel")
# 安装（如未安装）
if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  install.packages("BiocManager")  # 如果 BiocManager 也没有就先装这个
BiocManager::install("AnnotationDbi")
if (!requireNamespace("ReactomePA", quietly = TRUE)) BiocManager::install("ReactomePA")
library(ReactomePA)
# 加载
library(AnnotationDbi)
if (!requireNamespace("limma", quietly = TRUE))
  install.packages("BiocManager")  # 如果 BiocManager 也没有就先装这个
BiocManager::install("limma")
# ===== 环境与读取 =====
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(AnnotationDbi)
library(limma)

setwd("/Users/zhongjiayan/Documents/analysis_tool/bulkRNA_seq/Colorectal_Cancer")
counts <- read.delim("./GSE32323_series_matrix.txt", comment.char="#", check.names=FALSE)
counts_raw <- read.delim("GSE32323_series_matrix.txt", comment.char="!", check.names=FALSE)

# 确认第一列是探针 ID（名字通常叫 ID_REF）
head(colnames(counts_raw))  # 看看是否为 "ID_REF", "GSM800742", ...

# 设置 rownames 为探针 ID（第一列）
rownames(counts_raw) <- counts_raw$ID_REF
counts <- counts_raw[ , -1]  # 去掉第一列，只保留表达矩阵部分
head(rownames(counts))  # 应该是 "1007_s_at" 这样的 ID
## 安装注释包（只需一次）
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("hgu133plus2.db")

# 加载并转换探针 ID 为 gene symbol
library(hgu133plus2.db)
library(AnnotationDbi)
# 取交集确保合法 ID
valid_ids <- intersect(rownames(counts), keys(hgu133plus2.db, keytype = "PROBEID"))


# 做注释映射
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = valid_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# 应用到表达矩阵
valid <- !is.na(gene_symbols)
counts_clean <- counts[valid, ]
# 创建一个新的数据框，包含 symbol 注释
counts_clean_annotated <- counts_clean[valid, ]
gene_symbols_filtered <- gene_symbols[valid]

# 去除重复 symbol（只保留第一个）
non_dup <- !duplicated(gene_symbols_filtered)

counts_clean_annotated <- counts_clean_annotated[non_dup, ]
rownames(counts_clean_annotated) <- gene_symbols_filtered[non_dup]

# ===== 样本信息 =====

sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("Normal", 17), rep("Tumor", 17), rep("CellLine", 10))
)
rownames(sample_info)
# 读取表达矩阵（你已经有 counts，行为基因，列为样本）
expr_matrix <- counts_clean_annotated

# 构建样本信息表（确保顺序与列名一致）
group <- factor(c(rep("Normal", 17), rep("Tumor", 17), rep("CellLine", 10)))  # 自定义你的分组
design <- model.matrix(~ group)

# ===== 差异表达分析 =====
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, number=Inf, adjust="BH")

# 查看前几行结果
head(results)
write.csv(as.data.frame(results), "DEG_results.csv")
# ===== 火山图 =====
library(ggplot2)
results$gene <- rownames(results)
png("Volcanot_plot.png")
ggplot(results, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 1)) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-value")
dev.off()
# ===== complexheatmap =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("circlize")  # 配色支持

library(ComplexHeatmap)
library(circlize)
# 挑选你想要展示的 16 个 marker 基因（可以改成你感兴趣的）
marker_genes <- c("ASPP1", "C17ORF91", "CASZ1", "CLEC3B", "DUSP5", "FOXO1",
                  "GAS2L1", "KLF9", "NKX2-3", "PPP1R15A", "SCARA5",
                  "SEMA3B", "SESN2", "SMPD3", "WDR37", "ZFP36")

# 筛选表达矩阵
expr_markers <- counts_clean_annotated[rownames(counts_clean_annotated) %in% marker_genes, ]

# log2 转换 + z-score 标准化
expr_scaled <- t(scale(t(log2(expr_markers + 1))))
# 确保 colnames 已是你想要的顺序和格式
colnames(expr_scaled) <- c(
  paste0("Normal_", 1:17),
  paste0("Tumor_", 1:17),
  "COLO320_1", "COLO320_2",
  "HCT116_1", "HCT116_2",
  "HT29_1", "HT29_2",
  "RKO_1", "RKO_2",
  "SW480_1", "SW480_2"
)

# 分组信息：Normal / Tumor（这里使用你现有的 sample_info）
sample_info <- data.frame(
  row.names = colnames(expr_scaled),  # 保证行名是样本名
  status = c(rep("Normal", 17), rep("Tumor", 17), rep("CellLine", 10)),
  cell_line = c(rep(NA, 34), rep("COLO320", 2), rep("HCT116", 2), rep("HT29", 2), rep("RKO", 2), rep("SW480", 2))
)

sample_order <- c( "Normal_1", "Normal_2","Normal_3", "Normal_4","Normal_5","Normal_6","Normal_7",  "Normal_8", "Normal_9","Normal_10", "Normal_11","Normal_12",  "Normal_13","Normal_14","Normal_15","Normal_16","Normal_17",       
                   "Tumor_1", "Tumor_2", "Tumor_3","Tumor_4","Tumor_5","Tumor_6","Tumor_7","Tumor_8","Tumor_9","Tumor_10","Tumor_11","Tumor_12","Tumor_13","Tumor_14","Tumor_15","Tumor_16","Tumor_17",
                   "COLO320_1", "COLO320_2",
                  "HCT116_1", "HCT116_2",
                  "HT29_1", "HT29_2",
                  "RKO_1", "RKO_2",
                  "SW480_1", "SW480_2")
# 假设你知道每个 GSM 对应什么类型
colnames(expr_scaled) <- c(
  paste0("Normal_", 1:17),
  paste0("Tumor_", 1:17),
  "COLO320_1", "COLO320_2",
  "HCT116_1", "HCT116_2",
  "HT29_1", "HT29_2",
  "RKO_1", "RKO_2",
  "SW480_1", "SW480_2"
)
expr_scaled <- expr_scaled[, sample_order]
sample_info <- sample_info[sample_order, ]

sample_info <- sample_info[colnames(expr_scaled), ]
ann_colors <- list(
  Type = c(
    "Normal" = "#a1d99b",
    "Tumor" = "#de2d26",
    "COLO320" = "#fbb4ae",
    "HCT116" = "#b3cde3",
    "HT29" = "#ccebc5",
    "RKO" = "#decbe4",
    "SW480" = "#fed9a6"
  )
)
#group <- sample_info$condition
#names(group) <- rownames(sample_info)  # sample name 作为名字

# 转换为颜色标注
group_col <- ifelse(group == "Control", "skyblue", "firebrick")

# 检查 sample_info 的结构是否是 data.frame 且无 NA
str(sample_info)
# 应该输出 status 和 cell_line 两列，都是 factor 或 character
sample_info$cell_line <- NA
sample_info$cell_line[grep("COLO320", colnames(expr_scaled))] <- "COLO320"
sample_info$cell_line[grep("HCT116", colnames(expr_scaled))] <- "HCT116"
sample_info$cell_line[grep("HT29", colnames(expr_scaled))]  <- "HT29"
sample_info$cell_line[grep("RKO", colnames(expr_scaled))]   <- "RKO"
sample_info$cell_line[grep("SW480", colnames(expr_scaled))] <- "SW480"
# 创建注释
ha <- HeatmapAnnotation(
    Status = sample_info$status,
    CellLine = sample_info$cell_line,
    col = list(
    Status = c("Normal" = "#a1d99b", "Tumor" = "#de2d26", "CellLine" = "#9ecae1"),
    CellLine = c("COLO320" = "#fbb4ae", "HCT116" = "#b3cde3",
                 "HT29" = "#ccebc5", "RKO" = "#decbe4", "SW480" = "#fed9a6")
  ),
  annotation_name_side = "left"
)

# 行注释：添加每个基因表达均值的 barplot
row_means <- rowMeans(expr_scaled)

row_anno <- rowAnnotation(
  Expr = anno_barplot(row_means, 
                      gp = gpar(fill = "gray50", col = NA),
                      border = FALSE,
                      width = unit(2, "cm"))
)
pdf("marker_gene_heatmap.pdf", width = 10, height = 6)
Heatmap(
  expr_scaled,
  name = "z-score",
  top_annotation = ha,  # 用之前构建好的 HeatmapAnnotation
  right_annotation = NULL,
  left_annotation = row_anno,  # 如果你有行注释
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,       # 控制是否显示列名
  show_row_names = TRUE,           # 是否显示行名
  row_names_gp = gpar(fontsize = 10),
)
dev.off()
# ===== 富集分析（GO/KEGG）=====
res_df <- results
res_df$Geneid <- rownames(res_df)
degs <- res_df %>% filter(adj.P.Val  < 0.05 & abs(logFC) > 1)
degs <- na.omit(degs)
entrez_ids <- rownames(degs)
head(entrez_ids,)   
library(clusterProfiler)
library(org.Hs.eg.db)

deg_genes <- degs$Geneid  # 假设你差异分析结果里是 SYMBOL

deg_entrez <- bitr(
  deg_genes,
  fromType = "SYMBOL",        # 实际类型
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(deg_entrez)
go <- enrichGO(gene = deg_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", keyType="ENTREZID", readable=TRUE)
write.csv(as.data.frame(go), file = "go_enrichment.csv", row.names = FALSE)
dotplot(go, showCategory=20)
ggsave("go_dotplot.png")

# 读取富集结果
go_df <- read.csv("go_enrichment.csv")
head(go)
go_meta <- go_df[, c("ID", "Description", "p.adjust")]
head(go_meta)
colnames(go_meta) <- c("nodeId", "description", "adjustedPValue")
write.csv(go_meta, "go_node_meta_patch.csv", row.names = FALSE)
# 生成 GO 的边：gene — enriched_in — GO_term
go_edges <- do.call(rbind, lapply(1:nrow(go_df), function(i) {
  genes <- unlist(strsplit(go_df$geneID[i], "/"))
  data.frame(
    source = genes,
    target = go_df$ID[i],              # GO ID，如 "GO:0006281"
    relation = "enriched_in_go"
  )
}))
write.csv(go_edges, "go_edges.csv", row.names = FALSE)
head(go_edges)

kegg <- enrichKEGG(gene=deg_entrez$ENTREZID, organism="hsa")
write.csv(as.data.frame(kegg), "kegg_result.csv", row.names = FALSE)
dotplot(kegg, showCategory=20)
ggsave("kegg_dotplot.png")

# KEGG边：gene — part_of — KEGG_pathway
kegg_edges <- do.call(rbind, lapply(1:nrow(kegg), function(i) {
  genes <- unlist(strsplit(kegg$geneID[i], "/"))
  data.frame(
    source = genes,
    target = kegg$ID[i],           # KEGG pathway ID，如 "hsa04115"
    relation = "part_of_kegg"
  )
}))
write.csv(kegg_edges, "kegg_edges.csv", row.names = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)

# 转换 symbol → Entrez ID
deg_entrez <- bitr(degs$Geneid,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# 合并进原始 degs
degs_annotated <- merge(degs, deg_entrez, by.x = "Geneid", by.y = "SYMBOL")
head(degs_annotated)
#假设你有 deg_df，包括 SYMBOL 和 logFC
go_edges$logFC <- degs$logFC[match(go_edges$source, degs$ENSEMBL)]
go_edges$p.adjust <- go_df$p.adjust [match(go_edges$target, go_df$ID)]
head(go_edges)
write.csv(go_edges, "go_edges_with_attr.csv", row.names = FALSE)

kegg_edges$logFC <- degs_annotated$logFC[match(kegg_edges$source, degs_annotated$ENTREZID)]
kegg_df <- as.data.frame(kegg)
kegg_edges$p.adjust <- kegg_df$p.adjust [match(kegg_edges$target, kegg_df$ID)]
head(kegg_edges)
write.csv(kegg_edges, "kegg_edges_attr.csv", row.names = FALSE)

# 合并两种边
edges_all <- rbind(go_edges, kegg_edges)

# 提取所有节点
all_nodes <- unique(c(edges_all$source, edges_all$target))
# 简单地根据格式给出节点类型（可改进）
node_type <- sapply(all_nodes, function(x) {
  if (grepl("^GO:", x)) return("GO_term")
  else if (grepl("^hsa", x)) return("KEGG_pathway")
  else return("gene")
})

nodes_df <- data.frame(
  id = all_nodes,
  type = node_type
)
write.csv(nodes_df, "nodes.csv", row.names = FALSE)

# 表达矩阵（行为基因，列为样本），例如 limma分析前的log2表达量矩阵
dim(expr_matrix)
head(colnames(expr_matrix))  # 样本名
head(rownames(expr_matrix))  # 基因名
log2_expr_clean <- expr_matrix[complete.cases(expr_matrix), ]
# 过滤掉方差为 0 的基因
var_genes <- apply(log2_expr_clean, 1, var)
log2_expr_filtered <- log2_expr_clean[var_genes > 0, ]

# 再次计算 PCA
pca <- prcomp(t(log2_expr_filtered), scale. = TRUE)
# 取前两个主成分
group <- c(
  rep("Normal", 17),
  rep("Tumor", 17),
  rep("CellLine", 10)
)

cell_line <- c(
  rep(NA, 34),  # Normal 和 Tumor 没有 cell line 名
  rep("COLO320", 2),
  rep("HCT116", 2),
  rep("HT29", 2),
  rep("RKO", 2),
  rep("SW480", 2)
)
pca_df <- data.frame(
  Sample = colnames(log2_expr_filtered),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = factor(group),
  CellLine = factor(cell_line)
)

# ggplot2 绘图
library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")
ggsave("PCA of Samples.png")

# 安装并加载 KEGGREST 包
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}
library(KEGGREST)
kegg_pathways <- keggList("pathway", "hsa")
head(kegg_pathways)
kegg_df <- data.frame(
  id = names(kegg_pathways),
  description = gsub(" - Homo sapiens \\(human\\)", "", unname(kegg_pathways)),
  stringsAsFactors = FALSE
)
head(kegg_df)
write.csv(kegg_df, "kegg_patch.csv", row.names = FALSE)
