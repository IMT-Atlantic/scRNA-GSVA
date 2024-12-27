setwd("D:/Rproject/scRNA/Irg1KO")

library(Seurat)
library(RColorBrewer)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(pheatmap)
library(clusterProfiler)
library(tidyr)
library(limma)
library(edgeR)
library(grid)
library(ggplot2)

library(parallel)
library(BiocParallel)
library(snow)

register(SnowParam(workers = 30)) 

Data <- readRDS('E:/BioBank/scRNA/GSE215195ITAscRNA/DataIndex.rds')
Data <- JoinLayers(Data)
DimPlot(Data, reduction = "umap", group.by = 'cell_type_annotation', label=TRUE)

AM <- subset(Data, subset = cell_type_annotation == "Alveolar Macrophage")
DimPlot(AM, reduction = "umap", group.by = 'orig.ident', label=TRUE)

DiffGene <- readRDS("EnrichList.rds")
mouse <- msigdbr(species = "Mus musculus")
table(mouse$gs_cat)
table(mouse$gs_subcat)

# 针对AM中PBS_WT和PBS_IRG1KO组进行C8富集分析
c8_gene_sets_df <- msigdbr(species = "Mus musculus", category = "C8")
c8_gene_sets <- c8_gene_sets_df %>%
  split(x = .$gene_symbol, f = .$gs_name)
# 提取 PBS_WT 和 PBS_IRG1KO 组的细胞
unique(AM@meta.data$orig.ident)
AM <- SetIdent(AM, value = "orig.ident") # 指定metadata中的一部分作为细胞标签
PBS_WT_cells <- WhichCells(AM, idents = "PBS_WT")
PBS_IRG1KO_cells <- WhichCells(AM, idents = "PBS_IRG1KO")
LAC_WT <- WhichCells(AM, idents = "LAC_WT")
LAC_IRG1KO <- WhichCells(AM, idents = "LAC_IRG1KO")

# 提取基因表达数据（假设使用 RNA assay）
expr_data <- GetAssayData(AM, assay = "RNA", slot = "data")

# 提取两个组的表达数据
expr_data_PBS_WT <- expr_data[, PBS_WT_cells]
expr_data_PBS_IRG1KO <- expr_data[, PBS_IRG1KO_cells]
expr_data_PBS_WT_matrix <- as.matrix(expr_data_PBS_WT)
expr_data_PBS_IRG1KO_matrix <- as.matrix(expr_data_PBS_IRG1KO)
expr_data_LAC_WT <- expr_data[, LAC_WT]
expr_data_LAC_WT_matrix <- as.matrix(expr_data_LAC_WT)
expr_data_LAC_IRG1KO <- expr_data[, LAC_IRG1KO]
expr_data_LAC_IRG1KO_matrix <- as.matrix(expr_data_LAC_IRG1KO)





#———————————————————————————————————————————————————————————————————————————————use all cell calculate GSVA Score
# 输入矩阵为log转化后的连续表达矩阵指则设置kcdf参数为"Gaussian"，如果是counts矩阵则设置kcdf为"Poisson"
# 这里无法使用常规方法进行并行计算
GSVA_PBS_WT <- gsva(
  expr = expr_data_PBS_WT_matrix,
  gset.idx.list = c8_gene_sets,
  kcdf = "Poisson",  
  parallel.sz = 30    
)

GSVA_PBS_IRG1KO <- gsva(
  expr = expr_data_PBS_IRG1KO_matrix,
  gset.idx.list = c8_gene_sets,
  kcdf = "Poisson",  
  parallel.sz = 30    
)

GSVA_LAC_WT <- gsva(
  expr = expr_data_LAC_WT,
  gset.idx.list = c8_gene_sets,
  kcdf = "Poisson",  
  parallel.sz = 30   
)

GSVA_LAC_IRG1KO <- gsva(
  expr = expr_data_LAC_IRG1KO,
  gset.idx.list = c8_gene_sets,
  kcdf = "Poisson", 
  parallel.sz = 30   
)

# 确保转换为矩阵
GSVA_LAC_IRG1KO <- as.matrix(GSVA_LAC_IRG1KO)
GSVA_LAC_WT <- as.matrix(GSVA_LAC_WT)
GSVA_PBS_IRG1KO <- as.matrix(GSVA_PBS_IRG1KO)
GSVA_PBS_WT <- as.matrix(GSVA_PBS_WT)
# 获取所有数据集的通路名称的交集
common_pathways <- Reduce(intersect, list(
  rownames(GSVA_LAC_IRG1KO),
  rownames(GSVA_LAC_WT),
  rownames(GSVA_PBS_IRG1KO),
  rownames(GSVA_PBS_WT)
))
# 子集化每个数据集，只保留共同的通路
GSVA_LAC_IRG1KO <- GSVA_LAC_IRG1KO[common_pathways, , drop = FALSE]
GSVA_LAC_WT <- GSVA_LAC_WT[common_pathways, , drop = FALSE]
GSVA_PBS_IRG1KO <- GSVA_PBS_IRG1KO[common_pathways, , drop = FALSE]
GSVA_PBS_WT <- GSVA_PBS_WT[common_pathways, , drop = FALSE]
GSVA_LAC_IRG1KO <- as.data.frame(GSVA_LAC_IRG1KO)
GSVA_LAC_WT <- as.data.frame(GSVA_LAC_WT)
GSVA_PBS_IRG1KO <-as.data.frame(GSVA_PBS_IRG1KO)
GSVA_PBS_WT <- as.data.frame(GSVA_PBS_WT)
# 保存数据
saveRDS(GSVA_LAC_IRG1KO, 'GSVA_LAC_IRG1KO.rds')
saveRDS(GSVA_LAC_WT, 'GSVA_LAC_WT.rds')
saveRDS(GSVA_PBS_IRG1KO, 'GSVA_PBS_IRG1KO.rds')
saveRDS(GSVA_PBS_WT, 'GSVA_PBS_WT.rds')

# Limma in doubt
# 合并结果 cbind横向合并 rbind纵向合并
PBS_IRG1KO_WT <- cbind(GSVA_PBS_WT, GSVA_PBS_IRG1KO)
group <- c(rep("WT", 622), rep("IRG1KO", 425)) %>% as.factor()#设置分组，对照在前
desigN <- model.matrix(~ 0 + group) #构建比较矩阵
colnames(desigN) <- levels(group)
fit = lmFit(PBS_IRG1KO_WT, desigN)
fit2 <- eBayes(fit)
diff=topTable(fit2,adjust='fdr',coef=2,number=Inf)#校准按照需求修改 ？topTable
write.csv(diff, file = "Diff_IRG1KO_WT.csv")

#———————————————————————————————————————————————————————————————————————————————use CellAverage calculate GSVA Score
# Get average expression and use GSVA
mean_expr_PBS_WT <- rowMeans(expr_data_PBS_WT_matrix)
mean_expr_PBS_IRG1KO <- rowMeans(expr_data_PBS_IRG1KO_matrix)
mean_expr_LAC_WT <- rowMeans(expr_data_LAC_WT_matrix)
mean_expr_LAC_IRG1KO <- rowMeans(expr_data_LAC_IRG1KO_matrix)
mean_expr_df <- data.frame(
  PBS_WT = mean_expr_PBS_WT,
  PBS_IRG1KO = mean_expr_PBS_IRG1KO,
  LAC_WT = mean_expr_LAC_WT,
  LAC_IRG1KO = mean_expr_LAC_IRG1KO
)
mean_expr <- as.matrix(mean_expr_df)

# C8 Proliferation
GSVA <- gsva(
  expr = mean_expr,
  gset.idx.list = c8_gene_sets,
  kcdf = "Poisson", 
  parallel.sz = 30   
)
pheatmap(GSVA[1:30, ], show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# C5 GO
C5 <- msigdbr(species = "Mus musculus", category = "C5")
C5 <- C5 %>% split(x = .$gene_symbol, f = .$gs_name)
GSVA <- gsva(
  expr = mean_expr,
  gset.idx.list = C5,
  kcdf = "Poisson", 
  parallel.sz = 30   
)
GSVA_GO_Top20 <- pheatmap(GSVA[1:20, ], show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
ggsave("GSVA_GO_Top20.png", plot = GSVA_GO_Top20, width = 8, height = 6)
