setwd("D:/Rproject/scRNA/ZenghanAMU")

library(Seurat)
library(ggplot2)
library(GSVA)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(pheatmap)

Epithelial <- readRDS('E:/BioBank/scRNA/ZenghanAMU/Epithelial.rds')
Epithelial_MHC2high <- readRDS('E:/BioBank/scRNA/ZenghanAMU/Epithelial_MHC2high.rds')
Epithelial_MHC2low <- readRDS('E:/BioBank/scRNA/ZenghanAMU/Epithelial_MHC2low.rds')

gmt <- "E:/BioBank/Index/GSVA/v7.2.symbols_PCa_gene.sets.gmt"
gene_sets <- read.gmt(gmt)
# 定义感兴趣的通路名称
pathways_of_interest <- c(
  "Notch_signaling", "Hedgehog_signaling", "PI3K_AKT_mTOR_signaling", "mRTORC1_signaling", 
  "E2F_targets", "Myc_targets_V1", "Myc_targets_V2", "NMYC_01", "P53_pathway", "Loss of TP53 DN", 
  "Kras_signaling_UP", "Kras_signaling_DN", "RAS", "Loss of PTEN", "FOXA1", "Loss of SPOP", "EZH2_targets", 
  "Hippo", "SOX2_targets", "TMPRSS2 ERG fusion UP", "PRC", "Prostate_cancer_UP", "Loss of RB1 UP", "Proliferation", 
  "Cell_cycle", "G2M_checkpoint", "DNA_repair", "Apoptosis", "Stemness", "Stem.cell", "Senescence", 
  "Differentiation", "Mitotic_spindle", "Androgen_response", "AR_activity", "AR_V", "Lineage Plasticity Risk Signature", 
  "LPSig", "LPC_gene_signature", "CRPCsig51", "NE_differentiation", "NE_Sig_UP", "NE_Sig_DN", "PN", 
  "GB_plasticity_UP", "CCP", "Antigen_processing_and_presentation", "Interferon_alpha_response", 
  "Interferon_gamma_response", "TNFA_signaling_via_NFKB", "Inflammatory_response", "IL6_jak_stat3_signaling", 
  "IL2_stat5_signaling", "JAK_STAT", "JAK_STAT_FGFR", "WNT_beta_catenin_signaling", "TGF_beta_signaling", 
  "EMT", "EMT_PCa", "Angiogenesis", "Coagulation", "Myogenesis", "Hypoxia", "Glycolysis", 
  "Oxidative_phosphorylation", "Cholesterol_homeostasis", "Adipogenesis", "Xenobiotic_metabolism", 
  "Fatty_acid_metabolism", "Bile_acid_metabolism", "Reactive_oxygen_species_Pathway", "Heme_metabolism", 
  "Complement", "Peroxisome", "Allograft_rejection", "Stress"
)
gene_sets_filtered <- gene_sets[gene_sets$term %in% pathways_of_interest, ]
gene_sets_filtered_clean <- gene_sets_filtered %>%
  filter(!is.na(gene) & gene != "")
gene_sets <- gene_sets_filtered_clean
gene_sets_list <- split(gene_sets$gene, gene_sets$term)

head(Epithelial@meta.data[["predicted_celltype_BCL"]])
unique(Epithelial@meta.data[["predicted_celltype_BCL"]])

# get Matrix epithelial
Epithelial_MHC2high_expr <- GetAssayData(Epithelial_MHC2high, assay = "RNA", slot = "data")
Epithelial_MHC2low_expr <- GetAssayData(Epithelial_MHC2low, assay = "RNA", slot = "data")
Epithelial_MHC2high_matrix <- as.matrix(Epithelial_MHC2high_expr)
Epithelial_MHC2low_matrix <- as.matrix(Epithelial_MHC2low_expr)
meanExprHigh <- rowMeans(Epithelial_MHC2high_matrix)
meanExprLow <- rowMeans(Epithelial_MHC2low_matrix)
meanEpithelialExp <- data.frame(
  EpithelialMHC2High = meanExprHigh,
  EpithelialMHC2Low = meanExprLow
)
meanEpithelialExp <- as.matrix(meanEpithelialExp)
epithelialGSVA <- gsva(meanEpithelialExp, 
                       gset.idx.list = gene_sets_list,  
                     method = "gsva", 
                     kcdf = "Poisson", 
                     parallel.sz = 30)
MHC2 <- pheatmap(epithelialGSVA, show_colnames = T, 
         scale = "none",angle_col = "45",
         cluster_row = F,cluster_col = F,
         color = colorRampPalette(c("#7DA6C6", "white", "#E68B81"))(50))
ggsave("MHC2.png", plot = MHC2, width = 8, height = 20)



# get matrix 3 types
# get cell label
LuminalHigh <- WhichCells(Epithelial_MHC2high, expression = predicted_celltype_BCL == "Luminal")
ClubHigh <- WhichCells(Epithelial_MHC2high, expression = predicted_celltype_BCL == "Club")
BasalHigh <- WhichCells(Epithelial_MHC2high, expression = predicted_celltype_BCL == "Basal")
LuminalLow <- WhichCells(Epithelial_MHC2low, expression = predicted_celltype_BCL == "Luminal")
ClubLow <- WhichCells(Epithelial_MHC2low, expression = predicted_celltype_BCL == "Club")
BasalLow <- WhichCells(Epithelial_MHC2low, expression = predicted_celltype_BCL == "Basal")
# get cell matrix
exprLumHigh <- as.matrix(Epithelial_MHC2high_expr[, LuminalHigh])
exprClubHigh <- as.matrix(Epithelial_MHC2high_expr[, ClubHigh])
exprBasalHigh <- as.matrix(Epithelial_MHC2high_expr[, BasalHigh])
exprLumLow <- as.matrix(Epithelial_MHC2low_expr[, LuminalLow])
exprClubLow <- as.matrix(Epithelial_MHC2low_expr[, ClubLow])
exprBasalLow <- as.matrix(Epithelial_MHC2low_expr[, BasalLow])
# calculate mean
meanExprHighLuminal <- rowMeans(exprLumHigh)
meanExprHighClub <- rowMeans(exprClubHigh)
meanExprHighBasal <- rowMeans(exprBasalHigh)
meanExprLowLuminal <- rowMeans(exprLumLow)
meanExprLowClub <- rowMeans(exprClubLow)
meanExprLowBasal <- rowMeans(exprBasalLow)
# combine
meanExpr3type <- data.frame(
  LuminalHigh = meanExprHighLuminal,
  ClubHigh = meanExprHighClub,
  BasalHigh = meanExprHighBasal,
  LuminalLow = meanExprLowLuminal,
  ClubLow = meanExprLowClub,
  BasalLow = meanExprLowBasal
)
meanExpr3type <- as.matrix(meanExpr3type)
# GSVA
BCL_GSVA <- gsva(meanExpr3type, 
                       gset.idx.list = gene_sets_list,  
                       method = "gsva", 
                       kcdf = "Poisson", 
                       parallel.sz = 30)
BCL <- pheatmap(BCL_GSVA, show_colnames = T, 
                 scale = "row",angle_col = "45",
                 cluster_row = F,cluster_col = F,
                 color = colorRampPalette(c("#7DA6C6", "white", "#E68B81"))(50))
ggsave("BCL.png", plot = BCL, width = 10, height = 20)

