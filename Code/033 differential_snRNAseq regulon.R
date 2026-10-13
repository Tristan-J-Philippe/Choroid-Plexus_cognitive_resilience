###
### 033 differential_RNA regulon.R
###
# Purpose: Making regulon activity matrix and differential analysis thereof
# DependeCtrles:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(42)
#Scenic
##Prepare data for PyScenic

library(Matrix)
library(data.table)
library(Seurat)
library(SeuratDisk)


DefaultAssay(Ex_ChP) <- "RNA"
#remove genes that are not expressed
my_genes <- c(unique(limma_resilientvAD_Epithelial$gene, limma_resilientvAD_Fibroblast$gene, limma_resilientvAD_Mural$gene, limma_resilientvAD_Endothelial$gene, limma_resilientvAD_Immune$gene, limma_resilientvAD_Neuro$gene, limma_resilientvAD_Oligo$gene))
Ex_ChP <- Ex_ChP[my_genes, ]

SaveH5Seurat(Ex_ChP, filename = paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP_full.h5Seurat"), overwrite = TRUE)
Convert(paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP_full.h5Seurat"), assay = "RNA", dest = "h5ad", overwrite = T)

#Random downsampling for gnrboost2 function otherwise it may take weeks
set.seed(123)
Ex_ChP_sub <- subset(x = Ex_ChP, cells = sample(x = colnames(Ex_ChP), size = 100000, replace = FALSE))

SaveH5Seurat(Ex_ChP_sub, filename = paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP_sub.h5Seurat"), overwrite = TRUE)
Convert(paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP_sub.h5Seurat"), assay = "RNA", dest = "h5ad", overwrite = T)


#!!!!Go to 033 differential_snRNAseq regulon.py

##Load back in R and prep for limma

library(data.table)
library(Seurat)
library(Matrix)

dat <- read_csv(paste0(curr_dir, "/Scenic/results/auc_mtx_.csv"))
dat <- column_to_rownames(dat, "Cell")
meta <- data.frame(Ex_ChP@meta.data) #get metadata from seuratobject

unique(rownames(dat)==rownames(meta))
dat$cell <- rownames(dat)
meta$cell <- rownames(meta)

regulons <- colnames(dat)
regulons <- regulons[-length(regulons)] 

#Aggregate data per Lane_projid and celltype
auc_meta_Major <- left_join(dat, meta[,c("cell", "Lane_projid", "Majorcelltype")], by = "cell")
auc_meta_Major <- auc_meta_Major %>% group_by(Lane_projid, Majorcelltype) %>% filter(n() >= 50) %>% ungroup() #remove any combination with < 50 cells.
auc_Major <- auc_meta_Major %>% group_by(Lane_projid, Majorcelltype) %>%
  summarise(across(.cols = all_of(regulons), \(x) mean(x, na.rm = TRUE)), n_cells = n(), .groups = "drop") #get mean AUC across all cells for each sample

for(celltype in auc_Major$Majorcelltype){
  cts <- auc_Major[auc_Major$Majorcelltype==celltype,]
  cts <- data.frame(cts)
  cts$Majorcelltype <- NULL
  cts$n_cells <- NULL
  rownames(cts) <- cts$Lane_projid
  cts$Lane_projid <- NULL
  cts <- t(cts)
  rownames(cts) <- gsub("...", "", rownames(cts), fixed = TRUE)
  assign(paste("cts_TF_", "Mural", sep=""),  cts, envir = .GlobalEnv)
}

auc_meta_Sub <- left_join(dat, meta[,c("cell", "Lane_projid", "Subcelltype")], by = "cell")
auc_meta_Sub <- auc_meta_Sub %>% group_by(Lane_projid, Subcelltype) %>% filter(n() >= 50) %>% ungroup() #remove any combination with < 50 cells.
auc_Sub <- auc_meta_Sub %>% group_by(Lane_projid, Subcelltype) %>%
  summarise(across(.cols = all_of(regulons), \(x) mean(x, na.rm = TRUE)), n_cells = n(), .groups = "drop") #get mean AUC across all cells for each sample

for(celltype in auc_Sub$Subcelltype){
  cts <- auc_Sub[auc_Sub$Subcelltype==celltype,]
  cts <- data.frame(cts)
  cts$Subcelltype <- NULL
  cts$n_cells <- NULL
  rownames(cts) <- cts$Lane_projid
  cts$Lane_projid <- NULL
  cts <- t(cts)
  rownames(cts) <- gsub("...", "", rownames(cts), fixed = TRUE)
  assign(paste("cts_TF_", celltype, sep=""),  cts, envir = .GlobalEnv)
}


set.seed(42)

#Set and load dependent variables for limma/DESeq functions
Gene.search <- read.csv("/home/tristan-j-philippe/DEG to GO.csv")
pcutoff <- 0.05
LFcutoff <- 0.001

###cogdx
####Factors set up

#metadata
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly 3334.csv"))
mdat <- column_to_rownames(mdat, "projid")

msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))
Cog_Path = factor(paste(mdat$Cog_Path))
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)

####Resilience

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path + msex+ age + pmi + X10X_batch)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdxadnc_resilientvCtrl=Cog_Pathresilient-Cog_PathCtrl, cogdxadnc_ADdvresilient=Cog_PathADd-Cog_Pathresilient,  cogdxadnc_ADdvCtrl=Cog_PathADd-Cog_PathCtrl)

TFlimmafun_lane(cts_TF_Epithelial)
TFlimmafun_lane(cts_TF_Fibroblast)
TFlimmafun_lane(cts_TF_Endothelial)
TFlimmafun_lane(cts_TF_Immune)

TFlimmafun_lane(cts_TF_Epi_1a)
TFlimmafun_lane(cts_TF_Epi_1b)
TFlimmafun_lane(cts_TF_Epi_2a)
TFlimmafun_lane(cts_TF_Epi_2b)
TFlimmafun_lane(cts_TF_Fib_1)
TFlimmafun_lane(cts_TF_Fib_2)
TFlimmafun_lane(cts_TF_Fib_3)


##Figures

sig_TFs <- c(unique(c(limma_cogdxadnc_resilientvCtrl_TF_Epithelial_sig$gene, limma_cogdxadnc_resilientvCtrl_TF_Fibroblast_sig$gene, limma_cogdxadnc_resilientvCtrl_TF_Immune_sig$gene, limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig$gene, limma_cogdxadnc_ADdvCtrl_TF_Fibroblast_sig$gene, limma_cogdxadnc_ADdvCtrl_TF_Immune_sig$gene)))

lm_merge <- full_join(limma_cogdxadnc_resilientvCtrl_TF_Epithelial, limma_cogdxadnc_resilientvCtrl_TF_Fibroblast, by="gene", suffix =c("", "_Fibroblast_cogdxadnc_resilientvCtrl")) %>% full_join(limma_cogdxadnc_resilientvCtrl_TF_Endothelial, by="gene", suffix=c("", "_Endothelial_cogdxadnc_resilientvCtrl")) %>% full_join(limma_cogdxadnc_resilientvCtrl_TF_Immune, by="gene", suffix=c("", "_Immune_cogdxadnc_resilientvCtrl")) %>% full_join(limma_cogdxadnc_ADdvCtrl_TF_Epithelial, by="gene", suffix=c("", "_Epithelial_ADdvCtrl")) %>% full_join(limma_cogdxadnc_ADdvCtrl_TF_Fibroblast, by="gene", suffix=c("", "_Fibroblast_ADdvCtrl")) %>% full_join(limma_cogdxadnc_ADdvCtrl_TF_Endothelial, by="gene", suffix=c("", "_Endothelial_ADdvCtrl")) %>% full_join(limma_cogdxadnc_ADdvCtrl_TF_Immune, by="gene", suffix=c("", "_Immune_ADdvCtrl"))
lm_merge <- column_to_rownames(lm_merge, "gene")
lm_merge <- lm_merge[sig_TFs,]
mat <- lm_merge[,c(grep(pattern = "regulon_activity", x = colnames(lm_merge), value = T))]
colnames(mat) <- c("Epithelial_cogdxadnc_resilientvCtrl", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("regulon_activity_", "", colnames(mat))

gene_dist <- dist(mat)
gene_hc <- hclust(gene_dist, method = "complete")
gene_order <- rownames(mat)[gene_hc$order]

TF_tbb <- mat %>% rownames_to_column("gene") %>% pivot_longer(cols=colnames(mat), names_to="cell_group", values_to="change_AUC")
TF_tbb <- TF_tbb %>% separate(cell_group, into = c("Majorcelltype", "Comparison"), sep = "_", extra = "merge")

TF_tbb <- mutate(TF_tbb, gene = factor(gene, levels = gene_order),
                 Majorcelltype = factor(Majorcelltype, levels = c("Epithelial", "Fibroblast", "Endothelial", "Immune")),
                 Cog_Path = factor(Comparison, levels=c("cogdxadnc_resilientvCtrl", "ADdvCtrl")))

ggplot(TF_tbb, aes(x = Comparison, y = gene)) +
  geom_point(aes(color = change_AUC)) +
  scale_color_viridis_c(option = "plasma") + 
  facet_grid(. ~ Majorcelltype, scales = "free_x", space = "free_x") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.1, "lines"), strip.background = element_rect(fill = "white"), strip.text = element_text(size = 7)) +
  labs(x = "Major Cell Type by Comparison", y = "Transcription factor (Clustered across all groups)", color = "change_AUC")



auc_Major$projid <- substr(auc_Major$Lane_projid, 1, 8)
mdat$projid <- rownames(mdat)
TF_tbb <- left_join(auc_Major, data.frame("projid"=mdat$projid, "Cog_Path"=mdat$Cog_Path), by="projid")
TF_tbb$Lane_projid <- NULL
TF_tbb$projid <- NULL
TF_tbb$n_cells <- NULL

TF_tbb <- TF_tbb %>% group_by(Cog_Path, Majorcelltype) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
TF_tbb <- subset(TF_tbb, Majorcelltype=="Epithelial" | Majorcelltype=="Fibroblast" | Majorcelltype=="Endothelial" | Majorcelltype=="Immune")
TF_tbb <- subset(TF_tbb, Cog_Path=="Ctrl" | Cog_Path=="resilient" | Cog_Path=="ADd")
colnames(TF_tbb) <- gsub("(+)", "", colnames(TF_tbb), fixed = TRUE)

sig_genes <- c(unique(c(limma_cogdxadnc_resilientvCtrl_TF_Epithelial_sig$gene, limma_cogdxadnc_resilientvCtrl_TF_Fibroblast_sig$gene,  limma_cogdxadnc_resilientvCtrl_TF_Immune_sig$gene, limma_cogdxadnc_ADdvresilient_TF_Epithelial_sig$gene, limma_cogdxadnc_ADdvresilient_TF_Fibroblast_sig$gene,  limma_cogdxadnc_ADdvresilient_TF_Immune_sig$gene, limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig$gene, limma_cogdxadnc_ADdvCtrl_TF_Fibroblast_sig$gene, limma_cogdxadnc_ADdvCtrl_TF_Endothelial_sig$gene, limma_cogdxadnc_ADdvCtrl_TF_Immune_sig$gene)))
TF_tbb <- TF_tbb[,colnames(TF_tbb) %in% c(sig_genes, "Majorcelltype", "Cog_Path")]

TF_mtx <- mutate(TF_tbb, cell_group = paste(Majorcelltype, Cog_Path, sep = "_"))
TF_mtx$Cog_Path <- NULL
TF_mtx$Majorcelltype <- NULL
TF_mtx <- data.frame(TF_mtx)
rownames(TF_mtx) <- TF_mtx$cell_group
TF_mtx$cell_group <- NULL
TF_mtx <- t(TF_mtx)

gene_dist <- dist(TF_mtx)
gene_hc <- hclust(gene_dist, method = "complete")
gene_order <- rownames(TF_mtx)[gene_hc$order]

TF_tbb <- pivot_longer(TF_tbb, cols=-c("Majorcelltype", "Cog_Path"), names_to="gene", values_to="AUC")
TF_tbb <- mutate(TF_tbb, gene = factor(gene, levels = gene_order),
                 Majorcelltype = factor(Majorcelltype, levels = c("Epithelial", "Fibroblast", "Endothelial", "Immune")),
                 Cog_Path = factor(Cog_Path, levels=c("Ctrl", "resilient", "ADd")))

ggplot(TF_tbb, aes(x = Cog_Path, y = gene)) +
  geom_point(aes(color = AUC)) +
  scale_color_viridis_c(option = "plasma") + 
  facet_grid(. ~ Majorcelltype, scales = "free_x", space = "free_x") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.1, "lines"), strip.background = element_rect(fill = "white"), strip.text = element_text(size = 7)) +
  labs(x = "Major Cell Type by Pathology-Cognition", y = "Transcription factor (Clustered across all groups)", color = "AUC")



###Heatmap

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
mdat$projid <- rownames(mdat)


heatfun("Epithelial")
heatfun("Fibroblast")
heatfun("Immune")


library(ggraph)
library(igraph)
library(ggnewscale)
library(tidygraph)

edges <- read_csv("/Ex_ChP/Ex_ChP7/Scenic/results/SCENIC_regulon_edges.csv")
edges$NES <- NULL
edges$Regulon <- NULL

####Epithelial
#####AD

#network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Epithelial_sig, tflist=c(), name="Epithelial_AD_Cilia")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Epithelial_sig, tflist=c("ETV5", "FOSL2", "ETS2"), name="Epithelial_AD_Adhesion")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Epithelial_sig, tflist=c("ARID5B", "VDR"), name="Epithelial_AD_Ion_Transport")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Epithelial_sig, tflist=c("XBP1", "GATA3", "VDR", "RARB", "ETV5"), name="Epithelial_AD_Secretion")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Epithelial_sig, tflist=c("MLXIPL", "FOXO3", "VDR", "RARB", "ARID5B"), name="Epithelial_AD_lipid_Metabolism")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Epithelial_sig, tflist=c("FOXO3", "XBP1", "MAFF", "RELB", "IRF2", "ETS2"), name="Epithelial_AD_Stress_response_Immune")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Epithelial_sig, tflist=c("XBP1", "FOXO3", "MLXIPL", "ARID5B"), name="Epithelial_AD_protein_Metabolism")


#####Resilient

network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Epithelial_sig, tflist=c("XBP1", "GATA3", "VDR", "RARB", "ETV5"), name="Epithelial_resilient_Secretion")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Epithelial_sig, tflist=c("MLXIPL", "FOXO3", "VDR", "RARB", "ARID5B"), name="Epithelial_resilient_lipid_Metabolism")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Epithelial_sig, tflist=c("FOXO3", "XBP1", "MAFF", "RELB", "IRF2", "ETS2"), name="Epithelial_resilient_Stress_response")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Epithelial_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Epithelial_sig, tflist=c("XBP1", "FOXO3", "MLXIPL", "ARID5B"), name="Epithelial_resilient_protein_Metabolism")


####Fibroblast
#####AD

network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Fibroblast_sig, tflist=c("FOXC2", "TEAD2", "NFATC1"), name="Fibroblast_AD_Adhesion")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Fibroblast_sig, tflist=c("NFACT1", "VDR"), name="Fibroblast_AD_Ion_Transport")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Fibroblast_sig, tflist=c("HIF1A", "XBP1", "FOSL2"), name="Fibroblast_AD_Stress_response")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Fibroblast_sig, tflist=c("NFKB1", "RELB", "ETS2", "PRDM1"), name="Fibroblast_AD_Immune")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Fibroblast_sig, tflist=c("HIF1A", "VDR", "XBP1", "KLF5"), name="Fibroblast_AD_protein_Metabolism")


#####resilient

network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Fibroblast_sig, tflist=c("FOXC2", "TEAD2", "NFATC1"), name="Fibroblast_resilient_Adhesion")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Fibroblast_sig, tflist=c("NFACT1", "VDR"), name="Fibroblast_resilient_Ion_Transport")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Fibroblast_sig, tflist=c("HIF1A", "XBP1", "FOSL2"), name="Fibroblast_resilient_Stress_response")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Fibroblast_sig, tflist=c("NFKB1", "RELB", "ETS2", "PRDM1"), name="Fibroblast_resilient_Immune")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Fibroblast_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Fibroblast_sig, tflist=c("HIF1A", "VDR", "XBP1", "KLF5"), name="Fibroblast_resilient_protein_Metabolism")


####Immune
#####AD

network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Immune_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Immune_sig, tflist=c("HIF1A", "MAFF"), name="Immune_AD_lipid_Metabolism")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Immune_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Immune_sig, tflist=c("HIF1A", "NFKB1", "FOSL2", "STAT4", "FOSL2", "ETS2", "CHD1"), name="Immune_AD_Stress_response_Immune")
network_plot(lmTF=limma_cogdxadnc_ADdvCtrl_TF_Immune_sig, lmdeg=limma_cogdxadnc_ADdvCtrl_Immune_sig, tflist=c("HIF1A", "CHD1"), name="Immune_AD_protein_Metabolism")


#####Resilient

network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Immune_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Immune_sig, tflist=c("HIF1A", "MAFF"), name="Immune_resilient_lipid_Metabolism")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Immune_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Immune_sig, tflist=c("HIF1A", "NFKB1", "FOSL2", "STAT4", "FOSL2", "ETS2", "MAFF"), name="Immune_resilient_Stress_response_Immune")
network_plot(lmTF=limma_cogdxadnc_resilientvCtrl_TF_Immune_sig, lmdeg=limma_cogdxadnc_resilientvCtrl_Immune_sig, tflist=c("HIF1A", "CHD1"), name="Immune_resilient_protein_Metabolism")


