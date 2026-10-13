###
### 032 differential_ST.R
###
# Purpose: Perform differential analysis of ST dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder
set.seed(42)

pvalcutoff <- 0.05
LFcutoff <- 0.5
pct_cutoff <- 0.05

mdat <- read_csv("/Ex_ChP/CosMX2/CosMX-projID metadata.csv") %>% column_to_rownames("projID")


##Subset
###Final cell type
Cos_ChP <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Cos_ChP.RDS", sep=""))
for(i in names(Cos_ChP@images)){Cos_ChP@images[[i]]<- NULL} #remove images before subbing

###Major cell type
#!!!!Change to Majorcelltype in datprepfun
datprepfun(Cos_ChP, "Epithelial")
datprepfun(Cos_ChP, "Endothelial")
datprepfun(Cos_ChP, "Mural")
datprepfun(Cos_ChP, "Fibroblast")
datprepfun(Cos_ChP, "Immune")

##limma runs
###factors set up
#discrete var
msex = factor(paste(mdat$msex))
Cog_Path <- factor(mdat$Cog_Path, levels = c("Resilient", "AD"))
#continuous
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)

###cogdx
options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path + msex + age + pmi)
rownames(Maindm) <- paste0("g", rownames(mdat))
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, ADdvresilient=Cog_PathAD-Cog_PathResilient)

limmafun(cts_Epithelial, pct_Epithelial)
limmafun(cts_Fibroblast, pct_Fibroblast)
limmafun(cts_Endo, pct_Endo)
limmafun(cts_Immune, pct_Immune)
limmafun(cts_Mural, pct_Mural)



#FGSEA
Path_fun(limma_ADdvresilient_Endo)
Path_fun(limma_ADdvresilient_Mural)
Path_fun(limma_ADdvresilient_Epithelial)
Path_fun(limma_ADdvresilient_Fibroblast)
Path_fun(limma_ADdvresilient_Immune)

#Figures DEG-GSEA
Pvalcutoff <- 0.05
LFcutoff <- 0.2
Epi_genes <- c("TUBB4B", "CFAP43", "DNAH12", "FOXJ1", "SERPINA3", "IL6R", "CHI3L1", "UBB", "PLTP", "NQO1", "KCNN2", "CPE", "LAMB1", "CLDN5", "COL9A3")
Fib_genes <- c("COL8A1", "COL15A1", "ADAMTS9", "CD96", "APOE", "ACSL4", "APOLD1", "KCNMA1", "SLC38A2", "TFRC", "CD9", "GPX4", "C7", "GPX3")


###Load sn

##custom Volcano
Volcano_fun(limma_cogdx_ADvNCI_Epithelial, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_1, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_2, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_3, Epi_genes)
Volcano_fun(limma_cogdx_ADvNCI_Epi_4, Epi_genes)

Volcano_fun(limma_cogdx_ADvNCI_Fibroblast, Fib_genes)


##Boxplots
####Resilience
limma_ADdvresilient_Epithelial_sig <- limma_ADdvresilient_Epithelial[limma_ADdvresilient_Epithelial$P.Value < Pvalcutoff,]
Epi_genes <- c("NQO1", "ODC1", "CTSD", "LAMP2", "CLU", "PINK1", "SQSTM1", "P4HB", "PFKL")
Epi_genes <- Epi_genes[Epi_genes %in% limma_ADdvresilient_Epithelial_sig$gene]
plt <- NULL
for(gene in Epi_genes){
  plts <- box_fun(limma_ADdvresilient_Epithelial_sig, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_Res.pdf", sep=""), plt, width=.7, height=.5*length(Epi_genes))

gene_list <- c("EFNA5", "NRG3", "COL9A1", "NAMPT", "JAM3", "INSR", "ITGAV", "ITGB8", "EPHA4")
gene_list <- gene_list[gene_list %in% limma_ADdvresilient_Epithelial_sig$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_ADdvresilient_Epithelial_sig, "Cog_Path", c("skyblue", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_cellchat_Res.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))


##GO heatmap of overlap
library(ComplexHeatmap)
library(colorRamp2)
library(circlize)
library(DESeq2)
library(viridis)
library(export)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue","white", "red"))


##Resilience
###Epi

selection <- c(#Unique to resilience
  "GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOBP_PROTEIN_LOCALIZATION_TO_ORGANELLE", "GOCC_CATALYTIC_COMPLEX", "GOBP_NEGATIVE_REGULATION_OF_BIOSYNTHETIC_PROCESS", "GOMF_ATP_DEPENDENT_ACTIVITY",
  #Both
  "GOBP_SECRETION", "GOBP_EXOCYTOSIS", "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOBP_EXPORT_FROM_CELL",
  "GOBP_CELLULAR_GLUCOSE_HOMEOSTASIS", "GOBP_REGULATION_OF_INSULIN_SECRETION", "GOBP_LIPID_LOCALIZATION",
  "GOBP_REGULATION_OF_RESPONSE_TO_STRESS", "GOBP_RESPONSE_TO_CYTOKINE", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
  #Unique to ADd
  "GOCC_ION_CHANNEL_COMPLEX", "GOCC_CATION_CHANNEL_COMPLEX", "GOBP_REGULATION_OF_TRANSPORT", "GOMF_MOLECULAR_TRANSDUCER_ACTIVITY",
  "GOCC_ANCHORING_JUNCTION", "GOCC_CELL_CELL_JUNCTION", "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOBP_CELL_CELL_ADHESION", #"GOMF_CADHERIN_BINDING"
  "GOBP_INTRACILIARY_TRANSPORT", "GOCC_DYNEIN_COMPLEX", "GOBP_CILIUM_MOVEMENT", "GOBP_CEREBROSPINAL_FLUID_CIRCULATION"
)

##per contrast
GSEA_Epithelial_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Epithelial, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_Epithelial, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_ADdvresilient_Epithelial, by="pathway", suffix=c("", "_ADdvresilient"))

GSEA_Epithelial_join <- GSEA_Epithelial_join %>% column_to_rownames("pathway")
GSEA_Epithelial_join <- GSEA_Epithelial_join[selection,]
GSEA_Epithelial_join <- GSEA_Epithelial_join[match(selection, rownames(GSEA_Epithelial_join)),]
mat <- GSEA_Epithelial_join[,c("NES", "NES_adnc", "NES_ADdvresilient")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_Epithelial_join[,c("padj", "padj_adnc", "padj_ADdvresilient")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = colorRamp2(c(-2, 0, 2), c("blue","white", "red")), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Epithelial_Res.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_ADvResilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADdvresilient), col=col_fun(mat$ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
dev.off()


###Fib

selection <- c(#unique to resilience
  "GOBP_PEPTIDE_METABOLIC_PROCESS", "GOBP_MACROMOLECULE_CATABOLIC_PROCESS", "GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOMF_RRNA_BINDING", #"GOBP_AMIDE_BIOSYNTHETIC_PROCESS",
  "GOCC_CELL_SUBSTRATE_JUNCTION", "GOCC_ANCHORING_JUNCTION",
  "GOBP_TIGHT_JUNCTION_ORGANIZATION", "GOCC_APICAL_JUNCTION_COMPLEX", #"GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
  "GOBP_T_CELL_ACTIVATION", "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_REGULATION_OF_IMMUNE_RESPONSE", "GOBP_CELL_AGING",
  "GOCC_CATION_CHANNEL_COMPLEX", "GOCC_POTASSIUM_CHANNEL_COMPLEX", "GOBP_REGULATION_OF_ION_TRANSMEMBRANE_TRANSPORT", "GOMF_PASSIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY"
)

#Fibroblast
GSEA_Fibroblast_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_ADvMCI")) %>% full_join(GSEA_adnc_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_ADdvresilient_Fibroblast, by="pathway", suffix=c("", "_ADdvresilient"))

GSEA_Fibroblast_join <- GSEA_Fibroblast_join %>% column_to_rownames("pathway")
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[selection,]
GSEA_Fibroblast_join <- GSEA_Fibroblast_join[match(selection, rownames(GSEA_Fibroblast_join)),]
mat <- GSEA_Fibroblast_join[,c("NES", "NES_adnc", "NES_ADdvresilient")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_Epithelial_join[,c("padj", "padj_adnc", "padj_ADdvresilient")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)


plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = colorRamp2(c(-2, 0, 2), c("blue","white", "red")), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Fibroblast_Res.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_ADvResilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADdvresilient), col=col_fun(mat$ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
dev.off()


###BAM

selection <- c(#Both
  "GOMF_UNFOLDED_PROTEIN_BINDING", "GOMF_HSP90_PROTEIN_BINDING", "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_FATTY_ACID_DERIVATIVE_BIOSYNTHETIC_PROCESS", "GOBP_FATTY_ACID_DERIVATIVE_METABOLIC_PROCESS", "GOBP_MULTICELLULAR_ORGANISMAL_HOMEOSTASIS", "GOMF_PHOSPHATIDYLINOSITOL_3_KINASE_ACTIVITY")

##per contrast
GSEA_BAM_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_BAM, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_BAM, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_ADdvresilient_Immune, by="pathway", suffix=c("", "_ADdvresilient"))

GSEA_BAM_join <- GSEA_BAM_join %>% column_to_rownames("pathway")
GSEA_BAM_join <- GSEA_BAM_join[selection,]
GSEA_BAM_join <- GSEA_BAM_join[match(selection, rownames(GSEA_BAM_join)),]
mat <- GSEA_BAM_join[,c("NES", "NES_adnc", "NES_ADdvresilient")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_Epithelial_join[,c("padj", "padj_adnc", "padj_ADdvresilient")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)


plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = colorRamp2(c(-2, 0, 2), c("blue","white", "red")), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Immune_Res.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Immune_ADvResilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADdvresilient), col=col_fun(mat$ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
dev.off()


###Endo

selection <- c("GOMF_RRNA_BINDING", "GOBP_AMIDE_BIOSYNTHETIC_PROCESS", "GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOBP_PEPTIDE_METABOLIC_PROCESS", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS")

##per contrast
GSEA_Endo_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Endo, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_Endo, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_ADdvresilient_Endo, by="pathway", suffix=c("", "_ADdvresilient"))

GSEA_Endo_join <- GSEA_Endo_join %>% column_to_rownames("pathway")
GSEA_Endo_join <- GSEA_Endo_join[selection,]
GSEA_Endo_join <- GSEA_Endo_join[match(selection, rownames(GSEA_Endo_join)),]
mat <- GSEA_Endo_join[,c("NES", "NES_adnc", "NES_ADdvresilient")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_Epithelial_join[,c("padj", "padj_adnc", "padj_ADdvresilient")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)


plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = colorRamp2(c(-2, 0, 2), c("blue","white", "red")), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Endothelial_Res.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Endothelial_ADvResilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADdvresilient), col=col_fun(mat$ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
dev.off()


###Mural

selection <- c("GOBP_EXTRACELLULAR_TRANSPORT", "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOMF_INORGANIC_MOLECULAR_ENTITY_TRANSMEMBRANE_TRANSPORTER_ACTIVITY")

##per contrast
GSEA_Mural_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Mural, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_adnc_ADvNCI_Mural, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_ADdvresilient_Mural, by="pathway", suffix=c("", "_ADdvresilient"))

GSEA_Mural_join <- GSEA_Mural_join %>% column_to_rownames("pathway")
GSEA_Mural_join <- GSEA_Mural_join[selection,]
GSEA_Mural_join <- GSEA_Mural_join[match(selection, rownames(GSEA_Mural_join)),]
mat <- GSEA_Mural_join[,c("NES", "NES_adnc", "NES_ADdvresilient")]
colnames(mat) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_Epithelial_join[,c("padj", "padj_adnc", "padj_ADdvresilient")]
colnames(mat2) <- c("cogdx_ADvNCI", "ADNC_ADvNCI", "ADdvresilient")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)


plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = colorRamp2(c(-2, 0, 2), c("blue","white", "red")), cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Mural_Res.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Mural_ADvResilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat$ADdvresilient), col=col_fun(mat$ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
dev.off()

