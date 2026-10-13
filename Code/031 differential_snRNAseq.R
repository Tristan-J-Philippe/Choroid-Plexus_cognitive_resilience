###
### 031 differential_RNA.R
###
# Purpose: differential in single nuclei RNA dataset
# DependeCtrles:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(42)

#Set and load dependent variables for limma/DESeq functions
Gene.search <- read.csv(paste0(curr_dir, "/DEG to GO.csv"))
FDRcutoff <- 0.05
LFcutoff <- 0.25
pct_cutoff <- 0.05

#metadata
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly.csv"))
mdat <- column_to_rownames(mdat, "projid")

##Subset
###Final cell type
Ex_ChP <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS"))
#Subcelltype
datprepfun_lane(Ex_ChP, "Epi_1a")
datprepfun_lane(Ex_ChP, "Epi_1b")
datprepfun_lane(Ex_ChP, "Epi_2a")
datprepfun_lane(Ex_ChP, "Epi_2b")
datprepfun_lane(Ex_ChP, "Fib_1")
datprepfun_lane(Ex_ChP, "Fib_2")
datprepfun_lane(Ex_ChP, "Fib_3")
datprepfun_lane(Ex_ChP, "Endo_1")
datprepfun_lane(Ex_ChP, "Endo_3")
datprepfun_lane(Ex_ChP, "Endo_2")
datprepfun_lane(Ex_ChP, "Mural_1")
datprepfun_lane(Ex_ChP, "Mural_2")
datprepfun_lane(Ex_ChP, "Mural_3")
datprepfun_lane(Ex_ChP, "BAM_1")
datprepfun_lane(Ex_ChP, "BAM_2")
datprepfun_lane(Ex_ChP, "BAM_3")
datprepfun_lane(Ex_ChP, "T_Cell")
#datprepfun_lane(Ex_ChP, "BAM")

datprepfun_lane(Ex_ChP, "Neuro")
datprepfun_lane(Ex_ChP, "Astro")
datprepfun_lane(Ex_ChP, "Oligo")


###Major cell type

#!!!!Change to Majorcelltype in datprepfun_lane
datprepfun_lane(Ex_ChP, "Epithelial")
datprepfun_lane(Ex_ChP, "Fibroblast")
datprepfun_lane(Ex_ChP, "Endothelial")
datprepfun_lane(Ex_ChP, "Mural")
datprepfun_lane(Ex_ChP, "Immune")

datprepfun_lane(Ex_ChP, "Parenchyma")

pcts <- cbind(pct_Epi_1a, pct_Epi_1b, pct_Epi_2a, pct_Epi_2b, pct_Fib_1, pct_Fib_2, pct_Fib_3, pct_Endo_1, pct_Endo_2, pct_Endo_3, pct_Mural_1, pct_Mural_2, pct_Mural_3, pct_BAM_1, pct_BAM_2, pct_BAM_3, pct_T_Cell, pct_Epithelial, pct_Fibroblast, pct_Endothelial, pct_Mural, pct_Immune)
colnames(pcts) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b", "Fib_1", "Fib_2", "Fib_3", "Endo_1", "Endo_2", "Endo_3", "Mural_1", "Mural_2", "Mural_3", "BAM_1", "BAM_2", "BAM_3", "T_Cell", "Epithelial", "Fibroblast", "Endothelial", "Mural", "Immune")
write.csv(pcts, paste0(curr_dir, "/Percent_genes_in_cells.csv"))


##limma runs
###factors set up

#discrete var
msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))

#continuous
projid = as.numeric(mdat$projid)
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)

####Resilience

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path + msex + age + X10X_batch + pmi)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdxadnc_resilientvCtrl=Cog_Pathresilient-Cog_PathCtrl, cogdxadnc_ADdvresilient=Cog_PathADd-Cog_Pathresilient,  cogdxadnc_ADdvCtrl=Cog_PathADd-Cog_PathCtrl)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_Immune, pct_Immune)

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path + msex + age)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdxadnc_resilientvCtrl=Cog_Pathresilient-Cog_PathCtrl, cogdxadnc_ADdvresilient=Cog_PathADd-Cog_Pathresilient,  cogdxadnc_ADdvCtrl=Cog_PathADd-Cog_PathCtrl)
limmafun_lane(cts_Mural, pct_Mural)




#Pattern/Template matching

dir.create(paste0(curr_dir, "/limma/Pattern"))
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly.csv"))
mdat <- column_to_rownames(mdat, "projid")

#Make factors
Cog_Path = factor(paste(mdat$Cog_Path))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path)
rownames(Maindm) <- rownames(mdat)
Maindm <- Maindm[complete.cases(Maindm),]




##Epithelial sort

#Get data
limma_cogdxadnc_resilientvCtrl_Epithelial <- limma_cogdxadnc_resilientvCtrl_Epithelial[!duplicated(limma_cogdxadnc_resilientvCtrl_Epithelial$gene),]
limma_resilientvAD_Epithelial <- limma_resilientvAD_Epithelial[!duplicated(limma_resilientvAD_Epithelial$gene),]
limma_cogdxadnc_ADdvCtrl_Epithelial <- limma_cogdxadnc_ADdvCtrl_Epithelial[!duplicated(limma_cogdxadnc_ADdvCtrl_Epithelial$gene),]

SN <- limma_cogdxadnc_ADdvCtrl_Epithelial[,c(1,9:(length(limma_cogdxadnc_ADdvCtrl_Epithelial)-2))]
rownames(SN) <- NULL
SN <- column_to_rownames(SN, "gene")

SN <- data.frame(t(SN))
temp <- list(substr(rownames(SN),1,9))
SN <- aggregate(SN, temp, mean)
SN <- column_to_rownames(SN, "Group.1")
SN <- data.frame(t(SN))
rownames(SN) <- limma_cogdxadnc_ADdvCtrl_Epithelial$gene #put correct rownames, aggregate changes '-' to '.', thanks but no thanks?
colnames(SN) <- substr(colnames(SN),2,9)
SN <- SN[,colnames(SN) %in% rownames(mdat)]

vdatweights <- SN

#Templates
dm <- Maindm[rownames(Maindm) %in% colnames(vdatweights),]
dm <- data.frame(dm)
Temp_resilient <- dm$Cog_Pathresilient
Temp_ADd <- dm$Cog_PathADd
Temp_Both <- dm$Cog_PathCtrl


###TM

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.
SN_cogdxadnc_resilientvCtrl <- data.frame(limma_cogdxadnc_resilientvCtrl_Epithelial$gene, limma_cogdxadnc_resilientvCtrl_Epithelial$P.Value, limma_cogdxadnc_resilientvCtrl_Epithelial$logFC)
colnames(SN_cogdxadnc_resilientvCtrl) <- c("gene", "P.Value_cogdxadnc_resilientvCtrl", "LFC_cogdxadnc_resilientvCtrl")

SN_resilientvAD <- data.frame(limma_resilientvAD_Epithelial$gene, limma_resilientvAD_Epithelial$P.Value, limma_resilientvAD_Epithelial$logFC)
colnames(SN_resilientvAD) <- c("gene", "P.Value_resilientvAD", "LFC_resilientvAD")

SN_cogdxadnc_ADdvCtrl <- data.frame(limma_cogdxadnc_ADdvCtrl_Epithelial$gene, limma_cogdxadnc_ADdvCtrl_Epithelial$P.Value, limma_cogdxadnc_ADdvCtrl_Epithelial$logFC)
colnames(SN_cogdxadnc_ADdvCtrl) <- c("gene", "P.Value_cogdxadnc_ADdvCtrl", "LFC_cogdxadnc_ADdvCtrl")

stvals <- merge(merge(SN_cogdxadnc_resilientvCtrl, SN_resilientvAD), SN_cogdxadnc_ADdvCtrl)
stvals <- column_to_rownames(stvals, "gene")
stvals <- t(stvals)

rownames(stvals)

Epithelial_TM_resilient <- sort(stvals[1, stvals[1,]<Pvalcutoff & stvals[3,]<Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[4,])>LFcutoff])
Cor.fun(TM_select = Epithelial_TM_resilient, Template = Temp_resilient, tval = data.frame(limma_Pat_resilientvsCtrlandADd_Epithelial$gene, limma_Pat_resilientvsCtrlandADd_Epithelial$t))

Epithelial_TM_ADd <- sort(stvals[5, stvals[3,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])>LFcutoff])
Cor.fun(TM_select = Epithelial_TM_ADd, Template = Temp_ADd, tval = data.frame(limma_Pat_ADdvsresilientandCtrl_Epithelial$gene, limma_Pat_ADdvsresilientandCtrl_Epithelial$t))

Epithelial_TM_both <- sort(stvals[5, stvals[1,]<Pvalcutoff & stvals[5,]<Pvalcutoff & stvals[3,]>Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[6,])>LFcutoff & abs(stvals[4,])<LFcutoff])
Cor.fun(TM_select = Epithelial_TM_both, Template = Temp_Both, tval = data.frame(limma_Pat_resilientandADdvCtrl_Epithelial$gene, limma_Pat_resilientandADdvCtrl_Epithelial$t))


###Resilient

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Epithelial_TM_resilient_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("ADd")

rowMeans <- cbind(rowMeansresilient, rowMeansCtrl, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$resilient),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Epithelial resilient.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Epithelial_TM_resilient_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Epithelial_TM_resilient_GO <- as.data.frame(GO_results)
View(Pattern_Epithelial_TM_resilient_GO)
View(Pattern_Epithelial_TM_resilient_sig)


###ADd

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Epithelial_TM_ADd_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("ADd")

rowMeans <- cbind(rowMeansCtrl, rowMeansresilient, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Epithelial ADd.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Epithelial_TM_ADd_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Epithelial_TM_ADd_GO <- as.data.frame(GO_results)
View(Pattern_Epithelial_TM_ADd_GO)
View(Pattern_Epithelial_TM_ADd_sig)


###Both

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Epithelial_TM_both_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansCtrl, rowMeansresilient, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$Ctrl),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Epithelial Both.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Epithelial_TM_both_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1)
Pattern_Epithelial_TM_both_GO <- as.data.frame(GO_results)
View(Pattern_Epithelial_TM_both_GO)
View(Pattern_Epithelial_TM_both_sig)



##Fibroblast sort

#Get data
limma_cogdxadnc_resilientvCtrl_Fibroblast <- limma_cogdxadnc_resilientvCtrl_Fibroblast[!duplicated(limma_cogdxadnc_resilientvCtrl_Fibroblast$gene),]
limma_resilientvAD_Fibroblast <- limma_resilientvAD_Fibroblast[!duplicated(limma_resilientvAD_Fibroblast$gene),]
limma_cogdxadnc_ADdvCtrl_Fibroblast <- limma_cogdxadnc_ADdvCtrl_Fibroblast[!duplicated(limma_cogdxadnc_ADdvCtrl_Fibroblast$gene),]

SN <- limma_cogdxadnc_ADdvCtrl_Fibroblast[,c(1,9:(length(limma_cogdxadnc_ADdvCtrl_Fibroblast)-2))]
rownames(SN) <- NULL
SN <- column_to_rownames(SN, "gene")

SN <- data.frame(t(SN))
temp <- list(substr(rownames(SN),1,9))
SN <- aggregate(SN, temp, mean)
SN <- column_to_rownames(SN, "Group.1")
SN <- data.frame(t(SN))
rownames(SN) <- limma_cogdxadnc_ADdvCtrl_Fibroblast$gene #put correct rownames, aggregate changes '-' to '.', thanks but no thanks?
colnames(SN) <- substr(colnames(SN),2,9)
SN <- SN[,colnames(SN) %in% rownames(mdat)]

vdatweights <- SN

#Templates
dm <- Maindm[rownames(Maindm) %in% colnames(vdatweights),]
dm <- data.frame(dm)
Temp_resilient <- dm$Cog_Pathresilient
Temp_ADd <- dm$Cog_PathADd
Temp_Both <- dm$Cog_PathCtrl


###TM

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.
SN_cogdxadnc_resilientvCtrl <- data.frame(limma_cogdxadnc_resilientvCtrl_Fibroblast$gene, limma_cogdxadnc_resilientvCtrl_Fibroblast$P.Value, limma_cogdxadnc_resilientvCtrl_Fibroblast$logFC)
colnames(SN_cogdxadnc_resilientvCtrl) <- c("gene", "P.Value_cogdxadnc_resilientvCtrl", "LFC_cogdxadnc_resilientvCtrl")

SN_resilientvAD <- data.frame(limma_resilientvAD_Fibroblast$gene, limma_resilientvAD_Fibroblast$P.Value, limma_resilientvAD_Fibroblast$logFC)
colnames(SN_resilientvAD) <- c("gene", "P.Value_resilientvAD", "LFC_resilientvAD")

SN_cogdxadnc_ADdvCtrl <- data.frame(limma_cogdxadnc_ADdvCtrl_Fibroblast$gene, limma_cogdxadnc_ADdvCtrl_Fibroblast$P.Value, limma_cogdxadnc_ADdvCtrl_Fibroblast$logFC)
colnames(SN_cogdxadnc_ADdvCtrl) <- c("gene", "P.Value_cogdxadnc_ADdvCtrl", "LFC_cogdxadnc_ADdvCtrl")

stvals <- merge(merge(SN_cogdxadnc_resilientvCtrl, SN_resilientvAD), SN_cogdxadnc_ADdvCtrl)
stvals <- column_to_rownames(stvals, "gene")
stvals <- t(stvals)

rownames(stvals)

Fibroblast_TM_resilient <- sort(stvals[1, stvals[1,]<Pvalcutoff & stvals[3,]<Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[4,])>LFcutoff])
Cor.fun(TM_select = Fibroblast_TM_resilient, Template = Temp_resilient, tval = data.frame(limma_Pat_resilientvsCtrlandADd_Fibroblast$gene, limma_Pat_resilientvsCtrlandADd_Fibroblast$t))

Fibroblast_TM_ADd <- sort(stvals[5, stvals[3,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])>LFcutoff])
Cor.fun(TM_select = Fibroblast_TM_ADd, Template = Temp_ADd, tval = data.frame(limma_Pat_ADdvsresilientandCtrl_Fibroblast$gene, limma_Pat_ADdvsresilientandCtrl_Fibroblast$t))

Fibroblast_TM_both <- sort(stvals[5, stvals[1,]<Pvalcutoff & stvals[5,]<Pvalcutoff & stvals[3,]>Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[6,])>LFcutoff & abs(stvals[4,])<LFcutoff])
Cor.fun(TM_select = Fibroblast_TM_both, Template = Temp_Both, tval = data.frame(limma_Pat_resilientandADdvCtrl_Fibroblast$gene, limma_Pat_resilientandADdvCtrl_Fibroblast$t))


###Resilient

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Fibroblast_TM_resilient_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("ADd")

rowMeans <- cbind(rowMeansresilient, rowMeansCtrl, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$resilient),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Fibroblast resilient.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Fibroblast_TM_resilient_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Fibroblast_TM_resilient_GO <- as.data.frame(GO_results)
View(Pattern_Fibroblast_TM_resilient_GO)
View(Pattern_Fibroblast_TM_resilient_sig)


###ADd

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Fibroblast_TM_ADd_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("ADd")

rowMeans <- cbind(rowMeansCtrl, rowMeansresilient, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Fibroblast ADd.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Fibroblast_TM_ADd_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Fibroblast_TM_ADd_GO <- as.data.frame(GO_results)
View(Pattern_Fibroblast_TM_ADd_GO)
View(Pattern_Fibroblast_TM_ADd_sig)


###Both

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Fibroblast_TM_both_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansCtrl, rowMeansresilient, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$Ctrl),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Fibroblast Both.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Fibroblast_TM_both_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Fibroblast_TM_both_GO <- as.data.frame(GO_results)
View(Pattern_Fibroblast_TM_both_GO)
View(Pattern_Fibroblast_TM_both_sig)



##Immune sort

#Get data
limma_cogdxadnc_resilientvCtrl_Immune <- limma_cogdxadnc_resilientvCtrl_Immune[!duplicated(limma_cogdxadnc_resilientvCtrl_Immune$gene),]
limma_resilientvAD_Immune <- limma_resilientvAD_Immune[!duplicated(limma_resilientvAD_Immune$gene),]
limma_cogdxadnc_ADdvCtrl_Immune <- limma_cogdxadnc_ADdvCtrl_Immune[!duplicated(limma_cogdxadnc_ADdvCtrl_Immune$gene),]

SN <- limma_cogdxadnc_ADdvCtrl_Immune[,c(1,9:(length(limma_cogdxadnc_ADdvCtrl_Immune)-2))]
rownames(SN) <- NULL
SN <- column_to_rownames(SN, "gene")

SN <- data.frame(t(SN))
temp <- list(substr(rownames(SN),1,9))
SN <- aggregate(SN, temp, mean)
SN <- column_to_rownames(SN, "Group.1")
SN <- data.frame(t(SN))
rownames(SN) <- limma_cogdxadnc_ADdvCtrl_Immune$gene #put correct rownames, aggregate changes '-' to '.', thanks but no thanks?
colnames(SN) <- substr(colnames(SN),2,9)
SN <- SN[,colnames(SN) %in% rownames(mdat)]

vdatweights <- SN

#Templates
dm <- Maindm[rownames(Maindm) %in% colnames(vdatweights),]
dm <- data.frame(dm)
Temp_resilient <- dm$Cog_Pathresilient
Temp_ADd <- dm$Cog_PathADd
Temp_Both <- dm$Cog_PathCtrl


###TM

#pvalues from extracted comparisons after applying top table and sorting alphabetically into one sheet.
SN_cogdxadnc_resilientvCtrl <- data.frame(limma_cogdxadnc_resilientvCtrl_Immune$gene, limma_cogdxadnc_resilientvCtrl_Immune$P.Value, limma_cogdxadnc_resilientvCtrl_Immune$logFC)
colnames(SN_cogdxadnc_resilientvCtrl) <- c("gene", "P.Value_cogdxadnc_resilientvCtrl", "LFC_cogdxadnc_resilientvCtrl")

SN_resilientvAD <- data.frame(limma_resilientvAD_Immune$gene, limma_resilientvAD_Immune$P.Value, limma_resilientvAD_Immune$logFC)
colnames(SN_resilientvAD) <- c("gene", "P.Value_resilientvAD", "LFC_resilientvAD")

SN_cogdxadnc_ADdvCtrl <- data.frame(limma_cogdxadnc_ADdvCtrl_Immune$gene, limma_cogdxadnc_ADdvCtrl_Immune$P.Value, limma_cogdxadnc_ADdvCtrl_Immune$logFC)
colnames(SN_cogdxadnc_ADdvCtrl) <- c("gene", "P.Value_cogdxadnc_ADdvCtrl", "LFC_cogdxadnc_ADdvCtrl")

stvals <- merge(merge(SN_cogdxadnc_resilientvCtrl, SN_resilientvAD), SN_cogdxadnc_ADdvCtrl)
stvals <- column_to_rownames(stvals, "gene")
stvals <- t(stvals)

rownames(stvals)

Immune_TM_resilient <- sort(stvals[1, stvals[1,]<Pvalcutoff & stvals[3,]<Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[4,])>LFcutoff])
Cor.fun(TM_select = Immune_TM_resilient, Template = Temp_resilient, tval = data.frame(limma_Pat_resilientvsCtrlandADd_Immune$gene, limma_Pat_resilientvsCtrlandADd_Immune$t))

Immune_TM_ADd <- sort(stvals[5,  stvals[3,]<Pvalcutoff & stvals[5,]<Pvalcutoff & abs(stvals[4,])>LFcutoff & abs(stvals[6,])>LFcutoff])
Cor.fun(TM_select = Immune_TM_ADd, Template = Temp_ADd, tval = data.frame(limma_Pat_ADdvsresilientandCtrl_Immune$gene, limma_Pat_ADdvsresilientandCtrl_Immune$t))

Immune_TM_both <- sort(stvals[5, stvals[1,]<Pvalcutoff & stvals[5,]<Pvalcutoff & stvals[3,]>Pvalcutoff & abs(stvals[2,])>LFcutoff & abs(stvals[6,])>LFcutoff & abs(stvals[4,])<LFcutoff])
Cor.fun(TM_select = Immune_TM_both, Template = Temp_Both, tval = data.frame(limma_Pat_resilientandADdvCtrl_Immune$gene, limma_Pat_resilientandADdvCtrl_Immune$t))


###Resilient

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Immune_TM_resilient_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("ADd")

rowMeans <- cbind(rowMeansresilient, rowMeansCtrl, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$resilient),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Immune resilient.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Immune_TM_resilient_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Immune_TM_resilient_GO <- as.data.frame(GO_results)
View(Pattern_Immune_TM_resilient_GO)
View(Pattern_Immune_TM_resilient_sig)


###ADd

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Immune_TM_ADd_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("ADd", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("ADd")

rowMeans <- cbind(rowMeansCtrl, rowMeansresilient, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$AD),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Immune ADd.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Immune_TM_ADd_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Immune_TM_ADd_GO <- as.data.frame(GO_results)
View(Pattern_Immune_TM_ADd_GO)
View(Pattern_Immune_TM_ADd_sig)


###Both

SN_mdat <- data.frame(t(SN))
mdat_temp <- data.frame(rownames(mdat), mdat$Cog_Path)
mdat_temp <- column_to_rownames(mdat_temp, "rownames.mdat.")
SN_mdat <- merge(SN_mdat, mdat_temp, by='row.names')
SN_mdat$Row.names <- NULL
SN_mdat <- data.frame(t(SN_mdat))
colnames(SN_mdat) <- SN_mdat["mdat.Cog_Path",]
SN_mdat <- SN_mdat[-which(rownames(SN_mdat) == "mdat.Cog_Path"),]
SN_mdat <- SN_mdat[,order(colnames(SN_mdat), decreasing=T)]

#select
genes <- unique(c(Pattern_Immune_TM_both_sig$gene))

SN_mdat <- SN_mdat[which(rownames(SN_mdat) %in% genes),]
SN_mdatn <- apply(SN_mdat, 2, as.numeric)
rownames(SN_mdatn) <- rownames(SN_mdat)
Sc <- t(scale(t(SN_mdatn)))

rowMeans <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeans) <- c("Means")
rowMeans <- rowMeans[order(-rowMeans$Means),, drop = FALSE] #Order row means based on value
Sc <- Sc[order(match(rownames(Sc), rownames(rowMeans))),] #Order data based on row means

Heatmap(Sc, cluster_columns = F, cluster_rows=F)

rowMeansCtrl <- data.frame(rowMeans(Sc[,c(grep("Ctrl", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansCtrl) <- c("Ctrl")
rowMeansresilient <- data.frame(rowMeans(Sc[,c(grep("resilient", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansresilient) <- c("resilient")
rowMeansAD <- data.frame(rowMeans(Sc[,c(grep("AD", colnames(Sc)))])) #Row means match group of interest
colnames(rowMeansAD) <- c("AD")

rowMeans <- cbind(rowMeansCtrl, rowMeansresilient, rowMeansAD)
rowMeans <- rowMeans[order(-rowMeans$Ctrl),, drop = FALSE] #Order row means based on value
ht <- Heatmap(rowMeans, cluster_columns = F, cluster_rows=F)
graph2pdf(ht, paste0(curr_dir, "/limma/Pattern/Heatmap Immune Both.pdf"))
plot(ht)

GO_results <- enrichGO(gene = c(Pattern_Immune_TM_both_sig$gene), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
Pattern_Immune_TM_both_GO <- as.data.frame(GO_results)
View(Pattern_Immune_TM_both_GO)
View(Pattern_Immune_TM_both_sig)






#FGSEA
#vCtrl
Path_fun(limma_cogdxadnc_resilientvCtrl_Epi_1a)
Path_fun(limma_cogdxadnc_resilientvCtrl_Epi_1b)
Path_fun(limma_cogdxadnc_resilientvCtrl_Epi_2a)
Path_fun(limma_cogdxadnc_resilientvCtrl_Epi_2b)
Path_fun(limma_cogdxadnc_resilientvCtrl_Fib_1)
Path_fun(limma_cogdxadnc_resilientvCtrl_Fib_2)
Path_fun(limma_cogdxadnc_resilientvCtrl_Fib_3)
Path_fun(limma_cogdxadnc_resilientvCtrl_BAM_1)
Path_fun(limma_cogdxadnc_resilientvCtrl_BAM_2)

Path_fun(limma_cogdxadnc_resilientvCtrl_Mural)
Path_fun(limma_cogdxadnc_resilientvCtrl_Epithelial)
Path_fun(limma_cogdxadnc_resilientvCtrl_Fibroblast)
Path_fun(limma_cogdxadnc_resilientvCtrl_Endothelial)
Path_fun(limma_cogdxadnc_resilientvCtrl_Immune)

Path_fun(limma_cogdxadnc_resilientvCtrl_Neuro)
Path_fun(limma_cogdxadnc_resilientvCtrl_Oligo)

#cogdxadnc_ADdvCtrl
Path_fun(limma_cogdxadnc_ADdvCtrl_Epi_1a)
Path_fun(limma_cogdxadnc_ADdvCtrl_Epi_1b)
Path_fun(limma_cogdxadnc_ADdvCtrl_Epi_2a)
Path_fun(limma_cogdxadnc_ADdvCtrl_Epi_2b)
Path_fun(limma_cogdxadnc_ADdvCtrl_Fib_1)
Path_fun(limma_cogdxadnc_ADdvCtrl_Fib_2)
Path_fun(limma_cogdxadnc_ADdvCtrl_Fib_3)

Path_fun(limma_cogdxadnc_ADdvCtrl_Mural)
Path_fun(limma_cogdxadnc_ADdvCtrl_Epithelial)
Path_fun(limma_cogdxadnc_ADdvCtrl_Fibroblast)
Path_fun(limma_cogdxadnc_ADdvCtrl_Endothelial)
Path_fun(limma_cogdxadnc_ADdvCtrl_Immune)

Path_fun(limma_cogdxadnc_ADdvresilient_Neuro)
Path_fun(limma_cogdxadnc_ADdvresilient_Oligo)

#cogdxadnc_ADdvresilient
Path_fun(limma_cogdxadnc_ADdvresilient_Epi_1a)
Path_fun(limma_cogdxadnc_ADdvresilient_Epi_1b)
Path_fun(limma_cogdxadnc_ADdvresilient_Epi_2a)
Path_fun(limma_cogdxadnc_ADdvresilient_Epi_2b)
Path_fun(limma_cogdxadnc_ADdvresilient_Fib_1)
Path_fun(limma_cogdxadnc_ADdvresilient_Fib_2)
Path_fun(limma_cogdxadnc_ADdvresilient_Fib_3)

Path_fun(limma_cogdxadnc_ADdvresilient_Mural)
Path_fun(limma_cogdxadnc_ADdvresilient_Epithelial)
Path_fun(limma_cogdxadnc_ADdvresilient_Fibroblast)
Path_fun(limma_cogdxadnc_ADdvresilient_Endothelial)
Path_fun(limma_cogdxadnc_ADdvresilient_Immune)

Path_fun(limma_cogdxadnc_ADdvCtrl_Neuro)
Path_fun(limma_cogdxadnc_ADdvCtrl_Oligo)






library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
library(export)
col_fun <-  colorRamp2(c(-2, 0, 2), c("blue","white", "red"))

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
GSEA_join <- full_join(GO_labelled, GSEA_cogdxadnc_resilientvCtrl_Epithelial, by="pathway", suffix=c("", "_cogdxadnc_resilientvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvCtrl_Epithelial, by="pathway", suffix=c("", "_cogdxadnc_ADdvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvresilient_Epithelial, by="pathway", suffix=c("", "_cogdxadnc_ADdvresilient")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_Epi, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_gpath_gpathxCR_gpath_Epi, by="pathway", suffix=c("", "_gpath")) %>% full_join(GSEA_adnc_ADvCtrl_Epithelial, by="pathway", suffix=c("", "_ADNC")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Epi, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_CR_adnc_adncbyCR_CR_Epi, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_adnc_adncbyCR_adnc_Epi, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_CR_adnc_adncbyCR_adncbyCR_Epi, by="pathway", suffix=c("", "_CRxadnc"))

###Epithelial
GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("cogdxadnc_resilientvCtrl", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("cogdxadnc_resilientvCtrl", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Resilience heatmap Epithelial.pdf", sep = ""), font="Arial", width=w, height=h)

pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_cogdxadnc_resilientvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_resilientvCtrl), col=col_fun(mat$cogdxadnc_resilientvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_cogdxadnc_ADdvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvCtrl), col=col_fun(mat$cogdxadnc_ADdvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_cogdxadnc_ADdvresilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvresilient), col=col_fun(mat$cogdxadnc_ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial_adnc.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$ADNC), col=col_fun(mat$ADNC), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

###Fibroblast
selection <- c(#unique to resilience
  "GOBP_PEPTIDE_METABOLIC_PROCESS", "GOBP_MACROMOLECULE_CATABOLIC_PROCESS", "GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOMF_RRNA_BINDING", #"GOBP_AMIDE_BIOSYNTHETIC_PROCESS",
  "GOCC_CELL_SUBSTRATE_JUNCTION", "GOCC_ANCHORING_JUNCTION",
  "GOBP_TIGHT_JUNCTION_ORGANIZATION", "GOCC_APICAL_JUNCTION_COMPLEX", #"GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
  "GOBP_T_CELL_ACTIVATION", "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_REGULATION_OF_IMMUNE_RESPONSE", "GOBP_CELL_AGING",
  "GOCC_CATION_CHANNEL_COMPLEX", "GOCC_POTASSIUM_CHANNEL_COMPLEX", "GOBP_REGULATION_OF_ION_TRANSMEMBRANE_TRANSPORT", "GOMF_PASSIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY"
)

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxadnc_resilientvCtrl_Fibroblast, by="pathway", suffix=c("", "_cogdxadnc_resilientvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvCtrl_Fibroblast, by="pathway", suffix=c("", "_cogdxadnc_ADdvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvresilient_Fibroblast, by="pathway", suffix=c("", "_cogdxadnc_ADdvresilient")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_Fib, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_gpath_gpathxCR_gpath_Fib, by="pathway", suffix=c("", "_gpath")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Fib, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_adnc_ADvCtrl_Fibroblast, by="pathway", suffix=c("", "_ADNC")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Fib, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_CR_adnc_adncbyCR_CR_Fib, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_adnc_adncbyCR_adnc_Fib, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_CR_adnc_adncbyCR_adncbyCR_Fib, by="pathway", suffix=c("", "_CRxadnc"))
GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("cogdxadnc_resilientvCtrl", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("cogdxadnc_resilientvCtrl", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Resilience heatmap Fibroblast.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_cogdxadnc_resilientvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_resilientvCtrl), col=col_fun(mat$cogdxadnc_resilientvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_cogdxadnc_ADdvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvCtrl), col=col_fun(mat$cogdxadnc_ADdvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_cogdxadnc_ADdvresilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvresilient), col=col_fun(mat$cogdxadnc_ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast_adnc.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$ADNC), col=col_fun(mat$ADNC), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

###Immune
selection <- c(#Both
  "GOMF_UNFOLDED_PROTEIN_BINDING", "GOMF_HSP90_PROTEIN_BINDING", "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_FATTY_ACID_DERIVATIVE_BIOSYNTHETIC_PROCESS", "GOBP_FATTY_ACID_DERIVATIVE_METABOLIC_PROCESS", "GOBP_MULTICELLULAR_ORGANISMAL_HOMEOSTASIS", "GOMF_PHOSPHATIDYLINOSITOL_3_KINASE_ACTIVITY")

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxadnc_resilientvCtrl_Immune, by="pathway", suffix=c("", "_cogdxadnc_resilientvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvCtrl_Immune, by="pathway", suffix=c("", "_cogdxadnc_ADdvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvresilient_Immune, by="pathway", suffix=c("", "_cogdxadnc_ADdvresilient"))%>% full_join(GSEA_CR_gpath_gpathxCR_CR_Imm, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_gpath_gpathxCR_gpath_Imm, by="pathway", suffix=c("", "_gpath")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Imm, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_adnc_ADvCtrl_Immune, by="pathway", suffix=c("", "_ADNC")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Imm, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_CR_adnc_adncbyCR_CR_Imm, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_adnc_adncbyCR_adnc_Imm, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_CR_adnc_adncbyCR_adncbyCR_Imm, by="pathway", suffix=c("", "_CRxadnc"))

GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("cogdxadnc_resilientvCtrl", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("cogdxadnc_resilientvCtrl", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Resilience heatmap Immune.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Immune_cogdxadnc_resilientvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5) #if figure margins too large you have too few bars, set it to 2
plts <- barplot(abs(mat$cogdxadnc_resilientvCtrl), col=col_fun(mat$cogdxadnc_resilientvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Immune_cogdxadnc_ADdvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvCtrl), col=col_fun(mat$cogdxadnc_ADdvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Immune_cogdxadnc_ADdvresilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvresilient), col=col_fun(mat$cogdxadnc_ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Immune_adnc.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$ADNC), col=col_fun(mat$ADNC), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

###Endothelial
selection <- c("GOMF_RRNA_BINDING", "GOBP_AMIDE_BIOSYNTHETIC_PROCESS", "GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOBP_PEPTIDE_METABOLIC_PROCESS", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS")

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxadnc_resilientvCtrl_Endothelial, by="pathway", suffix=c("", "_cogdxadnc_resilientvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvCtrl_Endothelial, by="pathway", suffix=c("", "_cogdxadnc_ADdvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvresilient_Endothelial, by="pathway", suffix=c("", "_cogdxadnc_ADdvresilient")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_Endo, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_gpath_gpathxCR_gpath_Endo, by="pathway", suffix=c("", "_gpath")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Endo, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Endo, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_adnc_ADvCtrl_Endothelial, by="pathway", suffix=c("", "_ADNC")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Endo, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_CR_adnc_adncbyCR_CR_Endo, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_adnc_adncbyCR_adnc_Endo, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_CR_adnc_adncbyCR_adncbyCR_Endo, by="pathway", suffix=c("", "_CRxadnc"))
GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("cogdxadnc_resilientvCtrl_Endothelial", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("cogdxadnc_resilientvCtrl_Endothelial", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Resilience heatmap Endothelial.pdf", sep = ""), font="Arial", width=w, height=h)

pdf(paste(curr_dir, "/FGSEA/Bar_Endothelial_cogdxadnc_resilientvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_resilientvCtrl), col=col_fun(mat$cogdxadnc_resilientvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Endothelial_cogdxadnc_ADdvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvCtrl), col=col_fun(mat$cogdxadnc_ADdvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Endothelial_cogdxadnc_ADdvresilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvresilient), col=col_fun(mat$cogdxadnc_ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Endothelial_adnc.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$ADNC), col=col_fun(mat$ADNC), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

selection <- c("GOBP_EXTRACELLULAR_TRANSPORT", "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOMF_INORGANIC_MOLECULAR_ENTITY_TRANSMEMBRANE_TRANSPORTER_ACTIVITY")

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxadnc_resilientvCtrl_Mural, by="pathway", suffix=c("", "_cogdxadnc_resilientvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvCtrl_Mural, by="pathway", suffix=c("", "_cogdxadnc_ADdvCtrl")) %>% full_join(GSEA_cogdxadnc_ADdvresilient_Mural, by="pathway", suffix=c("", "_cogdxadnc_ADdvresilient")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_Mur, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_gpath_gpathxCR_gpath_Mur, by="pathway", suffix=c("", "_gpath")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Mur, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_adnc_ADvCtrl_Mural, by="pathway", suffix=c("", "_ADNC"))  %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Mur, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_adnc_ADvCtrl_Mural, by="pathway", suffix=c("", "_ADNC")) %>% full_join(GSEA_CR_gpath_gpathxCR_CR_BY_gpath_Mur, by="pathway", suffix=c("", "_CRxgpath")) %>% full_join(GSEA_CR_adnc_adncbyCR_CR_Mur, by="pathway", suffix=c("", "_CR")) %>% full_join(GSEA_CR_adnc_adncbyCR_adnc_Mur, by="pathway", suffix=c("", "_adnc")) %>% full_join(GSEA_CR_adnc_adncbyCR_adncbyCR_Mur, by="pathway", suffix=c("", "_CRxadnc"))
GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("cogdxadnc_resilientvCtrl_Mural", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("cogdxadnc_resilientvCtrl_Mural", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Resilience heatmap Mural.pdf", sep = ""), font="Arial", width=w, height=h)

pdf(paste(curr_dir, "/FGSEA/Bar_Mural_cogdxadnc_resilientvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_resilientvCtrl), col=col_fun(mat$cogdxadnc_resilientvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Mural_cogdxadnc_ADdvCtrl.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvCtrl), col=col_fun(mat$cogdxadnc_ADdvCtrl), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Mural_cogdxadnc_ADdvresilient.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$cogdxadnc_ADdvresilient), col=col_fun(mat$cogdxadnc_ADdvresilient), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()

pdf(paste(curr_dir, "/FGSEA/Bar_Mural_adnc.pdf", sep = ""), width=3, height=nrow(mat)*.5)
plts <- barplot(abs(mat$ADNC), col=col_fun(mat$ADNC), horiz=T, xlim=c(0,2), border=NA)
text(x = 0, y = plts, labels = rownames(mat), pos = 2, font = 7, xpd = TRUE)
dev.off()





#GSVA
###Make modules

#Go term matrices
Temp <- cbind(GSEA_cogdxadnc_resilientvCtrl_Epithelial$pathway, GSEA_cogdxadnc_resilientvCtrl_Epithelial$leadingEdge)
Temp2 <- cbind(GSEA_cogdxadnc_ADdvCtrl_Epithelial$pathway, GSEA_cogdxadnc_ADdvCtrl_Epithelial$leadingEdge)
Temp3 <- cbind(GSEA_cogdxadnc_ADdvresilient_Epithelial$pathway, GSEA_cogdxadnc_ADdvresilient_Epithelial$leadingEdge)
Cilia_terms <- c("GOBP_CILIUM_MOVEMENT", "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_DYNEIN_COMPLEX")
Barrier_terms <- c("GOCC_CELL_CELL_JUNCTION", "GOCC_ANCHORING_JUNCTION", "GOMF_CELL_ADHESION_MOLECULE_BINDING")
Transport_terms <- c("GOCC_ION_CHANNEL_COMPLEX", "GOBP_REGULATION_OF_TRANSPORT", "GOMF_MOLECULAR_TRANSDUCER_ACTIVITY")
Immune_terms <- c("GOBP_RESPONSE_TO_CYTOKINE", "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "GOBP_INFLAMMATORY_RESPONSE")
lipid_terms <- c("GOBP_CELLULAR_GLUCOSE_HOMEOSTASIS", "GOBP_INSULIN_SECRETION", "GOBP_LIPID_LOCALIZATION")
secretion_terms <- c("GOBP_SECRETION", "GOBP_EXPORT_FROM_CELL", "GOBP_VESICLE_MEDIATED_TRANSPORT")
metabolism_terms <- c("GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOCC_CATALYTIC_COMPLEX", "GOBP_NEGATIVE_REGULATION_OF_BIOSYNTHETIC_PROCESS")

SN_Epi_paths <- list("Epi_Cilia"=unique(c(unlist(list(Temp[Temp[,1] %in% Cilia_terms], Temp2[Temp2[,1] %in% Cilia_terms], Temp3[Temp3[,1] %in% Cilia_terms])))),
                     "Epi_barrier"=unique(c(unlist(list(Temp[Temp[,1] %in% Barrier_terms], Temp2[Temp2[,1] %in% Barrier_terms], Temp3[Temp3[,1] %in% Barrier_terms])))),
                     "Epi_Ion_Transport"=unique(c(unlist(list(Temp[Temp[,1] %in% Transport_terms], Temp2[Temp2[,1] %in% Transport_terms], Temp3[Temp3[,1] %in% Transport_terms])))),
                     "Epi_Immune"=unique(c(unlist(list(Temp[Temp[,1] %in% Immune_terms], Temp2[Temp2[,1] %in% Immune_terms], Temp3[Temp3[,1] %in% Immune_terms])))),
                     "Epi_lipid_glucose"=unique(c(unlist(list(Temp[Temp[,1] %in% lipid_terms], Temp2[Temp2[,1] %in% lipid_terms], Temp3[Temp3[,1] %in% lipid_terms])))),
                     "Epi_Secretion"=unique(c(unlist(list(Temp[Temp[,1] %in% secretion_terms], Temp2[Temp2[,1] %in% secretion_terms], Temp3[Temp3[,1] %in% secretion_terms])))),
                     "Epi_protein_metabolism"=unique(c(unlist(list(Temp[Temp[,1] %in% metabolism_terms], Temp2[Temp2[,1] %in% metabolism_terms], Temp3[Temp3[,1] %in% metabolism_terms]))))
)
saveRDS(SN_Epi_paths, paste0(curr_dir, "/GSVA_modules/SN_Epi_modules.rds"))


Temp <- cbind(GSEA_cogdxadnc_resilientvCtrl_Fibroblast$pathway, GSEA_cogdxadnc_resilientvCtrl_Fibroblast$leadingEdge)
Temp2 <- cbind(GSEA_cogdxadnc_ADdvCtrl_Fibroblast$pathway, GSEA_cogdxadnc_ADdvCtrl_Fibroblast$leadingEdge)
Temp3 <- cbind(GSEA_cogdxadnc_ADdvresilient_Fibroblast$pathway, GSEA_cogdxadnc_ADdvresilient_Fibroblast$leadingEdge)
Barrier_terms <- c("GOCC_CELL_SUBSTRATE_JUNCTION", "GOBP_TIGHT_JUNCTION_ORGANIZATION", "GOCC_ANCHORING_JUNCTION")
Transport_terms <- c("GOCC_CATION_CHANNEL_COMPLEX", "GOCC_POTASSIUM_CHANNEL_COMPLEX", "GOMF_PASSIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY")
Immune_terms <- c("GOBP_T_CELL_ACTIVATION", "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_REGULATION_OF_IMMUNE_RESPONSE")
metabolism_terms <- c("GOBP_PEPTIDE_METABOLIC_PROCESS", "GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOBP_MACROMOLECULE_CATABOLIC_PROCESS")

SN_Fib_paths <- list(
  "Fib_Barrier"=unique(c(unlist(list(Temp[Temp[,1] %in% Barrier_terms], Temp2[Temp2[,1] %in% Barrier_terms], Temp3[Temp3[,1] %in% Barrier_terms])))),
  "Fib_Transport"=unique(c(unlist(list(Temp[Temp[,1] %in% Transport_terms], Temp2[Temp2[,1] %in% Transport_terms], Temp3[Temp3[,1] %in% Transport_terms])))),
  "Fib_Immune"=unique(c(unlist(list(Temp[Temp[,1] %in% Immune_terms], Temp2[Temp2[,1] %in% Immune_terms], Temp3[Temp3[,1] %in% Immune_terms])))),
  "Fib_protein_metabolism"=unique(c(unlist(list(Temp[Temp[,1] %in% metabolism_terms], Temp2[Temp2[,1] %in% metabolism_terms], Temp3[Temp3[,1] %in% metabolism_terms]))))
)
saveRDS(SN_Fib_paths, paste0(curr_dir, "/GSVA_modules/SN_Fib_modules.rds"))


Temp <- cbind(GSEA_cogdxadnc_resilientvCtrl_Endothelial$pathway, GSEA_cogdxadnc_resilientvCtrl_Endothelial$leadingEdge)
Temp2 <- cbind(GSEA_cogdxadnc_ADdvCtrl_Endothelial$pathway, GSEA_cogdxadnc_ADdvCtrl_Endothelial$leadingEdge)
Temp3 <- cbind(GSEA_cogdxadnc_ADdvresilient_Endothelial$pathway, GSEA_cogdxadnc_ADdvresilient_Endothelial$leadingEdge)
metabolism_terms <- c("GOBP_AMIDE_BIOSYNTHETIC_PROCESS", "GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS", "GOBP_PEPTIDE_METABOLIC_PROCESS")
SN_Endo_paths <- list("Endo_Metabolism"=unique(c(unlist(list(Temp[Temp[,1] %in% metabolism_terms], Temp2[Temp2[,1] %in% metabolism_terms], Temp3[Temp3[,1] %in% metabolism_terms])))))
saveRDS(SN_Endo_paths, paste0(curr_dir, "/GSVA_modules/SN_Endo_modules.rds"))


Temp <- cbind(GSEA_cogdxadnc_resilientvCtrl_Mural$pathway, GSEA_cogdxadnc_resilientvCtrl_Mural$leadingEdge)
Temp2 <- cbind(GSEA_cogdxadnc_ADdvCtrl_Mural$pathway, GSEA_cogdxadnc_ADdvCtrl_Mural$leadingEdge)
Temp3 <- cbind(GSEA_cogdxadnc_ADdvresilient_Mural$pathway, GSEA_cogdxadnc_ADdvresilient_Mural$leadingEdge)
Transport_terms <- c("GOMF_INORGANIC_MOLECULAR_ENTITY_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_EXTRACELLULAR_TRANSPORT", "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY")
SN_Mural_paths <- list("Mural_Transport"=unique(c(unlist(list(Temp[Temp[,1] %in% Transport_terms], Temp2[Temp2[,1] %in% Transport_terms], Temp3[Temp3[,1] %in% Transport_terms])))))
saveRDS(SN_Mural_paths, paste0(curr_dir, "/GSVA_modules/SN_Mural_modules.rds"))


Temp <- cbind(GSEA_cogdxadnc_resilientvCtrl_Immune$pathway, GSEA_cogdxadnc_resilientvCtrl_Immune$leadingEdge)
Temp2 <- cbind(GSEA_cogdxadnc_ADdvCtrl_Immune$pathway, GSEA_cogdxadnc_ADdvCtrl_Immune$leadingEdge)
Temp3 <- cbind(GSEA_cogdxadnc_ADdvresilient_Immune$pathway, GSEA_cogdxadnc_ADdvresilient_Immune$leadingEdge)
HSP_terms <- c("GOMF_UNFOLDED_PROTEIN_BINDING", "GOMF_HSP90_PROTEIN_BINDING", "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE")
metabolism_terms <- c("GOBP_FATTY_ACID_DERIVATIVE_BIOSYNTHETIC_PROCESS", "GOBP_MULTICELLULAR_ORGANISMAL_HOMEOSTASIS", "GOMF_PHOSPHATIDYLINOSITOL_3_KINASE_ACTIVITY")
SN_Imm_paths <- list("Imm_Metabolism"=unique(c(unlist(list(Temp[Temp[,1] %in% HSP_terms], Temp2[Temp2[,1] %in% HSP_terms], Temp3[Temp3[,1] %in% HSP_terms])))),
                     "Imm_HSP"=unique(c(unlist(list(Temp[Temp[,1] %in% metabolism_terms], Temp2[Temp2[,1] %in% metabolism_terms], Temp3[Temp3[,1] %in% metabolism_terms]))))
)
saveRDS(SN_Imm_paths, paste0(curr_dir, "/GSVA_modules/SN_Imm_modules.rds"))


Paren_mods <- list("Par_HSP" = c("HSPA1A", "HSPA1B", "DNAJB1", "DNAJA4", "HSPA6", "HSPH", "PRKCA", "HSP90AB1", "HSP90AA1"),
                   "Par_Neuron" = c("MEF2C",  "RELN", "ATP8B1", "BDNF", "NTRK2", "SORT1", "NGFR", "ANGPT2", "TEK"),
                   "Par_AstD" = c("AP3B2", "AVIL", "CHI3L2", "GPCPD1", "GPX3", "HMGN2", "LGMN", "NQO1", "ODC1", "PNPLA6", "CHDH"),
                   "Par_AstU" = c("ADORA1", "CAPN2", "CHDH", "FBXO2", "KANK2", "MARCH6", "PHPLAB", "SLC39A11", "GPCPD1"))
saveRDS(Paren_mods, paste0(curr_dir, "/GSVA_modules/Paren_modules.rds"))


###Merge modules from all modalities

SN_Epi_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/SN_Epi_modules.rds"))
ST_Epi_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/ST_Epi_modules.rds"))
Epi_modules <- StratigrapheR::merge_list(SN_Epi_paths, ST_Epi_paths)
saveRDS(Epi_modules, paste0(curr_dir, "/GSVA_modules/combo_Epi_modules_Res.rds"))

SN_Fib_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/SN_Fib_modules.rds"))
ST_Fib_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/ST_Fib_modules.rds"))
Fib_modules <- StratigrapheR::merge_list(SN_Fib_paths, ST_Fib_paths)
saveRDS(Fib_modules, paste0(curr_dir, "/GSVA_modules/combo_Fib_modules_Res.rds"))

SN_Endo_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/SN_Endo_modules.rds"))
ST_Endo_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/ST_Endo_modules.rds"))
Endo_modules <- StratigrapheR::merge_list(SN_Endo_paths, ST_Endo_paths)
saveRDS(Endo_modules, paste0(curr_dir, "/GSVA_modules/combo_Endo_modules_Res.rds"))

SN_Mural_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/SN_Mural_modules.rds"))
ST_Mural_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/ST_Mural_modules.rds"))
Mural_modules <- StratigrapheR::merge_list(SN_Mural_paths, ST_Mural_paths)
saveRDS(Mural_modules, paste0(curr_dir, "/GSVA_modules/combo_Mural_modules_Res.rds"))

SN_Imm_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/SN_Imm_modules.rds"))
ST_Imm_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/ST_Imm_modules.rds"))
Imm_modules <- StratigrapheR::merge_list(SN_Imm_paths, ST_Imm_paths)
saveRDS(Imm_modules, paste0(curr_dir, "/GSVA_modules/combo_Imm_modules_Res.rds"))


###load modules

Epi_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Epi_modules_Res.rds"))
Fib_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Fib_modules_Res.rds"))
Endo_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Endo_modules_Res.rds"))
Mural_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Mural_modules_Res.rds"))
Imm_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Imm_modules_Res.rds"))
Paren_mods <- readRDS(paste0(curr_dir, "/GSVA_modules/Paren_modules.rds"))


##Dataprep

GSVA_fun(limma_cogdxadnc_resilientvCtrl_Epithelial, pathways=c(Epi_modules, Paren_mods))
GSVA_fun(limma_cogdxadnc_resilientvCtrl_Fibroblast, pathways=c(Fib_modules, Paren_mods))
GSVA_fun(limma_cogdxadnc_resilientvCtrl_Endothelial, pathways=c(Endo_modules, Paren_mods))
GSVA_fun(limma_cogdxadnc_resilientvCtrl_Mural, pathways=c(Mural_modules, Paren_mods))
GSVA_fun(limma_cogdxadnc_resilientvCtrl_Immune, pathways=c(Imm_modules, Paren_mods))


###Factors set up

mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly.csv"))
mdat <- column_to_rownames(mdat, "projid")

#discrete var
msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))
library_batch = factor(paste(mdat$library_batch))
Cog_Path = factor(paste(mdat$Cog_Path))

#continuous
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)


####limma Res

options(na.action='na.pass')
Maindm <- model.matrix(~0 + Cog_Path + msex + age + X10X_batch + pmi)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdxadnc_resilientvCtrl=Cog_Pathresilient-Cog_PathCtrl, cogdxadnc_ADdvresilient=Cog_PathADd-Cog_Pathresilient,  cogdxadnc_ADdvCtrl=Cog_PathADd-Cog_PathCtrl)

GSVA_limmafun_lane(GSVA_Epithelial)
GSVA_limmafun_lane(GSVA_Fibroblast)
GSVA_limmafun_lane(GSVA_Endothelial)
GSVA_limmafun_lane(GSVA_Mural)
GSVA_limmafun_lane(GSVA_Immune)


###Load limma ST

limma_GSVA_ST_Epithelial <- read.csv("/Ex_ChP/CosMX2/limma/ADdvresilient/Spreadsheet_GSVA_Epithelial_ADdvresilient_uncorrected.csv", row.names=1)
limma_GSVA_ST_Fibroblast <- read.csv("/Ex_ChP/CosMX2/limma/ADdvresilient/Spreadsheet_GSVA_Fibroblast_ADdvresilient_uncorrected.csv", row.names=1)
limma_GSVA_ST_Endo <- read.csv("/Ex_ChP/CosMX2/limma/ADdvresilient/Spreadsheet_GSVA_Endo_ADdvresilient_uncorrected.csv", row.names=1)
limma_GSVA_ST_Mural <- read.csv("/Ex_ChP/CosMX2/limma/ADdvresilient/Spreadsheet_GSVA_Mural_ADdvresilient_uncorrected.csv", row.names=1)
limma_GSVA_ST_Immune <- read.csv("/Ex_ChP/CosMX2/limma/ADdvresilient/Spreadsheet_GSVA_Immune_ADdvresilient_uncorrected.csv", row.names=1)


#Epithelial
zScat_fun(limma_cogdxadnc_ADdvresilient_GSVA_Epithelial, limma_cogdxadnc_resilientvCtrl_GSVA_Epithelial)
zScat_fun(limma_cogdxadnc_ADdvresilient_GSVA_Epithelial, limma_cogdxadnc_ADdvCtrl_GSVA_Epithelial)
zScat_fun(limma_cogdxadnc_ADdvresilient_GSVA_Epithelial, limma_GSVA_ST_Epithelial)
