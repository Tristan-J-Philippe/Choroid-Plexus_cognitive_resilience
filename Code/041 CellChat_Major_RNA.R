###
### 041 CellChat_Major_RNA.R
###
# Purpose: Produce ligand-receptor interactions and compare groups.
# Dependencies:
source("/040 CellChat_util_RNA.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder

Ex_ChP <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS", sep=""))
Ex_ChP <- subset(Ex_ChP, Majorcelltype!="Parenchyma")
Idents(Ex_ChP) <- Ex_ChP@meta.data$Majorcelltype
levels(Ex_ChP) <- c("Epithelial", "Fibroblast", "Mural", "Endothelial", "Immune")
Idents(Ex_ChP) -> Ex_ChP@meta.data$Majorcelltype


##RunCC
!!!Make sure to change output folder

CCfun(Ex_ChP, "Majorcelltype")
gc()


##RunCC per
###Preprocess

unique(Ex_ChP$Cog_Path)
Ex_ChPNCI <- subset(Ex_ChP, subset=(Cog_Path=="NCI"))
Ex_ChPres <- subset(Ex_ChP, subset=Cog_Path=="resilient")
Ex_ChPAD <- subset(Ex_ChP, subset=Cog_Path=="AD")
rm(Ex_ChP)
gc()



CCfun(Ex_ChPNCI, "Majorcelltype")
rm(Ex_ChPNCI)
gc()
CCfun(Ex_ChPres, "Majorcelltype")
rm(Ex_ChPres)
gc()
CCfun(Ex_ChPAD, "Majorcelltype")
rm(Ex_ChPAD)
gc()


###Venns
####Interaction

interact_NCI <- unique(paste(Ex_ChPNCI_df.net$interaction_name))
interact_res <- unique(paste(Ex_ChPres_df.net$interaction_name))
interact_AD <- unique(paste(Ex_ChPAD_df.net$interaction_name))
plts <- ggVennDiagram(list(interact_NCI, interact_res, interact_AD), category.names = c("NCI", "res", "AD")) + theme(legend.position = "none",  text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_Interaction_Comparison.png", sep=""), width=4, height=2.5)

NCI_res_AD <- Reduce(intersect, list(interact_NCI, interact_res, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_res_AD)]
A_res <- interact_res[!(interact_res %in% NCI_res_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_res_AD)]
Intersects <- list(NCI_res_AD, A_NCI, A_res, A_AD)
names(Intersects) <- c("NCI_res_AD", "NCI", "res", "AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_Interaction_NCI_AD.txt", sep=""))


####Pathway

interact_NCI <- unique(paste(Ex_ChPNCI_df.net$pathway_name))
interact_res <- unique(paste(Ex_ChPres_df.net$pathway_name))
interact_AD <- unique(paste(Ex_ChPAD_df.net$pathway_name))
plts <- ggVennDiagram(list(interact_NCI, interact_res, interact_AD), category.names = c("NCI", "res", "AD")) + theme(legend.position = "none", text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_pathway_Comparison.png", sep=""), width=4, height=2.5)

NCI_res_AD <- Reduce(intersect, list(interact_NCI, interact_res, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_res_AD)]
A_res <- interact_res[!(interact_res %in% NCI_res_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_res_AD)]
Intersects <- list(NCI_res_AD, A_NCI, A_res, A_AD)
names(Intersects) <- c("NCI_res_AD", "NCI", "res", "AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_pathway_NCI_AD.txt", sep=""))


####ligand-receptor

interact_NCI <- unique(c(paste(Ex_ChPNCI_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
interact_res <- unique(c(paste(Ex_ChPres_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
interact_AD <- unique(c(paste(Ex_ChPAD_df.net$ligand), paste(Ex_ChPNCI_df.net$receptor)))
plts <- ggVennDiagram(list(interact_NCI, interact_res, interact_AD), category.names = c("NCI", "res", "AD")) + theme(legend.position = "none", text=element_text(size=7))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellChat/Comparison/Venn_ligand-receptor_Comparison.png", sep=""), width=4, height=2.5)

NCI_res_AD <- Reduce(intersect, list(interact_NCI, interact_res, interact_AD))
A_NCI <- interact_NCI[!(interact_NCI %in% NCI_res_AD)]
A_res <- interact_res[!(interact_res %in% NCI_res_AD)]
A_AD <- interact_AD[!(interact_AD %in% NCI_res_AD)]
Intersects <- list(NCI_res_AD, A_NCI, A_res, A_AD)
names(Intersects) <- c("NCI_res_AD", "A_NCI", "A_res", "A_AD")
capture.output(Intersects, file = paste(curr_dir, "/CellChat/Comparison/Intersect_ligand-receptor_NCI_AD.txt", sep=""))



##Merge, compare and DE
###Merge and compare per condition

#merge results
Ex_ChP_merge_CC <- mergeCellChat(list(NCI=Ex_ChPNCI_CC, res=Ex_ChPres_CC, AD=Ex_ChPAD_CC), add.names=names(list(NCI=Ex_ChPNCI_CC, res=Ex_ChPres_CC, AD=Ex_ChPAD_CC)))

# perform differential expression analysis
Ex_ChP_merge_CC <- identifyOverExpressedGenes(Ex_ChP_merge_CC, group.dataset = "datasets", pos.dataset = "AD", features.name = "AD", only.pos = FALSE) #!!!!pos.dataset set AD as default dataset for positive fold change #Keep all results
net <- netMappingDEG(Ex_ChP_merge_CC, features.name = "AD") #select DE

#select up and down regulated pathways with a logFC of .1
net.up <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "AD", ligand.logFC = 0.1, ligand.pvalues=0.01, receptor.logFC = 0.1, receptor.pvalues=0.01) #!!!!up in "ad"
net.down <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "AD", ligand.logFC = -0.1, ligand.pvalues=0.01, receptor.logFC = -0.1, receptor.pvalues=0.01) #!!!!up in "NCI" = down in AD
net_sig <- rbind(net.up, net.down)

# perform differential expression analysis
Ex_ChP_merge_CC <- identifyOverExpressedGenes(Ex_ChP_merge_CC, group.dataset = "datasets", pos.dataset = "res", features.name = "res", only.pos = FALSE) #!!!!pos.dataset set AD as default dataset for positive fold change 
net <- netMappingDEG(Ex_ChP_merge_CC, features.name = "res") #select DE

#select up and down regulated pathways with a logFC of .1
net.up <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "res", ligand.logFC = 0.1, ligand.pvalues=0.01, receptor.logFC = 0.1, receptor.pvalues=0.01) #!!!!up in "ad"
gene.up <- extractGeneSubsetFromPair(net.up, Ex_ChP_merge_CC) #get genes
net.down <- subsetCommunication(Ex_ChP_merge_CC, net = net, datasets = "res", ligand.logFC = -0.1, ligand.pvalues=0.01, receptor.logFC = -0.1, receptor.pvalues=0.01) #!!!!up in "NCI" = down in AD
gene.down <- extractGeneSubsetFromPair(net.down, Ex_ChP_merge_CC)
net_sig <- rbind(net_sig, net.up, net.down)

saveRDS(Ex_ChP_merge_CC, paste(curr_dir, "/CellChat/Ex_ChP_merge_CC_ADNCI.RDS", sep=""))



##DEG with limma instead
###function and prep CC

#Construct ADvres results DEG table
DEG_source_fun <- function(limma_res, net){
  #Merges CC net results to limma results table
  #limma_res: limma result
  #net: output from subsetCommunications
  
  #Add source limma results
  net_sig <- data.frame(t(data.frame(sapply(
    unique(net$source), function(j){
      net_S <- net[net$source==j,]
      limma_S <- limma_res[limma_res$Celltype==j,]
      colnames(limma_S) <- paste0(colnames(limma_S), "_Ligand")
      sapply(1:nrow(net_S), function(i) bind_cols(net_S[i,], limma_S[limma_S$gene==net_S[i,]$ligand,]))
    }
  ))))
  rownames(net_sig) <- NULL
  net_sig <- unnest(net_sig, cols = colnames(net_sig))
  
  #Add receptorA limma results
  net_sig <- data.frame(t(data.frame(sapply(
    unique(net_sig$target), function(j){
      net_S <- net_sig[net_sig$target==j,]
      limma_S <- limma_res[limma_res$Celltype==j,]
      colnames(limma_S) <- paste0(colnames(limma_S), "_ReceptorA")
      sapply(1:nrow(net_S), function(i) bind_cols(net_S[i,], limma_S[limma_S$gene==net_S[i,]$receptorA,]))
    }
  ))))
  rownames(net_sig) <- NULL
  net_sig <- unnest(net_sig, cols = colnames(net_sig))
  
  #Add receptorB limma results
  net_sigB <- data.frame(t(data.frame(sapply(
    unique(net_sig$target), function(j){
      net_S <- net_sig[net_sig$target==j,]
      limma_S <- limma_res[limma_res$Celltype==j,]
      colnames(limma_S) <- paste0(colnames(limma_S), "_ReceptorB")
      sapply(1:nrow(net_S), function(i) bind_cols(net_S[i,], limma_S[limma_S$gene==net_S[i,]$receptorB,]))
    }
  ))))
  rownames(net_sigB) <- NULL
  net_sigB <- unnest(net_sigB, cols = colnames(net_sigB))
  net_sigB  <- net_sigB[complete.cases(net_sigB), ]
  
  net_sig  <- net_sig[!complete.cases(net_sig), ]
  net_sig <- bind_rows(net_sig, net_sigB)
  return(net_sig)
}


#get net
net <- subsetCommunication(Ex_ChP_merge_CC)


net <- rbind(net$NCI, net$res, net$AD)
net <- net[,c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name", "pathway_name", "annotation", "evidence")]

#get genes
net <- bind_cols(net, data.frame(t(sapply(str_split(net$receptor,"_",),'[',1:2))))
colnames(net) <- c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name", "pathway_name", "annotation", "evidence", "receptorA", "receptorB")


genes <- c(net$ligand, net$receptorA, net$receptorB)
genes <- unique(genes)
genes <- genes[!is.na(genes)]

cutoff_pval <- 0.05
cutoff_logFC <- 0.2


###NCIvAD

#make limma_res
limma_cogdxadnc_ADdvNCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdxadnc_ADdvNCI_Epithelial)
limma_cogdxadnc_ADdvNCI_Fibroblast <- cbind(Celltype="Fibroblast", limma_cogdxadnc_ADdvNCI_Fibroblast)
limma_cogdxadnc_ADdvNCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdxadnc_ADdvNCI_Endothelial)
limma_cogdxadnc_ADdvNCI_Mural <- cbind(Celltype="Mural", limma_cogdxadnc_ADdvNCI_Mural)
limma_cogdxadnc_ADdvNCI_Immune <- cbind(Celltype="Immune", limma_cogdxadnc_ADdvNCI_Immune)
limma_cogdxadnc_ADdvNCI_Oligo <- cbind(Celltype="Oligo", limma_cogdxadnc_ADdvNCI_Oligo)
limma_cogdxadnc_ADdvNCI_Neuro <- cbind(Celltype="Neuro", limma_cogdxadnc_ADdvNCI_Neuro)

limma_res <- bind_rows(limma_cogdxadnc_ADdvNCI_Epithelial, limma_cogdxadnc_ADdvNCI_Fibroblast, limma_cogdxadnc_ADdvNCI_Endothelial, limma_cogdxadnc_ADdvNCI_Mural, limma_cogdxadnc_ADdvNCI_Immune, limma_cogdxadnc_ADdvNCI_Oligo, limma_cogdxadnc_ADdvNCI_Neuro)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdx_ADvNCI <- DEG_source_fun(limma_res, net)
net_cogdxadnc_ADdvNCI_sig <- net_cogdx_ADvNCI[net_cogdx_ADvNCI$P.Value_Ligand<cutoff_pval & abs(net_cogdx_ADvNCI$logFC_Ligand)>cutoff_logFC,]  #& c(net_cogdxadnc_ADdvresilient$P.Value_ReceptorA<cutoff_pval | net_cogdxadnc_ADdvresilient$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdxadnc_ADdvresilient$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdxadnc_ADdvresilient$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdxadnc_ADdvNCI_sig <- net_cogdxadnc_ADdvNCI_sig[complete.cases(net_cogdxadnc_ADdvNCI_sig$source),]
view(net_cogdxadnc_ADdvNCI_sig)



###NCIvres

#make limma_res
limma_cogdxadnc_resilientvNCI_Epithelial <- cbind(Celltype="Epithelial", limma_cogdxadnc_resilientvNCI_Epithelial)
limma_cogdxadnc_resilientvNCI_Fibroblast <- cbind(Celltype="Fibroblast", limma_cogdxadnc_resilientvNCI_Fibroblast)
limma_cogdxadnc_resilientvNCI_Endothelial <- cbind(Celltype="Endothelial", limma_cogdxadnc_resilientvNCI_Endothelial)
limma_cogdxadnc_resilientvNCI_Mural <- cbind(Celltype="Mural", limma_cogdxadnc_resilientvNCI_Mural)
limma_cogdxadnc_resilientvNCI_Immune <- cbind(Celltype="Immune", limma_cogdxadnc_resilientvNCI_Immune)
limma_cogdxadnc_resilientvNCI_Oligo <- cbind(Celltype="Oligodendrocyte", limma_cogdxadnc_resilientvNCI_Oligo)
limma_cogdxadnc_resilientvNCI_Neuro <- cbind(Celltype="Neuron", limma_cogdxadnc_resilientvNCI_Neuro)

limma_res <- bind_rows(limma_cogdxadnc_resilientvNCI_Epithelial, limma_cogdxadnc_resilientvNCI_Fibroblast, limma_cogdxadnc_resilientvNCI_Endothelial, limma_cogdxadnc_resilientvNCI_Mural, limma_cogdxadnc_resilientvNCI_Immune, limma_cogdxadnc_resilientvNCI_Oligo, limma_cogdxadnc_resilientvNCI_Neuro)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdxadnc_resilientvNCI <- DEG_source_fun(limma_res, net)
net_cogdxadnc_resilientvNCI_sig <- net_cogdxadnc_resilientvNCI[net_cogdxadnc_resilientvNCI$P.Value_Ligand<cutoff_pval & abs(net_cogdxadnc_resilientvNCI$logFC_Ligand)>cutoff_logFC,] #& c(net_cogdxadnc_ADdvresilient$P.Value_ReceptorA<cutoff_pval | net_cogdxadnc_ADdvresilient$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdxadnc_ADdvresilient$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdxadnc_ADdvresilient$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdxadnc_resilientvNCI_sig <- net_cogdxadnc_resilientvNCI_sig[complete.cases(net_cogdxadnc_resilientvNCI_sig$source),]
view(net_cogdxadnc_resilientvNCI_sig)



###resvAD

#make limma_res
limma_cogdxadnc_ADdvresilient_Epithelial <- cbind(Celltype="Epithelial", limma_cogdxadnc_ADdvresilient_Epithelial)
limma_cogdxadnc_ADdvresilient_Fibroblast <- cbind(Celltype="Fibroblast", limma_cogdxadnc_ADdvresilient_Fibroblast)
limma_cogdxadnc_ADdvresilient_Endothelial <- cbind(Celltype="Endothelial", limma_cogdxadnc_ADdvresilient_Endothelial)
limma_cogdxadnc_ADdvresilient_Mural <- cbind(Celltype="Mural", limma_cogdxadnc_ADdvresilient_Mural)
limma_cogdxadnc_ADdvresilient_Immune <- cbind(Celltype="Immune", limma_cogdxadnc_ADdvresilient_Immune)
limma_cogdxadnc_ADdvresilient_Oligo <- cbind(Celltype="Oligodendrocyte", limma_cogdxadnc_ADdvresilient_Oligo)
limma_cogdxadnc_ADdvresilient_Neuro <- cbind(Celltype="Neuron", limma_cogdxadnc_ADdvresilient_Neuro)

limma_res <- bind_rows(limma_cogdxadnc_ADdvresilient_Epithelial, limma_cogdxadnc_ADdvresilient_Fibroblast, limma_cogdxadnc_ADdvresilient_Endothelial, limma_cogdxadnc_ADdvresilient_Mural, limma_cogdxadnc_ADdvresilient_Immune, limma_cogdxadnc_ADdvresilient_Oligo, limma_cogdxadnc_ADdvresilient_Neuro)
limma_res <- limma_res[limma_res$gene %in% genes, 1:8] #extract result

net_cogdxadnc_ADdvresilient <- DEG_source_fun(limma_res, net)
net_cogdxadnc_ADdvresilient_sig <- net_cogdxadnc_ADdvresilient[net_cogdxadnc_ADdvresilient$P.Value_Ligand<cutoff_pval & abs(net_cogdxadnc_ADdvresilient$logFC_Ligand)>cutoff_logFC ,] #& c(net_cogdxadnc_ADdvresilient$P.Value_ReceptorA<cutoff_pval | net_cogdxadnc_ADdvresilient$P.Value_ReceptorB<cutoff_pval) & c(abs(net_cogdxadnc_ADdvresilient$logFC_ReceptorA)>cutoff_logFC | abs(net_cogdxadnc_ADdvresilient$logFC_ReceptorB)>cutoff_logFC) #Don't exclude cases where the receptor is not changing or is changing in the opposite direction. Receptors can be up or down regulated (or nothing) due to changes in the signal.
net_cogdxadnc_ADdvresilient_sig <- net_cogdxadnc_ADdvresilient_sig[complete.cases(net_cogdxadnc_ADdvresilient_sig$source),]
view(net_cogdxadnc_ADdvresilient_sig)



###Aggregate net_sigs

net <- rbind(net_cogdx_ADvNCI, net_cogdxadnc_ADdvresilient, net_cogdxadnc_resilientvNCI)
write.csv(net, paste0(curr_dir, "/CellChat/Comparison/net_pathways.csv"))
net_sig <- rbind(net_cogdxadnc_ADdvNCI_sig, net_cogdxadnc_ADdvresilient_sig, net_cogdxadnc_resilientvNCI_sig)
write.csv(net_sig, paste0(curr_dir, "/CellChat/Comparison/net_sig_pathways.csv"))


###Aggregate by pathway

ligandlogFC <- aggregate(net_cogdx_ADvNCI$logFC_Ligand, by=list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdx_ADvNCI$logFC_ReceptorA, list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdx_ADvNCI$logFC_ReceptorB, list(net_cogdx_ADvNCI$pathway_name, net_cogdx_ADvNCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_ADvNCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_ADvNCI <- cbind(combinedlogFC_ADvNCI[,1:2], rowMeans(combinedlogFC_ADvNCI[3:5], na.rm = TRUE))
view(combinedlogFC_ADvNCI)

ligandlogFC <- aggregate(net_cogdxadnc_ADdvresilient$logFC_Ligand, by=list(net_cogdxadnc_ADdvresilient$pathway_name, net_cogdxadnc_ADdvresilient$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdxadnc_ADdvresilient$logFC_ReceptorA, list(net_cogdxadnc_ADdvresilient$pathway_name, net_cogdxadnc_ADdvresilient$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdxadnc_ADdvresilient$logFC_ReceptorB, list(net_cogdxadnc_ADdvresilient$pathway_name, net_cogdxadnc_ADdvresilient$Celltype_ReceptorB), mean, all=T)
combinedlogFC_ADvres <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_ADvres <- cbind(combinedlogFC_ADvres[,1:2], rowMeans(combinedlogFC_ADvres[3:5], na.rm = TRUE))
view(combinedlogFC_ADvres)

ligandlogFC <- aggregate(net_cogdxadnc_resilientvNCI$logFC_Ligand, by=list(net_cogdxadnc_resilientvNCI$pathway_name, net_cogdxadnc_resilientvNCI$Celltype_Ligand), FUN=mean, all=T)
receptorAlogFC <- aggregate(net_cogdxadnc_resilientvNCI$logFC_ReceptorA, list(net_cogdxadnc_resilientvNCI$pathway_name, net_cogdxadnc_resilientvNCI$Celltype_ReceptorA), mean, all=T)
receptorBlogFC <- aggregate(net_cogdxadnc_resilientvNCI$logFC_ReceptorB, list(net_cogdxadnc_resilientvNCI$pathway_name, net_cogdxadnc_resilientvNCI$Celltype_ReceptorB), mean, all=T)
combinedlogFC_resvNCI <- merge(receptorAlogFC, merge(receptorBlogFC, ligandlogFC, by=c("Group.1", "Group.2"), all=T), by=c("Group.1", "Group.2"), all=T)
combinedlogFC_resvNCI <- cbind(combinedlogFC_resvNCI[,1:2], rowMeans(combinedlogFC_resvNCI[3:5], na.rm = TRUE))
view(combinedlogFC_resvNCI)


##Plot differences
###heatmap comps

#barplots
plts <- compareInteractions(Ex_ChP_merge_CC, show.legend = F, group=c(1,2,3), color.use=c("skyblue", "#F86247", "red3")) + theme(text=element_text(size=7)) #Barplot of counts
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_bar.pdf", sep=""), width=2, height=2.5)

plts <- compareInteractions(Ex_ChP_merge_CC, show.legend = F, group=c(1,2,3), measure="weight", color.use=c("skyblue", "#F86247", "red3")) + theme(text=element_text(size=7)) #Barplot of weights
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_bar.pdf", sep=""), width=2, height=2.5)

#Network
plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,2)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_AD-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,2), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_AD-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(2,3)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_AD-res.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(2,3), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_AD-res.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,3)) #networkplot of counts
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Count_network_res-NCI.pdf", sep=""))

plts <- netVisual_diffInteraction(Ex_ChP_merge_CC, weight.scale = T, comparison=c(1,3), measure = "weight") #networkplot of weights
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Weight_network_res-NCI.pdf", sep=""))

signaling_NCI <- dimnames(Ex_ChP_merge_CC@net$NCI$prob)[[3]]
signaling_res <- dimnames(Ex_ChP_merge_CC@net$res$prob)[[3]]
signaling_AD <- dimnames(Ex_ChP_merge_CC@net$AD$prob)[[3]]
signaling_shared <- intersect(signaling_AD, intersect(signaling_res, signaling_NCI))
signaling_sig <- unique(net_cogdxadnc_ADdvNCI_sig$interaction_name) # 284 sig pathways generated using CC_limma.R script
signaling_final <- intersect(signaling_sig, signaling_shared)

Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI","AD"), signaling=signaling_final, n_sample=20, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"), group.color=c("skyblue", "red3"), col.max=0.2)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI","AD"), signaling=signaling_final, n_sample=20, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"), group.color=c("skyblue", "red3"), col.max=8)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("NCI","res"), signaling=signaling_final, n_sample=15, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"), group.color=c("skyblue", "#F86247"), col.max=0.2)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("NCI","res"), signaling=signaling_final, n_sample=15, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"), group.color=c("skyblue", "#F86247"), col.max=8)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="weight", comparison=c("res","AD"), signaling=signaling_final, n_sample=17, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"), group.color=c("#F86247", "red3"), col.max=0.2)
Compare_interaction_heatbar_fun(Ex_ChP_merge_CC, measure="count", comparison=c("res","AD"), signaling=signaling_final, n_sample=17, cell.color=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"), group.color=c("#F86247", "red3"), col.max=8)


###stackedbar sig paths

net_sig <- rbind(net_cogdxadnc_ADdvNCI_sig, net_cogdxadnc_resilientvNCI_sig, net_cogdxadnc_ADdvresilient_sig)
net_sig <- net_sig[order(net_sig$source, decreasing=T),]
view(net_sig)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_count.pdf", sep=""), width=3, height=6)
remove(plts)
Count_ADNC_sig <- rankNet(Ex_ChP_merge_CC, mode = "comparison", measure = "count", do.stat = TRUE, return.data=T)$signaling.contribution

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_weight.pdf", sep=""), width=3, height=6)
remove(plts)
Weight_ADNC_sig <- rankNet(Ex_ChP_merge_CC, mode = "comparison", measure = "weight", do.stat = TRUE, return.data=T)$signaling.contribution

#Heatmap of differences
CC_dat <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "count", signaling = unique(net$pathway_name), do.stat = TRUE, return.data=T)$signaling.contribution
CC_dat$contribution.scaled <- NULL
CC_dat$contribution.relative.1 <- NULL
CC_dat$contribution.relative.2 <- NULL
CC_dat_val <- as.data.frame(pivot_wider(CC_dat, names_from = "group",  values_from = "contribution"))
rownames(CC_dat_val) <- CC_dat_val$name
CC_dat_val$name <- NULL
mat <- t(scale(t(CC_dat_val)))
mat <- mat[order(match(rownames(mat), rev(unique(net$pathway_name)))),]

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="Relative Number of Interactions", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = colorRamp2(c(ceiling(min(mat, na.rm=T)), floor(max(mat, na.rm=T))), c("white", "black")), na_col = "white",)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_count_heatmap.pdf", font="Arial", sep = ""), width=w, height=h)

CC_dat <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2,3), measure = "weight", signaling = unique(net$pathway_name), do.stat = TRUE, return.data=T)$signaling.contribution
CC_dat$contribution.scaled <- NULL
CC_dat$contribution.relative.1 <- NULL
CC_dat$contribution.relative.2 <- NULL
CC_dat_val <- as.data.frame(pivot_wider(CC_dat, names_from = "group",  values_from = "contribution"))
rownames(CC_dat_val) <- CC_dat_val$name
CC_dat_val$name <- NULL
mat <- t(scale(t(CC_dat_val)))
mat <- mat[order(match(rownames(mat), rev(unique(net$pathway_name)))),]

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="Relative Number of Interactions", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = colorRamp2(c(ceiling(min(mat, na.rm=T)), floor(max(mat, na.rm=T))), c("white", "black")), na_col = "white",)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_weight_heatmap.pdf", font="Arial", sep = ""), width=w, height=h)

#ADvNCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvNCI_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvNCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)

#ADvres
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(2,3), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvres_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(2,3), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("#F86247", "red3"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_ADvres_weight.pdf", sep=""), width=3, height=6)
remove(plts)

#resvNCI
plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2), measure = "count", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_resvNCI_count.pdf", sep=""), width=3, height=6)
remove(plts)

plts <- rankNet(Ex_ChP_merge_CC, mode = "comparison", comparison=c(1,2), measure = "weight", signaling = unique(net$pathway_name), bar.w = 0.65, font.size = 7, color.use=c("skyblue", "#F86247"), stacked = T, do.stat = TRUE)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/stackedbar_sig_resvNCI_weight.pdf", sep=""), width=3, height=6)
remove(plts)


###Heatmap of sig paths

#preprocess
mat1 <- prep_mat(object=Ex_ChPNCI_CC, slot.name="netP", thresh=0.05, pattern="all")
mat2 <- prep_mat(object=Ex_ChPres_CC, slot.name="netP", thresh=0.05, pattern="all")
mat3 <- prep_mat(object=Ex_ChPAD_CC, slot.name="netP", thresh=0.05, pattern="all")

mat <- rbind.fill(mat1, mat2, mat3)
mat <- mat[,order(match(colnames(mat), rev(unique(net$pathway_name))))]

netAnalysis_signalingRole_heatmap_combofun(mat, pattern = "all", filename = "NCI-res-AD", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna", "#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna", "#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))


###Communication bubble plot of sig paths

pairLR.use <- rbind(net_sig[, "interaction_name", drop = F])
plts <- netVisual_bubble(Ex_ChP_merge_CC, pairLR.use = pairLR.use, sources.use = net_sig$source, comparison = c(1, 2, 3),  angle.x = 90, remove.isolate = T, color.heatmap = c("white", "black", title.name = paste("Differential signaling"), font.size = 7, font.size.title = 7, grid.on = F))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellChat/Comparison/Dot_LRpair.pdf", sep=""), width=3+length(unique(net_sig$source))*30*.2, height=3+length(unique(pairLR.use$interaction_name))*.2)


###Chord plot

net_path <- unique(c(net$pathway_name))
chords1 <- function(i){
  #pdf(paste(curr_dir, "/CellChat/Comparison/chord_NCI_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_chord_cell(Ex_ChPNCI_CC, signaling = i, scale=F, remove.isolate=T, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
  #dev.off()
  #netVisual_chord_cell(Ex_ChPNCI_CC, signaling = i, scale=F, remove.isolate=T, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B"))
}
chords2 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_res_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_chord_cell(Ex_ChPres_CC, signaling = i, scale=F, remove.isolate=T, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
  dev.off()
  netVisual_chord_cell(Ex_ChPres_CC, signaling = i, scale=F, remove.isolate=T, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
}
chords3 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_AD_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_chord_cell(Ex_ChPAD_CC, signaling = i, scale=F, remove.isolate=T, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
  dev.off()
  netVisual_chord_cell(Ex_ChPAD_CC, signaling = i, scale=F, remove.isolate=T, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
}


for (i in net_path) {
  skip<-FALSE
  tryCatch(chords1(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords2(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords3(i), error=function(e){skip<-TRUE})
  if(skip){next}
}
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()



pdf(paste(curr_dir, "/CellChat/Comparison/chord_NCI_allcells.pdf", sep=""), height=4, width=4)
netVisual_chord_gene(Ex_ChPNCI_CC, slot.name="netP", signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_chord_gene(Ex_ChPNCI_CC, slot.name="netP", signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))

pdf(paste(curr_dir, "/CellChat/Comparison/chord_res_allcells.pdf", sep=""), height=4, width=4)
netVisual_chord_gene(Ex_ChPres_CC, slot.name="netP", signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_chord_gene(Ex_ChPres_CC, slot.name="netP", signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))

pdf(paste(curr_dir, "/CellChat/Comparison/chord_AD_allcells.pdf", sep=""), height=4, width=4)
netVisual_chord_gene(Ex_ChPAD_CC, slot.name="netP", signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_chord_gene(Ex_ChPAD_CC, slot.name="netP", signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))


net_path <- unique(c(net$pathway_name))
chords1 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_NCI_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_chord_gene(Ex_ChPNCI_CC, slot.name="netP", sources.use=i, signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
  dev.off()
  netVisual_chord_gene(Ex_ChPNCI_CC, slot.name="netP", sources.use=i, signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
}
chords2 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_res_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_chord_gene(Ex_ChPres_CC, slot.name="netP", sources.use=i, signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
  dev.off()
  netVisual_chord_gene(Ex_ChPres_CC, slot.name="netP", sources.use=i, signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
}
chords3 <- function(i){
  pdf(paste(curr_dir, "/CellChat/Comparison/chord_AD_", i, ".pdf", sep=""), height=4, width=4)
  netVisual_chord_gene(Ex_ChPAD_CC, slot.name="netP", sources.use=i, signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
  dev.off()
  netVisual_chord_gene(Ex_ChPAD_CC, slot.name="netP", sources.use=i, signaling = net_path, color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
}


for (i in c(1,2,3,4,5)) {
  skip<-FALSE
  tryCatch(chords1(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords2(i), error=function(e){skip<-TRUE})
  if(skip){next}
  tryCatch(chords3(i), error=function(e){skip<-TRUE})
  if(skip){next}
}
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()


###Circle plots

#NCI
Ex_ChPNCI_CC <- aggregateNet(Ex_ChPNCI_CC,  signaling=unique(net_sig$pathway_name), remove.isolate=F)
pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Weight_NCI.pdf"))
netVisual_circle(Ex_ChPNCI_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_circle(Ex_ChPNCI_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))

pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Count_NCI.pdf"))
netVisual_circle(Ex_ChPNCI_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_circle(Ex_ChPNCI_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))

#res
Ex_ChPres_CC <- aggregateNet(Ex_ChPNCI_CC,  signaling=unique(net_sig$pathway_name), remove.isolate=F)
pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Weight_res.pdf"))
netVisual_circle(Ex_ChPres_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_circle(Ex_ChPres_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))

pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Count_res.pdf"))
netVisual_circle(Ex_ChPres_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_circle(Ex_ChPres_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))

#AD
Ex_ChPNCI_CC <- aggregateNet(Ex_ChPAD_CC,  signaling=unique(net_sig$pathway_name), remove.isolate=F)
pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Weight_AD.pdf"))
netVisual_circle(Ex_ChPAD_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_circle(Ex_ChPAD_CC@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))

pdf(paste0(curr_dir, "/CellChat/Comparison/Circle_Count_AD.pdf"))
netVisual_circle(Ex_ChPAD_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))
dev.off()
netVisual_circle(Ex_ChPAD_CC@net$count, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "lightsalmon", "sienna"))




