###
### 051 differential_plots.R
###
# Purpose: Produce spatial molecular expression visualizations
# Dependencies:
source("/050 ST_specific_plots_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder

##Resilience
###Epithelial
####Metabolism
major_cols <- c("Epithelial"="#C48F45", "Fibroblast"="white", "Endothelial"="white", "Mural"="white", "Immune"="white", "Parenchyma"="white", "low-quality"="white")
DefaultAssay(Cos_ChP) <- "Nanostring"

genelist <- c("NQO1", "ODC1", "CTSD", "LAMP2", "CLU", "PINK1", "SQSTM1", "P4HB", "PFKL")
valid_FOVs <- get_valid_FOVs(Cos_ChP, major_cell_type = "Epithelial", genelist, min_cells=400, max_cells=600, min_total_cells=400)
interest_FOVs <- get_valid_FOVs_for_category(Cos_ChP, genelist, valid_FOVs, decreasing=T)

# Sort the dataframe by Difference and select the top 3 unique FOVs for AD and NCI (adjust the '1:3' if you want to plot more)
top_ad_fovs <- unique(interest_FOVs$ad_fov_names[1:4])
top_nci_fovs <- unique(interest_FOVs$nci_fov_names[1:4])

# Combine the top AD and NCI FOVs
all_fovs <- c(top_ad_fovs, top_nci_fovs, "Slide.1.89")

# Plot first 3 unique AD FOVs and first 3 unique NCI FOVs and save to PDF
pdf(paste0(curr_dir, "/ST_images/Epi_Cilia.pdf"))  # Open a PDF device to save the plots
for (fov in all_fovs) {
  # Ensure segmentation:
  DefaultBoundary(Cos_ChP@images[[fov]]) <- 'segmentation'
  cat("Plotting FOV:", fov, "\n")
  plot <- ImageDimPlot(
    Cos_ChP, 
    fov = fov, 
    group.by = "Majorcelltype", 
    cols = major_cols, 
    alpha = 0.8, 
    molecules = genelist, 
    nmols = 10000, 
    #mols.cols = c("#E41A1C", "#F781BF", "green3", "lightgreen"), 
    dark.background = FALSE, 
    border.color = "black", 
    coord.fixed = FALSE, 
    border.size = 0.2, 
    mols.size = 0.3
  ) + ggtitle(paste("Cilia; FOV:", fov))  # Add the FOV name as title
  print(plot)  # Explicitly print the plot
}
dev.off()  # Close the PDF device


scaled_imagedimplot(Cos_ChP, gene="TUBB4B", fov1="Slide.2.23", fov2="Slide.2.124", major_cell_type="Epithelial")
scaled_imagedimplot(Cos_ChP, "FOXJ1", "Slide.2.23", "Slide.2.124", "Epithelial")
