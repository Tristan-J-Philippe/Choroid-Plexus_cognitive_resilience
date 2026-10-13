###
### 030 limma_pattern_fgsea_util.R 
###
# Purpose: limma-Voom analyses, fgsea GO term analyses, pattern matching analysis functions and GSCA module analysis functions
# Dependencies:
library(lattice)
library(limma)
library(edgeR)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(org.Hs.eg.db)
library(fgsea)
library(msigdbr)
library(ggrepel)
library(splines)
library(stringr)
library(tidyverse)
library(ggVennDiagram)
library(dplyr)
library(clusterProfiler)
set.seed(1234)

#Set and load dependent variables for limma/DESeq functions
Gene.search <- read.csv("/DEG to GO.csv")
GO_labelled <- read.csv("/GO labelled.csv")

FDRcutoff <- 0.01
LFcutoff <- 0.25
pct_cutoff <- 0.05

db_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
db_CC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
db_MC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
db_list <- rbind(db_BP, db_CC, db_MC)

pctincell <- function(object){
  #Function to calculate the percentage of cells that express all genes in object.
  result <- data.frame(object@assays[["RNA"]]@counts)
  result <- data.frame(rowSums(result>0)/ncol(result))
  colnames(result) <- "Percent_in_cells"
  return(result)
}

aggcols_fun <- function(cts, min_to_sum=100000){
  #function to aggregate Lanes (0 and 1) from projids in dataset when the number of counts from BOTH Lanes is too low to yield meaningful results.
  #min_to_sum: Minimum number of reads that need to be reached for columns to not be summed together. 
  low_sum <- colSums(cts)<min_to_sum #Which cols fall below or above min
  cts_s <- cts[,low_sum==FALSE] #save cols above min
  cts_l <- cts[,low_sum==TRUE] #cols below min
  if(is.null(ncol(cts_l))==FALSE){
    
    name_cols_projid <- sapply(str_split(colnames(cts_l),"_ChP",),'[',1) #get projid names
    name_cols_Lane <- sapply(str_split(colnames(cts_l),"_ChP",),'[', 2) #get Lane names
    
    #set up while/for loop
    cts_b <- matrix(nrow=nrow(cts_l))
    cts_b <- cts_b[,-1]
    i=1
    
    #for every col in cts_l see if they can be combined else keep it
    while(i<length(name_cols_projid)){
      if(name_cols_projid[i]==name_cols_projid[i+1]){
        name_col1 <- colnames(cts_l)[i]
        name_col2 <- name_cols_Lane[i+1]
        name_cols <- colnames(cts_b)
        cts_b <- cbind(cts_b, rowSums(cbind(cts_l[,i], cts_l[,i+1])))#add them to the df
        colnames(cts_b) <- c(name_cols, paste(name_col1, name_col2, sep="-"))
        i <- i+2
      }else{
        name_col1 <- colnames(cts_l)[i]
        name_cols <- colnames(cts_b)
        cts_b <- cbind(cts_b, cts_l[,i])#add them to the df
        colnames(cts_b) <- c(name_cols, paste(name_col1))
        i <- i+1
      }
    }
    
    #if land on the last col make sure to keep it.
    if(i==length(name_cols_projid)){
      name_cols <- colnames(cts_b)
      name_col1 <- colnames(cts_l)[i]
      cts_b <- cbind(cts_b, cts_l[,i])
      colnames(cts_b) <- c(name_cols, name_col1)
      i <- i+1
    }
    
    cts <- cbind(cts_b, cts_s)
    cts <- cts[,order(colnames(cts))]
  }
  return(cts)
}

aggcols_fun_b <- function(cts, min_to_sum=100000, search_direction){
  #function to aggregate Lane (0 and 1) from projids in dataset when the number of counts in ANY Lane below threshold.
  #min_to_sum: Minimum number of reads that need to be reached for columns to not be summed together.
  #search_direction: 0 (forwards) or 1 (backwards)
  low_sum <- colSums(cts)<min_to_sum #Which cols fall below or above min
  name_cols_projid <- sapply(str_split(colnames(cts),"_ChP",),'[',1) #get projid names
  name_cols_Lane <- sapply(str_split(colnames(cts),"_ChP",),'[', 2) #get Lane names
  
  #set up while/for loop
  cts_b <- matrix(nrow=nrow(cts))
  cts_b <- cts_b[,-1]
  
  i <- 1
  
  #for every col in cts that is low see if it can be combined else keep it
  while(i<length(name_cols_projid)){
    if(name_cols_projid[i]==name_cols_projid[i+1] & low_sum[i+search_direction]==TRUE){
      name_col1 <- colnames(cts)[i]
      name_col2 <- name_cols_Lane[i+1]
      name_cols <- colnames(cts_b)
      cts_b <- cbind(cts_b, rowSums(cbind(cts[,i], cts[,i+1])))#add them to the df
      colnames(cts_b) <- c(name_cols, paste(name_col1, name_col2, sep="-"))
      i <- i+2
    }else{
      name_col1 <- colnames(cts)[i]
      name_cols <- colnames(cts_b)
      cts_b <- cbind(cts_b, cts[,i])#add them to the df
      colnames(cts_b) <- c(name_cols, paste(name_col1))
      i <- i+1
    }
  }
  
  #if land on the last col make sure to keep it.
  if(i==length(name_cols_projid)){
    name_cols <- colnames(cts_b)
    name_col1 <- colnames(cts)[i]
    cts_b <- cbind(cts_b, cts[,i])
    colnames(cts_b) <- c(name_cols, name_col1)
    i <- i+1
  }
  
  cts <- cts_b
  cts <- cts[,order(colnames(cts))]
  return(cts)
}

datprepfun_lane <- function(srt, celltype){
  #Function to select and aggregate data by cell type and projid*Lane
  sub <- subset(srt, Subcelltype==celltype) #Majorcelltype/Subcelltype/EpiFibBAMcelltype
  gc()
  #Aggregate
  cts <- AggregateExpression(sub, group.by = "Lane_projid", assay = "RNA", slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  cts <- cts[,order(colnames(cts))]
  
  #Aggregate or combine Lanes that have too few counts
  min_to_sum <- 100000 #set minimum number of counts that need to be reached for columns to not be summed together.
  cts <- aggcols_fun(cts, min_to_sum) #run aggcols_fun the number of Lanes you have divided by 2; rounding up!
  cts <- aggcols_fun(cts, min_to_sum)
  cts <- aggcols_fun(cts, min_to_sum)
  
  cts <- aggcols_fun_b(cts, min_to_sum, 0) #Aggregate even if any are <min_to_sum, forwards
  cts <- aggcols_fun_b(cts, min_to_sum , 1) #Aggregate even if any are <min_to_sum, backwards
  
  #Exclude Lanes that have too few counts
  min_to_exclude <- 60000
  low_sum <- colSums(cts)<min_to_exclude #Which cols fall below or above min
  cts <- cts[,low_sum==FALSE] #keep cols above min
  
  assign(paste("cts_", celltype, sep=""),  cts, envir = .GlobalEnv)
  
  #Find percentage of cells expressing genes
  pct <- pctincell(sub)
  assign(paste("pct_", celltype, sep=""),  pct, envir = .GlobalEnv)
  gc()
}

datprepfun <- function(srt, celltype){
  #Function to select and aggregate data by cell type and projid*Lane
  sub <- subset(srt, Subcelltype==celltype) #Majorcelltype/Subcelltype
  #Aggregate
  gc()
  cts <- AggregateExpression(sub, group.by = "projid", assay = "RNA", slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  
  cts <- cts[,order(colnames(cts))]
  
  #Exclude Lanes that have too few counts
  min_to_exclude <- 60000
  low_sum <- colSums(cts)<min_to_exclude #Which cols fall below or above min
  cts <- cts[,low_sum==FALSE] #keep cols above min
  
  assign(paste0("cts_", celltype),  cts, envir = .GlobalEnv)
  
  #Find percentage of cells expressing genes
  pct <- pctincell(sub)
  assign(paste0("pct_", celltype),  pct, envir = .GlobalEnv)
  #write.csv(cts, paste0(curr_dir, "/limma/Pseudobulk_per_projid/cst_", celltype, ".csv"))
  gc()
}

limmafun_lane <- function(cts, pct){
  celltype <- deparse(substitute(cts))
  celltype <- substr(celltype,5,30)
  
  #Select relevant rows from Design Matrix  
  dm <- matrix(ncol=ncol(Maindm))
  dm <- dm[-1,]
  colnames(dm) <- colnames(Maindm)
  
  for(i in colnames(cts)){
    base_name <- sapply(str_split(i,"_",),'[',1) #get projid names
    if(base_name %in% rownames(Maindm)){
      name_rows <- rownames(dm)
      row_add <- Maindm[base_name,]
      dm <- rbind(dm, row_add)
      rownames(dm) <- c(name_rows, paste(i))
    }else{print(paste("Warning colname of data:", i, "not in rowname of Design Matrix"))}
  }
  
  #Enforce Maindm and cts name match
  if(all(rownames(dm)==colnames(cts))==FALSE){
    dm <- dm[rownames(dm) %in% colnames(cts),]      
    cts <- cts[,colnames(cts) %in% rownames(dm)]
  }
  
  #order dm and cts based on first variable of interest
  dm <- dm[order(dm[,1]),]
  cts <- cts[,order(match(colnames(cts), rownames(dm)))]
  
  #Gene expression cutoffs
  cts <- cts[rowSums(cts)>2*ncol(cts),]
  
  #Gene expression cutoff for % in cells
  pct <- pct %>% rownames_to_column("gene")
  pct <- pct[pct$Percent_in_cells>pct_cutoff,]
  cts <- cts[rownames(cts) %in% pct$gene,]
  
  #voom
  vdat <- voom(cts, design = dm, plot=T)
  heatmap(cor(vdata), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")
  
  #technical replicate
  projid <- factor(sapply(str_split(colnames(cts),"_",),'[', 1)) #get projid names
  corfit <- duplicateCorrelation(vdat, block=projid)
  vdat <- voom(cts, design = dm, plot=T, block=projid, correlation=corfit$consensus.correlation)
  vdata <- data.frame(vdat)
  
  #limma
  Mainf <- lmFit(vdat, dm, weights = vdat$weights, block=projid, correlation=corfit$consensus.correlation)
  Mainf <- contrasts.fit(Mainf, Maincm)
  Mainf <- eBayes(Mainf)
  
  assign("vdata", data.frame(vdat), envir = .GlobalEnv)
  assign("Mainf",  Mainf, envir = .GlobalEnv)
  #assign("cts", data.frame(cts[["counts"]]), envir = .GlobalEnv)
  assign("pct", pct, envir = .GlobalEnv)
  
  for(i in colnames(Mainf$coefficients)){
    limma_resfun(coef=i, celltype=celltype)
  }
}

limma_resfun <- function(coef, celltype){
  #Small function to extract coefficients (and save as csv) and display a p value distribution
  #Fill in coef = whatever the coefficient is
  #Expects to receive Mainf from global environment
  cutoff = 1
  df.1a <- topTable(Mainf, coef = coef, p.value=cutoff, n=length(Mainf[["Amean"]]))
  df.1a <- df.1a[order(row.names(df.1a)), ]
  vdata <- vdata[order(row.names(vdata)), ]
  df.1a <- cbind(df.1a, vdata)
  
  path <- paste(curr_dir, "/limma/", coef, "/", sep="") #specify dir
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  hist(df.1a[,"P.Value"], main = paste(coef, celltype), ylab = "Number of genes", xlab = "P value")
  
  limma::plotMA(Mainf, coef = coef, main = coef)
  
  df.1b <-df.1a
  df.1b$diffexpressed <- "NS"
  df.1b$diffexpressed[df.1b$logFC > LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "UP"
  df.1b$diffexpressed[df.1b$logFC < -LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "DOWN"
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NS")
  
  plts <- ggplot(data=df.1b, aes(x=logFC, y=-log10(adj.P.Val), col = diffexpressed)) +
    geom_point() + theme_minimal() +
    scale_color_manual(values=mycolors) +
    geom_vline(xintercept=c(-LFcutoff, LFcutoff), col="red") +
    geom_hline(yintercept=-log10(FDRcutoff), col="red")
  plot(plts)
  graph2png(plts, paste(path, celltype, "_",  coef, "_Volcano.png", sep=""), width=4, height=2.5)
  
  df.2a <- df.1a %>% rownames_to_column("gene")
  df.2a <- left_join(df.2a, Gene.search, by = "gene")
  df.2a <- left_join(df.2a, pct, by = "gene")
  #cts <- cts %>% rownames_to_column("gene")
  #df.2a <- merge(df.2a, cts, by = "gene")
  df.2a <- df.2a[order(df.2a$adj.P.Val),, drop = FALSE]
  write.csv(df.2a, file = paste(path, "Spreadsheet_", celltype, "_",  coef, "_uncorrected.csv", sep = ""))
  assign(c(paste("limma", coef, celltype, sep = "_")), df.2a, envir = .GlobalEnv) #Save in Global environment
  
  df.2c <- subset(df.2a, adj.P.Val < FDRcutoff)
  df.2c <- subset(df.2c, abs(logFC) > LFcutoff)
  if(nrow(df.2c)>1){
    write.csv(df.2c, file = paste(path, "Spreadsheet_", celltype, "_",  coef, "_corrected.csv", sep = ""))
    assign(c(paste("limma", coef, celltype, "sig", sep = "_")), df.2c, envir = .GlobalEnv) #Save in Global environment
  }
}

#pvalue scatterplot function requires the toptable dfs for both male and female and a filename.
PvalScat_fun <- function(datA, datB, filename, xlab, ylab){
  #Genes in one in other !!!!(this may throw out some genes only expressed in one)
  datA <- datA[order(datA$gene),]
  datB <- datA[order(datB$gene),]
  datA <- datA[rownames(datA) %in% rownames(datB),]
  datB <- datB[rownames(datB) %in% rownames(datA),]
  
  #Sorting and log transforming
  pval <- data.frame(datA$adj.P.Val, datB$adj.P.Val)
  pval <- log10(pval)
  pval <- abs(pval)
  colnames(pval) <- c("pdatA","pdatB")
  LFDR <- abs(log10(FDRcutoff))
  
  #Determining which are significant in both or either
  pval <- mutate(pval, AandB = ifelse(pval$pdatA >= LFDR & pval$pdatB >= LFDR, 1, 0))
  pval <- mutate(pval, AorB = ifelse(pval$pdatA >= LFDR | pval$pdatB >= LFDR, 1, 0))
  pval <- mutate(pval, Sums = rowSums(pval[, c(3,4)]))
  
  #plot for use in illustrator
  plts <- ggplot(pval, aes(x = pdatA, y = pdatB)) +
    xlim(0,15) + ylim(0,15) +
    geom_point(aes(colour = cut(Sums, c(-Inf, 0, 1, 1.9, Inf))), show.legend = FALSE, size = 2, shape = 15) +
    scale_colour_manual(name = "Sums",
                        values = c("(-Inf,0]" = "black",
                                   "(0,1]" = "darkred",
                                   "(1.9, Inf]" = "red")) +
    coord_fixed(ratio = 1) +
    theme_void()
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/limma/P Value scatterplot", filename, ".pdf", sep=""))
}


tcutoff=3.2

#logfc scatterplot function requires the toptable dfs for both male and female and a filename.
tScat_fun <- function(datA, datB, filename){
  #Genes in one in other !!!!(this may throw out some genes only expressed in one)
  datA <- datA[order(datA$gene),]
  datB <- datA[order(datB$gene),]
  datA <- datA[rownames(datA) %in% rownames(datB),]
  datB <- datB[rownames(datB) %in% rownames(datA),]
  
  #Sorting and log transforming
  tval <- data.frame(datA$t, datB$t)
  colnames(tval) <- c("ldatA","ldatB")
  
  #Determining which are significant in both or either
  tval <- mutate(tval, AandB = ifelse(tval$ldatA >= tcutoff & tval$ldatB >= tcutoff, 1, 0) + ifelse(tval$ldatA <= -tcutoff & tval$ldatB <= -tcutoff, 1, 0))
  tval <- mutate(tval, AorB = ifelse(tval$ldatA >= tcutoff | tval$ldatB >= tcutoff, 1, 0) + ifelse(tval$ldatA <= -tcutoff | tval$ldatB <= -tcutoff, 1, 0))
  tval <- mutate(tval, Sums = rowSums(tval[, c(3,4)]))
  
  
  #plot for use in illustrator
  plts <- ggplot(tval, aes(x = ldatA, y = ldatB)) +
    xlim(-10,10) + ylim(-10,10) +
    geom_point(aes(colour = cut(Sums, c(-Inf, 0, 1, 1.9, Inf))), show.legend = FALSE, size = 1) +
    scale_colour_manual(name = "Sums",
                        values = c("(-Inf,0]" = "black",
                                   "(0,1]" = "darkred",
                                   "(1.9, Inf]" = "red")) +
    coord_fixed(ratio = 1) +
    theme_void(base_size=7) 
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/limma/t scatterplot", filename, ".pdf", sep=""))
}







###Pattern
library(clusterProfiler)

#cutoff value can be changed as required
Pvalcutoff <- 0.05
LFcutoff <- 0.25

#cor.test for df
cor.test.df <- function(x, Template){
  k<-cor.test(x,Template)$p.value
}

#function to plot selected genes for Cog_Path and NCI, MCI, AD variable only
Pattern_plot_fun_n <- function(TM.dat){
  rownames(TM.dat) <- NULL
  avg_weights_TM <- column_to_rownames(TM.dat, "gene")
  avg_weights_TM$TM_select <- NULL
  avg_weights_TM$t <- NULL
  avg_weights_TM$Correlation <- NULL
  #Take mean per group
  avg_weights_TM <- t(avg_weights_TM)
  temp <- mdat[rownames(mdat) %in% rownames(avg_weights_TM),]
  temp <- list(temp$Cog_Path)
  avg_weights_TM <- aggregate(avg_weights_TM[, 1:(ncol(avg_weights_TM)-1)], temp, mean)
  avg_weights_TM <- avg_weights_TM %>% column_to_rownames("Group.1") %>% scale() %>% t() %>% data.frame() %>% rownames_to_column("gene") %>% pivot_longer(cols=c(NCI, MCI, AD), names_to="Group") %>% mutate(Group = fct_relevel(Group, "NCI", "MCI", "AD"))
  
  df_summary <- data.frame(Group=avg_weights_TM$Group, n=tapply(avg_weights_TM$value, avg_weights_TM$Group, length), mean=tapply(avg_weights_TM$value, avg_weights_TM$Group, mean))
  df_summary$up <- df_summary$mean+df_summary$n/3/500
  df_summary$low <- df_summary$mean-df_summary$n/3/500
  
  plts <- ggplot(df_summary, aes(x=Group, y=mean, group=1)) +
    geom_ribbon(aes(ymin=up, ymax=low), fill="blue", alpha=0.5) +
    geom_line(data=avg_weights_TM, aes(Group, value, group=gene), color="gray", alpha=0.5) +
    theme_classic() + ylab("Expression")
  plot(plts)
  return(plts)
}

Pattern_plot_fun <- function(TM.dat){
  rownames(TM.dat) <- NULL
  avg_weights_TM <- column_to_rownames(TM.dat, "gene")
  avg_weights_TM$TM_select <- NULL
  avg_weights_TM$t <- NULL
  avg_weights_TM$Correlation <- NULL
  #Take mean per group
  avg_weights_TM <- t(avg_weights_TM)
  temp <- mdat[rownames(mdat) %in% rownames(avg_weights_TM),]
  temp <- list(temp$Cog_Path)
  avg_weights_TM <- aggregate(avg_weights_TM[, 1:(ncol(avg_weights_TM)-1)], temp, mean)
  avg_weights_TM <- avg_weights_TM %>% column_to_rownames("Group.1") %>% scale() %>% t() %>% data.frame() %>% rownames_to_column("gene") %>% pivot_longer(cols=c(NCI, MCI, AD), names_to="Group") %>% mutate(Group = fct_relevel(Group, "NCI", "MCI", "AD"))
  
  df_summary <- data.frame(Group=avg_weights_TM$Group, n=tapply(avg_weights_TM$value, avg_weights_TM$Group, length), mean=tapply(avg_weights_TM$value, avg_weights_TM$Group, mean))
  df_summary$sd <- tapply(avg_weights_TM$value, avg_weights_TM$Group, sd)
  df_summary$sem <- df_summary$sd/sqrt(df_summary$n-1)
  df_summary$CI_lower <- df_summary$mean + qt((1-0.95)/2, df=df_summary$n-1)*df_summary$sem
  df_summary$CI_upper <- df_summary$mean - qt((1-0.95)/2, df=df_summary$n-1)*df_summary$sem
  
  plts <- ggplot(df_summary, aes(x=Group, y=mean, group=1)) +
    geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper), fill="blue", alpha=0.5) +
    geom_line(data=avg_weights_TM, aes(Group, value, group=gene), color="gray", alpha=0.5) +
    theme_classic() + ylab("Expression")
  plot(plts)
  return(plts)
}

#Small function to make sure the sorted data matches the desired template. Enter TM (sorted data), Sex ("M" or "F"), Template (as above)
Cor.fun <- function(TM_select, Template, tval){
  #Extract significant limma DEG results
  nameTM <- deparse(substitute(TM_select))
  
  colnames(tval) <- c("gene", "t")
  tval <- column_to_rownames(tval, "gene")
  TM.dat <- vdatweights #make sure you have vdatweights in your environment
  TM.dat <- TM.dat[order(match(rownames(TM.dat),rownames(tval))),]
  
  #Run correlation
  Correlation <- apply(TM.dat, 1, cor.test.df, Template)
  Correlation <- data.frame(Correlation)
  Cor.Hist <- hist(Correlation$Correlation, main = nameTM)
  pdf(paste0(curr_dir, "/limma/Pattern/P Value distribution correlation_", nameTM, ".pdf"))
  plot(Cor.Hist)
  dev.off()
  
  #merge dataframes
  TM.dat <- cbind(TM.dat, tval, Correlation)
  TM.dat$t <- TM.dat$t*(1-TM.dat$Correlation) #keep t value for fgsea and account for correlation
  
  TM.dat <- rownames_to_column(TM.dat, "gene")
  TM.dat <- TM.dat[order(TM.dat$t),, drop=FALSE]
  assign(paste0("Pattern_", nameTM), TM.dat, .GlobalEnv)
  write.csv(TM.dat, paste0(curr_dir, "/limma/Pattern/P Value distribution correlation_", nameTM, ".csv"))
  
  #Select only significant
  TM.dat <- TM.dat[TM.dat$Correlation<CORcutoff,]
  TM_select <- data.frame(TM_select)
  TM.dat <- subset(TM.dat, TM.dat$gene %in% rownames(TM_select))
  
  TM.dat_pos <- TM.dat[TM.dat$t>0,]
  if(nrow(TM.dat_pos)>2){
    plts <- Pattern_plot_fun(TM.dat_pos)
    graph2pdf(plts, paste0(curr_dir, "/limma/Pattern/P Value distribution correlation", nameTM, "_pos.pdf"))
    plts <- Pattern_plot_fun_n(TM.dat_pos)
    graph2pdf(plts, paste0(curr_dir, "/limma/Pattern/P Value distribution correlation N", nameTM, "_pos.pdf"))
  }
  
  TM.dat_neg <- TM.dat[TM.dat$t<0,]
  if(nrow(TM.dat_neg)>2){
    plts <- Pattern_plot_fun(TM.dat_neg)
    graph2pdf(plts, paste0(curr_dir, "/limma/Pattern/P Value distribution correlation", nameTM, "_neg.pdf"))
    plts <- Pattern_plot_fun_n(TM.dat_neg)
    graph2pdf(plts, paste0(curr_dir, "/limma/Pattern/P Value distribution correlation N", nameTM, "_neg.pdf"))
  }
  
  assign(paste0("Pattern_", nameTM, "_sig"), TM.dat, .GlobalEnv)
  write.csv(TM.dat, paste0(curr_dir, "/limma/Pattern/P Value distribution correlation", nameTM, "_corrected.csv"))
}


Volcano_fun <- function(TT_res, genes){
  name <- deparse(substitute(TT_res))
  name <- str_sub(name, 7)
  df.1b <- TT_res
  df.1b <- df.1b[!duplicated(df.1b$gene),]
  rownames(df.1b) <- NULL
  df.1b <- df.1b %>% column_to_rownames("gene")
  df.1b$diffexpressed <- "NS"
  df.1b$diffexpressed[df.1b$logFC > LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "UP"
  df.1b$diffexpressed[df.1b$logFC < -LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "DOWN"
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NS")
  
  lbl_gene <- df.1b[rownames(df.1b) %in% genes,]
  
  plts <- ggplot(data=df.1b, aes(x=logFC, y=-log10(adj.P.Val), col = diffexpressed)) +
    ggrastr::rasterise(geom_point(size=.1), dpi=600) + theme_bw() + xlim(-2,2) + ylim(0,8) + 
    scale_color_manual(values=mycolors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 7), axis.title=element_text(size=7)) +
    geom_vline(xintercept=c(-LFcutoff, LFcutoff), col="red") +
    geom_hline(yintercept=-log10(FDRcutoff), col="red") + geom_text_repel(data=lbl_gene, mapping=aes(label=rownames(lbl_gene)), size=2, max.overlaps=1000)
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/limma/volcano_", name, ".pdf", sep=""), width=1.574, height=1.574)
}

tScat_fun <- function(datA, datB, genes=NULL){
  nameA <- deparse(substitute(datA))
  nameA <- str_sub(nameA, 7)
  nameB <- deparse(substitute(datB))
  nameB <- str_sub(nameB, 7)
  
  cutoff <- tcutoff
  
  #Genes in one in other
  datA <- datA[!duplicated(datA$gene),]
  datB <- datB[!duplicated(datB$gene),]
  datA <- data.frame(datA$gene, datA$t)
  datB <- data.frame(datB$gene, datB$t)
  colnames(datA) <- c("gene", "datA")
  colnames(datB) <- c("gene", "datB")
  
  #Sorting and log transforming
  val <- merge(datA, datB, by="gene", all.x=T, all.y=T)
  val[is.na(val)] <- runif(sum(is.na(val)), min = -cutoff, max = cutoff)
  val <- column_to_rownames(val, "gene")
  
  #Determining which are significant in both or either
  val <- mutate(val, AandB = ifelse(val$datA >= cutoff & val$datB >= cutoff, 1, 0) + ifelse(val$datA <= -cutoff & val$datB <= -cutoff, 1, 0))
  val <- mutate(val, AorB = ifelse(val$datA >= cutoff | val$datB >= cutoff, 1, 0) + ifelse(val$datA <= -cutoff | val$datB <= -cutoff, 1, 0))
  val <- mutate(val, Sums = rowSums(val[, c(3,4)]))
  
  lbl_gene <- val[rownames(val) %in% genes,]
  
  #plot for use in illustrator
  plts <- ggplot(val, aes(x = datA, y = datB)) +
    xlim(-8,8) + ylim(-8,8) +
    xlab(nameA) + ylab(nameB) +
    geom_point(aes(colour = cut(Sums, c(-Inf, 0, 1, 1.9, Inf))), show.legend = FALSE, size = .1) +
    scale_colour_manual(name = "Sums",
                        values = c("(-Inf,0]" = "black",
                                   "(0,1]" = "darkred",
                                   "(1.9, Inf]" = "red")) +
    coord_fixed(ratio = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 7), axis.title=element_text(size=7)) +
    geom_vline(xintercept=c(-tcutoff, tcutoff), col="red") +
    geom_hline(yintercept=c(-tcutoff, tcutoff), col="red") +
    geom_text_repel(data=lbl_gene, mapping=aes(label=rownames(lbl_gene)), size=2, seed = 42, box.padding = 0.5, max.overlaps = Inf, arrow = arrow(length = unit(0.010, "npc")))
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/limma/scatterplot_", nameA, "_", nameB, "tval.pdf", sep=""))
  
  assign(paste(nameA, nameB, sep="_"), val, envir = .GlobalEnv)
}

box_fun <- function(limma_res, metadat, colour){
  #Takes a column of gene expression values and merges it with metadata from mdat in the environment to produce plot and return a box + point plot.
  #limma result (will trim the first 8 and last 2 cols since those are usually stats)
  #metadat: Column from mdat to aggregate data
  #Set of colours to use for box plot
  #Select relevant rows from Design Matrix 
  df <- data.frame(t(limma_res[limma_res$gene==gene, c(8:(length(colnames(limma_res))-2))])) #!!!!Trims the first 8 and last 2 cols since those are stats
  colnames(df) <- "Expression"
  dm <- matrix(ncol=ncol(mdat))
  dm <- dm[-1,]
  colnames(dm) <- colnames(mdat)
  
  for(i in rownames(df)){
    base_name <- sapply(str_split(i,"_",),'[',1) #get projid names
    base_name <- str_remove(base_name, "X")
    if(base_name %in% rownames(mdat)){
      name_rows <- rownames(dm)
      row_add <- mdat[base_name,]
      dm <- rbind(dm, row_add)
      rownames(dm) <- c(name_rows, paste(i))
    }else{print(paste("Warning colname of data:", i, "not in rowname of Design Matrix"))}
  }
  
  df <- merge(df, dm, by='row.names')
  df <- data.frame(df$Expression, df[[metadat]])
  colnames(df) <- c("Expression", "Group")
  
  df <- df[order(df$Group, decreasing=T),]
  df$Group <- factor(df$Group, levels=unique(df$Group))
  plts <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
    geom_boxplot(fatten=1, colour="black") +
    scale_fill_manual(values=colour) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), legend.position="none", axis.line = element_line(colour = "black"), text=element_text(size=7, color="black"),
          axis.title.x=element_blank(), axis.title.y=element_text(face="italic", size=7)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    ylab(gene)
  plot(plts)
  return(plts)
}


Path_fun <- function(TT_res, invert=F, scoreType="std"){
  #TT_res: limma top table result (not cutoff) for fgsea to analyse
  #invert: invert direction of change
  path <- paste(curr_dir, "/FGSEA/", sep="") #specify dir
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  name <- deparse(substitute(TT_res))
  name <- str_sub(name, 7)
  
  DEG <- TT_res %>% dplyr::select(gene, t) %>%  na.omit() %>% distinct() %>% deframe()
  if(invert==T){DEG <- -DEG}
  
  fgseaRes <- fgsea(pathways = pathway, stats = DEG, minSize = 10, scoreType = scoreType)
  fgseaResSig <- subset(fgseaRes, padj < 0.1 & log2err<0.6)
  
  #plts <- ggplot(fgseaResSig, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=NES<0)) + coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA") + theme_minimal() + theme(text = element_text(size = 3))
  #graph2pdf(plts, paste(path, name, ".pdf", sep=""), width=4, height=2.5)
  #plot(plts)
  
  #label pathways and save
  fgseaResSig <- fgseaResSig[order(fgseaResSig$padj),, drop=FALSE]
  fgseaResSig <- right_join(GO_labelled, fgseaResSig, by="pathway")
  write.csv(apply(fgseaResSig, 2, as.character), file=paste(path, name, ".csv", sep=""))
  
  fgseaRes <- fgseaRes[order(fgseaRes$padj),, drop=FALSE]
  fgseaRes <- right_join(GO_labelled, fgseaRes, by="pathway")
  write.csv(apply(fgseaRes, 2, as.character), file=paste(path, name, "_all.csv", sep=""))
  assign(paste("GSEA", name, "sig", sep="_"), fgseaResSig, envir = .GlobalEnv)
  assign(paste("GSEA", name, sep="_"), fgseaRes, envir = .GlobalEnv)
}


library(purrr)
library(GSVA)
library(msigdbr)
library(stringr)

FDRcutoff <- 0.01
LFcutoff <- 0.25

GSVA_fun <- function(limma_res, pathways){
  #any limma result will do, you want a normalized expression matrix
  name <- deparse(substitute(limma_res))
  name <- sub(".*_", "", name)
  
  limma_res <- limma_res[!duplicated(limma_res$gene),]
  rownames(limma_res) <- limma_res$gene
  exp_matrix <- limma_res[, grepl("^X[0-9]{8}", colnames(limma_res))]
  pathway_scores <- GSVA::gsva(expr = as.matrix(exp_matrix), gset.idx.list = pathways, method = "gsva", kcdf = "Gaussian", parallel.sz = 32, min.sz=1)
  colnames(pathway_scores) <- sub("^X", "", colnames(pathway_scores))
  colnames(pathway_scores) <- sub("\\.", "-", colnames(pathway_scores))
  write.csv(pathway_scores, file = paste(curr_dir, "/GSVA_modules/GSVA_snRNA_", name, ".csv", sep = ""))
  assign(paste("GSVA", name, sep="_"), pathway_scores, envir = .GlobalEnv)
}

GSVA_limmafun_lane <- function(cts){
  #cts is a GSVA_matrix
  celltype <- deparse(substitute(cts))
  
  #Select relevant rows from Design Matrix  
  dm <- matrix(ncol=ncol(Maindm))
  dm <- dm[-1,]
  colnames(dm) <- colnames(Maindm)
  
  for(i in colnames(cts)){
    base_name <- sapply(str_split(i,"_",),'[',1) #get projid names
    if(base_name %in% rownames(Maindm)){
      name_rows <- rownames(dm)
      row_add <- Maindm[base_name,]
      dm <- rbind(dm, row_add)
      rownames(dm) <- c(name_rows, paste(i))
    }else{print(paste("Warning colname of data:", i, "not in rowname of Design Matrix"))}
  }
  
  #Enforce Maindm and cts name match
  if(all(rownames(dm)==colnames(cts))==FALSE){
    dm <- dm[rownames(dm) %in% colnames(cts),]      
    cts <- cts[,colnames(cts) %in% rownames(dm)]
  }
  
  #order dm and cts based on first variable of interest
  dm <- dm[order(dm[,1]),]
  cts <- cts[,order(match(colnames(cts), rownames(dm)))]
  
  projid <- factor(sapply(str_split(colnames(cts),"_",),'[', 1)) #get projid names
  corfit <- duplicateCorrelation(cts, block=projid)
  
  #limma
  Mainf <- lmFit(cts, dm, block=projid, correlation=corfit$consensus.correlation)
  Mainf <- contrasts.fit(Mainf, Maincm)
  Mainf <- eBayes(Mainf, trend = TRUE) #not Homoskedastic among others
  
  assign("vdata", data.frame(cts), envir = .GlobalEnv)
  assign("Mainf",  Mainf, envir = .GlobalEnv)
  
  for(i in colnames(Mainf$coefficients)){
    limma_resfun(coef=i, celltype=celltype)
  }
}


library(ggplot2)
library(ggrepel)
library(ggpubr)

z_cutoff=.25

#z scatterplot function 
zScat_fun <- function(datA, datB){
  nameA <- deparse(substitute(datA))
  nameA <- str_sub(nameA, 7)
  nameB <- deparse(substitute(datB))
  nameB <- str_sub(nameB, 7)
  
  cutoff=z_cutoff
  
  #Genes in one in other
  datA <- datA[!duplicated(datA$gene),]
  datB <- datB[!duplicated(datB$gene),]
  datA <- data.frame(datA$gene, datA$effect_z, datA$CI.R)
  datB <- data.frame(datB$gene, datB$effect_z, datB$CI.R)
  colnames(datA) <- c("gene", "datA", "xCIL")
  colnames(datB) <- c("gene", "datB", "yCIL")
  
  #add CI to effect score
  datA$xCIH <- datA$datA-datA$xCIL
  datA$xCIL <- datA$datA+datA$xCIL
  datB$yCIH <- datB$datB-datB$yCIL
  datB$yCIL <- datB$datB+datB$yCIL
  
  #Sorting and log transforming
  val <- merge(datA, datB, by="gene", all.x=F, all.y=F)
  val[is.na(val)] <- runif(sum(is.na(val)), min = -cutoff, max = cutoff)
  val <- column_to_rownames(val, "gene")
  
  #Determining which are significant in both or either
  val <- mutate(val, AandB = ifelse(val$datA >= cutoff & val$datB >= cutoff, 1, 0) + ifelse(val$datA <= -cutoff & val$datB <= -cutoff, 1, 0))
  val <- mutate(val, AorB = ifelse(val$datA >= cutoff | val$datB >= cutoff, 1, 0) + ifelse(val$datA <= -cutoff | val$datB <= -cutoff, 1, 0))
  val <- mutate(val, Sums = rowSums(val[, c(3,4)]))
  
  
  #plot for use in illustrator
  plts <- ggplot(val, aes(x=datA, y=datB)) +
    xlim(-2,2) + ylim(-2,2) +
    xlab(nameA) + ylab(nameB) +
    geom_errorbar(aes(xmin=xCIL, xmax=xCIH), orientation="y", height=.01, show.legend=FALSE) +
    geom_errorbar(aes(ymin=yCIL, ymax=yCIH), width=.01, show.legend=FALSE) +
    geom_point(show.legend = FALSE, size = 1) +
    coord_fixed(ratio = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 7), axis.title=element_text(size=7)) +
    geom_vline(xintercept=c(-cutoff, cutoff), col="red") +
    geom_hline(yintercept=c(-cutoff, cutoff), col="red") +
    geom_text_repel(data=val, mapping=aes(label=rownames(val)), size=2, seed = 42, box.padding = 0.5, max.overlaps = Inf, arrow = arrow(length = unit(0.010, "npc"))) +
    geom_smooth(method = "lm", col = "black", se = FALSE) +
    stat_cor(method = "pearson")
  plot(plts)
  graph2pdf(plts, paste0(curr_dir, "/GSVA_modules/scatterplot_", nameA, "_", nameB, "_Zval.pdf"))
  
  assign(paste(nameA, nameB, sep="_"), val, envir = .GlobalEnv)
}


pcutoff <- 0.05
LFcutoff <- 0.001

TFlimmafun_lane <- function(cts){
  celltype <- deparse(substitute(cts))
  celltype <- substr(celltype,5,30)
  
  #Select relevant rows from Design Matrix  
  dm <- matrix(ncol=ncol(Maindm))
  dm <- dm[-1,]
  colnames(dm) <- colnames(Maindm)
  
  #Do not voom or scale, data is already normalized  cts <- t(scale(t(cts))) #Scale AUC scores
  
  for(i in colnames(cts)){
    base_name <- sapply(str_split(i,"_",),'[',1) #get projid names
    if(base_name %in% rownames(Maindm)){
      name_rows <- rownames(dm)
      row_add <- Maindm[base_name,]
      dm <- rbind(dm, row_add)
      rownames(dm) <- c(name_rows, paste(i))
    }else{print(paste("Warning colname of data:", i, "not in rowname of Design Matrix"))}
  }
  
  #Enforce Maindm and cts name match
  if(all(rownames(dm)==colnames(cts))==FALSE){
    dm <- dm[rownames(dm) %in% colnames(cts),]      
    cts <- cts[,colnames(cts) %in% rownames(dm)]
  }
  
  #order dm and cts based on first variable of interest
  dm <- dm[order(dm[,1]),]
  cts <- cts[,order(match(colnames(cts), rownames(dm)))]
  
  #technical replicate
  projid <- factor(sapply(str_split(colnames(cts),"_",),'[', 1)) #get projid names
  corfit <- duplicateCorrelation(cts, block=projid)
  
  #limma
  Mainf <- lmFit(cts, dm, block=projid, correlation=corfit$consensus.correlation)
  Mainf <- contrasts.fit(Mainf, Maincm)
  Mainf <- eBayes(Mainf, trend=T)
  
  assign("vdata", data.frame(cts), envir = .GlobalEnv)
  assign("Mainf",  Mainf, envir = .GlobalEnv)
  
  for(i in colnames(Mainf$coefficients)){
    limma_resfun(coef=i, celltype=celltype)
  }
}

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
mdat$projid <- rownames(mdat)

heatfun <- function(Celltype){
  lmr1 <- paste0("limma_cogdxadnc_resilientvNCI_TF_", Celltype, "_sig")
  lmr2 <- paste0("limma_cogdxadnc_ADdvNCI_TF_", Celltype, "_sig")
  cts <- paste0("cts_TF_", Celltype)
  
  # Fetch the actual objects from the R environment using get()
  lmr1 <- get(lmr1)
  lmr2 <- get(lmr2)
  cts <- get(cts)
  
  sig_genes <- c(unique(c(lmr1$gene, lmr2$gene)))
  print(paste0("Sig genes for ", Celltype, ": "))
  print(sig_genes)
  mat <- cts[rownames(cts) %in% sig_genes,]
  mat <- data.frame(t(mat))
  mat$projid <- substr(rownames(mat), 1, 8)
  mat <- merge(mat, data.frame("projid"=mdat$projid, "Cog_Path"=mdat$Cog_Path), by = "projid", all.x = TRUE)
  mat$projid <- NULL
  mat <- mat %>% group_by(Cog_Path) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  mat <- subset(mat, Cog_Path=="NCI" | Cog_Path=="resilient" | Cog_Path=="ADd") %>% as.data.frame()
  rownames(mat) <- mat$Cog_Path
  mat$Cog_Path <- NULL
  mat <- t(scale(mat))
  mat <- mat[, c("NCI", "resilient", "ADd")]
  plts <- Heatmap(mat, cluster_rows = T, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="Transcription Factor AUC", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), col = colorRamp2(c(-1, 0, 1), c("blue","white", "red")))
  plot(plts)
  graph2pdf(plts, file=paste(curr_dir, "/Differential_Regulon/", Celltype, "_expression_heatmap.pdf", sep = ""), font="Arial", width=4, height=7)
}

library(ggraph)
library(igraph)
library(ggnewscale)
library(tidygraph)

edges <- read_csv("/Ex_ChP/Ex_ChP7/Scenic/results/SCENIC_regulon_edges.csv")
edges$NES <- NULL
edges$Regulon <- NULL

network_plot <- function(lmTF, lmdeg, tflist, name){
  lmTF <- data.frame("gene"=lmTF$gene, "logFC"=lmTF$regulon_activity)
  lmdeg <- data.frame("gene"=lmdeg$gene, "logFC"=lmdeg$logFC)
  
  edges_cell <- edges[edges$TF %in% tflist,]
  edges_cell <- edges_cell[edges_cell$Target %in% c(lmdeg$gene),]
  edges_cell <- edges_cell[edges_cell$TF %in% c(lmTF$gene),]
  edges_cell <- edges_cell[edges_cell$Weight>1,]
  edge_TF <- edges_cell[edges_cell$Target %in% tflist,]
  edges_cell <- edges_cell[order(edges_cell$Weight, decreasing = TRUE)[1:30],]
  edges_cell <- rbind(edges_cell, edge_TF)
  edges_cell <- edges_cell[!is.na(edges_cell$Weight),]
  
  nodes <- tibble(name = unique(c(edges_cell$TF, edges_cell$Target))) %>%
    left_join(lmTF, by=c("name"="gene")) %>%
    left_join(lmdeg, by=c("name"="gene"))
  nodes$type <- ifelse(nodes$name %in% edges_cell$TF, "TF", "Target")
  nodes <- nodes[!duplicated(nodes$name),]
  
  graph <- graph_from_data_frame(d = edges_cell, vertices = nodes, directed = TRUE)
  
  plts <- ggraph(graph, layout = "stress") +
    geom_edge_link(aes(width = Weight), alpha = 0.5) + scale_edge_width(range = c(0.1, .5)) +
    geom_node_point(data = . %>% filter(type == "Target"), aes(color = logFC.y), size = 2) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1), oob = scales::squish, name = "gene logFC") +
    new_scale_color() +
    geom_node_point(data = . %>% filter(type == "TF"), aes(color = logFC.x), size = 5) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.01, 0.01), oob = scales::squish, name = "TF ΔAUC") +
    geom_node_text(aes(label = name), repel = TRUE)
  plot(plts)
  graph2pdf(plts, file=paste(curr_dir, "/Differential_Regulon/", name, "_networkplot.pdf", sep = ""), font="Arial", width=5, height=4)
}
