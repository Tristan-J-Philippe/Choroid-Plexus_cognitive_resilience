###
### 021 Cell_type_proportion_RNA.R
###
# Purpose: Cell type proportion analysis in single nuclei RNA dataset
# Dependencies:
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(1234)
##Cell proportion
####Per subtype

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)

#load up to date meta
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly.csv", col_type=cols(.default="c")))

#cell counts
md <- data.frame(Ex_ChP@meta.data$Subcelltype, Ex_ChP@meta.data$Lane_projid, Ex_ChP@meta.data$projid)
colnames(md) <- c("active.ident", "Lane_projid", "projid")
prop <- tabyl(md, active.ident, Lane_projid)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Neuro", "Astro"),]
parenchyma <- colSums(parenchyma)
prop <- prop[!rownames(prop)%in%c("Neuro", "Astro", "Oligo"),]

#discrete metadata
Ex_ChP_md <- left_join(md, mdat, by="projid")
Ex_ChP_md$active.ident <- NULL
Ex_ChP_md <- distinct(Ex_ChP_md, .keep_all = T)

#Factors
msex = factor(paste(Ex_ChP_md$msex))
X10X_batch = factor(paste(Ex_ChP_md$X10X_batch))
library_batch = factor(paste(Ex_ChP_md$library_batch))
apoe = factor(paste(Ex_ChP_md$apoe_genotype))
cogdx = factor(paste(Ex_ChP_md$cogdx))
Cog_Path = factor(paste(Ex_ChP_md$Cog_Path))
ad_adnc <- factor(paste(Ex_ChP_md$ad_adnc))

#continuous
age = as.numeric(Ex_ChP_md$age_death)
pmi = as.numeric(Ex_ChP_md$pmi)
educ = as.numeric(Ex_ChP_md$educ)
cogn_global = as.numeric(Ex_ChP_md$cogn_global_last)
parenchyma = as.numeric(parenchyma)
tdp_st4 = as.numeric(Ex_ChP_md$tdp_st4)

Lane <- factor(paste(sapply(str_split(colnames(prop), "_"), '[', 2), sapply(str_split(colnames(prop), "_"), '[', 3), sep = "_"))


library(lme4)
library(lmerTest)
library(emmeans)

Ex_ChP_md$Cog_Path <- factor(ifelse(Ex_ChP_md$cogdx == 1 & Ex_ChP_md$ad_adnc==0, "NCI",
                                    ifelse(Ex_ChP_md$cogdx == 1 & Ex_ChP_md$ad_adnc==1, "resilient",
                                           ifelse(Ex_ChP_md$cogdx %in% c(2, 3), "MCI",
                                                  ifelse(Ex_ChP_md$cogdx %in% c(4, 5), "ADd", NA)))),
                             levels = c("NCI", "resilient", "MCI", "ADd"))

cols <- c("contrast", "estimate", "SE", "df", "CI.L", "CI.R", "t.ratio", "p.value")
c_ResvNCI <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(c_ResvNCI) <- cols
c_ADvNCI <- c_ResvNCI
c_ADvRes <- c_ResvNCI

prop_percent <- as.data.frame(proportions(as.matrix(prop), margin=2))*100
Lane <- factor(paste(sapply(str_split(colnames(prop), "_"), '[', 2), sapply(str_split(colnames(prop), "_"), '[', 3), sep = "_"))

for(celltype in rownames(prop_percent)){
  dat <- data.frame(proportion = as.numeric(prop_percent[celltype, ]),
                    projid = Ex_ChP_md$projid,
                    Lane = Lane,
                    Cog_Path = Ex_ChP_md$Cog_Path,
                    msex = Ex_ChP_md$msex,
                    age = as.numeric(Ex_ChP_md$age_death),
                    pmi = as.numeric(Ex_ChP_md$pmi))
  
  fit <- lmer(proportion ~ Cog_Path + msex + age + pmi + parenchyma + (1 | projid) + (1 | Lane), data=dat)
  
  emm <- emmeans(fit, ~Cog_Path, at = list(Cog_Path=c("resilient", "NCI")))
  int <- contrast(emm, interaction = c("pairwise"))
  c_ResvNCI[celltype,] <- data.frame(summary(int, infer = c(TRUE)))
  
  emm <- emmeans(fit, ~Cog_Path, at = list(Cog_Path=c("ADd", "NCI")))
  int <- contrast(emm, interaction = c("pairwise"))
  c_ADvNCI[celltype,] <- data.frame(summary(int, infer = c(TRUE)))
  
  emm <- emmeans(fit, ~Cog_Path, at = list(Cog_Path=c("ADd", "resilient")))
  int <- contrast(emm, interaction = c("pairwise"))
  c_ADvRes[celltype,] <- data.frame(summary(int, infer = c(TRUE)))
}

c_ResvNCI$FDR <- ave(c_ResvNCI$p.value, c_ResvNCI$contrast, FUN = p.adjust, method = "BH")
c_ADvNCI$FDR <- ave(c_ADvNCI$p.value, c_ADvNCI$contrast, FUN = p.adjust, method = "BH")
c_ADvRes$FDR <- ave(c_ADvRes$p.value, c_ADvRes$contrast, FUN = p.adjust, method = "BH")
c_ResvNCI
c_ADvNCI
c_ADvRes

library(ggplot2)
library(patchwork)
plts1 <- ggplot(c_ResvNCI, aes(x = estimate, y = rownames(c_ResvNCI))) + geom_pointrange(aes(xmin = CI.L, xmax = CI.R)) + theme_minimal() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") + labs(x="Res vs NCI", y="Sub Cell Type", title="Percent proportion")
plts2 <- ggplot(c_ADvNCI, aes(x = estimate, y = rownames(c_ADvNCI))) + geom_pointrange(aes(xmin = CI.L, xmax = CI.R)) + theme_minimal() + geom_vline(xintercept = 0, linetype = "dotted", color = "black") + labs(x="AD vs NCI", y="Sub Cell Type", title="Percent proportion")
plts3 <- ggplot(c_ADvRes, aes(x = -estimate, y = rownames(c_ADvRes))) + geom_pointrange(aes(xmin = -CI.L, xmax = -CI.R)) + theme_minimal() + geom_vline(xintercept = 0, linetype = "dotted", color = "black") + labs(x="Res vs AD", y="Sub Cell Type", title="Percent proportion")

plts1+plts2+plts3

ggsave2(paste0(curr_dir, "/CellType/Publish/Proportions/Resilience subtypes.pdf"), plts1+plts2+plts3)
write.csv(c_ResvNCI, paste0(curr_dir, "/CellType/Publish/Proportions/ResvNCI.cvs"))
write.csv(c_ADvNCI, paste0(curr_dir, "/CellType/Publish/Proportions/ADvNCI.csv"))
write.csv(c_ADvRes, paste0(curr_dir, "/CellType/Publish/Proportions/ADvRes.csv"))
