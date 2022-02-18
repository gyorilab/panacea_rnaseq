library(xlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library("org.Mm.eg.db")

# Redo heatmaps by adding annotation col for spec. excitibility related genes,
# RAGS, immune-related

# Work on re-doing volcano by fixing x-axis and revert back arrows


# set wd
setwd('~/gitHub/panacea-rnaseq/cold_slurry_drg/')
source('./scripts/functions.R')

# functions
get_colors <- function(){
  my_colors = brewer.pal(n = 10, name = "RdBu")
  my_colors = colorRampPalette(my_colors)(50)
  my_colors = rev(my_colors)
  return(my_colors)
}



# Read the FPKM counts
fpkm <- read.csv('./counts/FPKM.csv')

# clean FPKM table
fpkm_clean = fpkm %>%
  dplyr::select(GeneID, GeneName, ends_with("FPKM")) %>%
  as.data.frame()

rownames(fpkm_clean) <- fpkm_clean$GeneID
fpkm_clean$GeneID <-  NULL
colnames(fpkm_clean)[2:ncol(fpkm_clean)] <- gsub(".FPKM", "", colnames(fpkm_clean)[2:ncol(fpkm_clean)])
colnames(fpkm_clean)

# boxplot of the fpkm counts
png('./output/FPKM_BOXPLOT.png', width=1500, height=1800, res=150)
print(boxplot(fpkm_clean[,-1], las = 2))
dev.off()

# log2 transform
fpkm_clean[c(2:ncol(fpkm_clean))] <- log2(fpkm_clean[,c(2:ncol(fpkm_clean))] + 1)

png('./output/FPKM_LOG_BOXPLOT.png', width=1500, height=1800, res=150)
print(boxplot(fpkm_clean[, -1], las = 2))
dev.off()

# read human rat orthologue
ortho <- read.csv('./gene_lists/human_rat_ortho.tsv', sep='\t')
ortho <- ortho[, c(2,4)]
colnames(ortho) <- c('RAT_GENES', 'HUMAN_GENES')

# get the differentially expressed genes
diff_genes <- read.csv('./output/DESEQ/diff_CS_ipsi_3D_7D_vs_CS_contra_3D_7D.csv')
diff_genes <- diff_genes$rat_gene_name

# Read annotation file
annot_df <- read.csv('./counts/annotation.csv', header=F, strip.white = T, stringsAsFactors = F)
rownames(annot_df) <- colnames(fpkm_clean)[2:ncol(fpkm_clean)]
colnames(annot_df) <- 'type'
sample_type <- list("CTL_DRG_Con"="naive (C)",
                    "CTL_DRG_Ipsi"="naive (I)",
                    "RT_6h_DRG_Con"="6h RT (C)",
                    "RT_6h_DRG_Ipsi"="6h RT (I)",
                    "CS6h_DRG_Con"="6h CS (C)",
                    "CS6h_DRG_Ipsi"="6h CS (I)",
                    "CS1d_DRG_Con"="1D CS (C)",
                    "CS1d_DRG_Ipsi"="1D CS (I)",
                    "CS3d_DRG_Con"="3D CS (C)",
                    "CS3d_DRG_Ipsi"="3D CS (I)",
                    "CS7d_DRG_Con"="7D CS (C)",
                    "CS7d_DRG_Ipsi"="7D CS (I)",
                    "CS14d_DRG_Con"="14D CS (C)",
                    "CS14d_DRG_Ipsi"="14D CS (I)"
                    )
annot_df$type <- sapply(annot_df$type, FUN = function(x)sample_type[[x]])
annot_df$type <- factor(annot_df$type, levels = unname(unlist(sample_type)))


# read tf
treg <- read.xlsx('./gene_lists/transcription regulator.xlsx', sheetIndex = 1)
colnames(treg)[3] <- 'HUMAN_GENES'
treg <- merge(treg, ortho, by='HUMAN_GENES')
tregs <- fpkm_clean[fpkm_clean$GeneName %in% treg$RAT_GENES,]
tregs <- tregs[!duplicated(tregs$GeneName), ]
rownames(tregs) <- tregs$GeneName
tregs$GeneName <- NULL
x <-  t(scale(t(tregs)))
x_ipsi <- na.omit(x[,str_detect(colnames(x), '.*Ipsi')])
x_contra <- na.omit(x[,str_detect(colnames(x), '.*Con')])
draw_heatmap(x_contra, 'tregs_contra',res=300, height = nrow(x)*26, width=2500)
draw_heatmap(x_ipsi, 'tregs_ipsi',res=300, height = nrow(x)*26, width=2500)
dev.off()

## Plot specific genes
rag_specific <- read.xlsx('./gene_lists/specific_genes/specificRAGs.xlsx', sheetIndex = 1, header = F)
rag_specific$X1 <- str_to_title(rag_specific$X1) 
rag_specific_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% rag_specific$X1, ]
rownames(rag_specific_fpkm) <- rag_specific_fpkm$GeneName
rag_specific_fpkm$GeneName <- NULL
x <- t(scale(t(rag_specific_fpkm)))
x <- na.omit(x)
draw_heatmap(x, 'rag_specific_targets',res=300, height = nrow(x)*26, width=2500,
             annot_df)

targets <- read.xlsx('./gene_lists/specific_genes/specificAtf3targets.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
annot_df <- data.frame('all_samples' = colnames(x),
                       'samples' = factor(annot_df$V1)) %>% column_to_rownames('all_samples')
draw_heatmap(x, 'atf3_specific_targets',res=300, height = nrow(x)*29, width=2300)

targets <- read.xlsx('./gene_lists/specific_genes/specificExcitabilitytargets.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
targets_fpkm <- targets_fpkm[!duplicated(targets_fpkm$GeneName),]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
draw_heatmap(x, 'excitability_specific_targets',res=300, height = nrow(x)*30, width=2300,
             annot_df)


targets <- read.xlsx('./gene_lists/specific_genes/specificGPCRs.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
targets_fpkm <- targets_fpkm[!duplicated(targets_fpkm$GeneName),]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
annot_df <- data.frame('all_samples' = colnames(x),
                       'samples' = factor(annot_df$V1)) %>% column_to_rownames('all_samples')
draw_heatmap(x, 'gpcr_specific_targets',res=300, height = nrow(x)*32, width=2300)

targets <- read.xlsx('./gene_lists/specific_genes/specificimmunegenesxlsx.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
targets_fpkm <- targets_fpkm[!duplicated(targets_fpkm$GeneName),]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
draw_heatmap(x, 'immune_specific_targets',res=300, height = nrow(x)*80, width=3000,
             annot_df)


targets <- read.xlsx('./gene_lists/specific_genes/specificionchannels.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
targets_fpkm <- targets_fpkm[!duplicated(targets_fpkm$GeneName),]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
annot_df <- data.frame('all_samples' = colnames(x),
                       'samples' = factor(annot_df$V1)) %>% column_to_rownames('all_samples')
draw_heatmap(x, 'ion_channel_specific_targets',res=300, height = nrow(x)*32, width=2300)


targets <- read.xlsx('./gene_lists/specific_genes/specifickinases.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
targets_fpkm <- targets_fpkm[!duplicated(targets_fpkm$GeneName),]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
annot_df <- data.frame('all_samples' = colnames(x),
                       'samples' = factor(annot_df$V1)) %>% column_to_rownames('all_samples')
draw_heatmap(x, 'kinases_specific_targets',res=300, height = nrow(x)*32, width=3300)


targets <- read.xlsx('./gene_lists/specific_genes/specificRAGs.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
targets_fpkm <- targets_fpkm[!duplicated(targets_fpkm$GeneName),]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
annot_df <- data.frame('all_samples' = colnames(x),
                       'samples' = factor(annot_df$V1)) %>% column_to_rownames('all_samples')
draw_heatmap(x, 'rags_specific_targets',res=300, height = nrow(x)*32, width=3300, annot_df)


targets <- read.xlsx('./gene_lists/specific_genes/specificTubulinlinkedgenes.xlsx', sheetIndex = 1, header = F)
targets$X1 <- str_to_title(targets$X1) 
targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% targets$X1, ]
targets_fpkm <- targets_fpkm[!duplicated(targets_fpkm$GeneName),]
rownames(targets_fpkm) <- targets_fpkm$GeneName
targets_fpkm$GeneName <- NULL
x <- t(scale(t(targets_fpkm)))
x <- na.omit(x)
annot_df <- data.frame('all_samples' = colnames(x),
                       'samples' = factor(annot_df$V1)) %>% column_to_rownames('all_samples')
draw_heatmap(x, 'tubulin_specific_targets',res=300, height = nrow(x)*150, width=2200)

################################################################################

# Read excitibility targets
ex_targets <- read.xlsx('./gene_lists/EXtargets2021-07-17 2.0.xlsx', sheetIndex = 1)
colnames(ex_targets)[1] <- 'HUMAN_GENES'
ex_targets_merged <- merge(ex_targets, ortho, by.x='HUMAN_GENES')
ex_targets_merged <- na.omit(ex_targets_merged[ex_targets_merged$RAT_GENES %in% diff_genes,])
ex_targets_fpkm <- fpkm_clean[fpkm_clean$GeneName %in% ex_targets_merged$RAT_GENES,]
ex_targets_fpkm <- ex_targets_fpkm[!duplicated(ex_targets_fpkm$GeneName), ]
rownames(ex_targets_fpkm) <- ex_targets_fpkm$GeneName
ex_targets_fpkm$GeneName <- NULL
x <- t(scale(t(ex_targets_fpkm)))
x <- na.omit(x)
draw_heatmap(x, 'excitibility_targets', res=300, height = 1000, width=2500)



# read the list of input files
atf3_targets <- read.xlsx('./gene_lists/Targets from ATF3 paper (annotated).xlsx',
                                          sheetIndex = 1)
atf3_targets <- na.omit(atf3_targets[atf3_targets$ID %in% diff_genes,])
atf3_targets <- fpkm_clean[fpkm_clean$GeneName %in% atf3_targets$ID,]
atf3_targets <- atf3_targets[!duplicated(atf3_targets$GeneName),]
rownames(atf3_targets) <- atf3_targets$GeneName
atf3_targets$GeneName <- NULL
x <- t(scale(t(atf3_targets)))
x <- na.omit(x)
draw_heatmap(x, 'atf3_targets',res=300, height = 6500, width=2500)



# common up regulated RAGs
commonupregulatedRAGs_in_DRG <- read.xlsx('./gene_lists/commonupregulatedRAGs in DRG.xlsx',
                                          sheetIndex = 1)
commonupregulatedRAGs_in_DRG <- na.omit(commonupregulatedRAGs_in_DRG[commonupregulatedRAGs_in_DRG$ID %in% diff_genes,])
common_up_reg_rags <- fpkm_clean[fpkm_clean$GeneName %in% commonupregulatedRAGs_in_DRG$ID,]
common_up_reg_rags <- common_up_reg_rags[!duplicated(common_up_reg_rags$GeneName),]
rownames(common_up_reg_rags) <- common_up_reg_rags$GeneName
common_up_reg_rags$GeneName <- NULL
x <- t(scale(t(common_up_reg_rags)))
x <- na.omit(x)
draw_heatmap(x, 'common_up_rags',res=300, height = 1500, width=2500)


# read gprotein  coupled receptors
g_protein <- read.xlsx('./gene_lists/G-protein coupled receptor.xlsx', sheetIndex = 1)
colnames(g_protein)[3] <- 'HUMAN_GENES' 
g_protein <- merge(g_protein, ortho, by.x='HUMAN_GENES')
g_protein <- na.omit(g_protein[g_protein$RAT_GENES %in% diff_genes,])
g_protein <- fpkm_clean[fpkm_clean$GeneName %in% g_protein$RAT_GENES,]
g_protein <- g_protein[!duplicated(g_protein$GeneName),]
rownames(g_protein) <- g_protein$GeneName
g_protein$GeneName <- NULL
x <- t(scale(t(g_protein)))
x <- na.omit(x)
draw_heatmap(x, 'g_protein_coupled_receptor', res=300, height = 1500, width=2500)


# read ion channel
ion_channel <- read.xlsx('./gene_lists/ion channel.xlsx', sheetIndex = 1)
colnames(ion_channel)[3] <- 'HUMAN_GENES'
ion_channel <- merge(ion_channel, ortho, by.x='HUMAN_GENES')
ion_channel <- na.omit(ion_channel[ion_channel$RAT_GENES %in% diff_genes,])
ion_channel <- fpkm_clean[fpkm_clean$GeneName %in% ion_channel$RAT_GENES,]
ion_channel <- ion_channel[!duplicated(ion_channel$GeneName),]
rownames(ion_channel) <- ion_channel$GeneName
ion_channel$GeneName <- NULL
x <- t(scale(t(ion_channel)))
x <- na.omit(x)
draw_heatmap(x, 'ion_channels', res=300, height = 1000, width=2500)

# read kinases
kinases <- read.xlsx('./gene_lists/kinase.xlsx', sheetIndex = 1)
colnames(kinases)[3] <- 'HUMAN_GENES'
kinases <- merge(kinases, ortho, by.x='HUMAN_GENES')
kinases <- na.omit(kinases[kinases$RAT_GENES %in% diff_genes,])
kinases <- fpkm_clean[fpkm_clean$GeneName %in% kinases$RAT_GENES,]
kinases <- kinases[!duplicated(kinases$GeneName),]
rownames(kinases) <- kinases$GeneName
kinases$GeneName <- NULL
x <- t(scale(t(kinases)))
x <- na.omit(x)
draw_heatmap(x, 'kinases', res=300, height = 1000, width=2500)

# read transcription regulators
treg <- read.xlsx('./gene_lists/transcription regulator.xlsx', sheetIndex = 1)
colnames(treg)[3] <- 'HUMAN_GENES'
treg <- merge(treg, ortho, by.x='HUMAN_GENES')
treg <- na.omit(treg[treg$RAT_GENES %in% diff_genes,])
treg <- fpkm_clean[fpkm_clean$GeneName %in% treg$RAT_GENES,]
treg <- treg[!duplicated(treg$GeneName),]
rownames(treg) <- treg$GeneName
treg$GeneName <- NULL
x <- t(scale(t(treg)))
x <- na.omit(x)
draw_heatmap(x, 'tregs', res=300, height = 1200, width=3000)


# read transmembrane receptor
trans_rec <- read.xlsx('./gene_lists/transmembrane receptor.xlsx', sheetIndex = 1)
colnames(trans_rec)[3] <- 'HUMAN_GENES'
trans_rec <- merge(trans_rec, ortho, by.x='HUMAN_GENES')
trans_rec <- na.omit(trans_rec[trans_rec$RAT_GENES %in% diff_genes,])
trans_rec <- fpkm_clean[fpkm_clean$GeneName %in% trans_rec$RAT_GENES,]
trans_rec <- trans_rec[!duplicated(trans_rec$GeneName),]
rownames(trans_rec) <- trans_rec$GeneName
trans_rec$GeneName <- NULL
x <- t(scale(t(trans_rec)))
x <- na.omit(x)
draw_heatmap(x, 'trans_membrane_receptor', res=300, height = 1000, width=2500)

# read immune related genes
immune_genes <- read.xlsx('./gene_lists/immuerelatedgenes (10-27-21).xlsx', 
                          sheetIndex = 1)
colnames(immune_genes)[3] <- 'HUMAN_GENES'
immune_genes <- merge(immune_genes, ortho, by.x='HUMAN_GENES')
immune_genes <- na.omit(immune_genes[immune_genes$RAT_GENES %in% diff_genes,])
immune_genes <- unique(immune_genes$RAT_GENES)
immune_genes <- fpkm_clean[fpkm_clean$GeneName %in% immune_genes,]
immune_genes <- immune_genes[!duplicated(immune_genes$GeneName),]
rownames(immune_genes) <- immune_genes$GeneName
immune_genes$GeneName <- NULL
x <- t(scale(t(immune_genes)))
x <- na.omit(x)
draw_heatmap(x, 'immune_related_genes', res=300, height = 1000, width=2500)


# read tubulin linked genes
human_mouse <- read.csv('./scripts/HOM_MouseHuman.txt', sep='\t')
human_genes <- subset(human_mouse, Common.Organism.Name == 'human')
mouse_genes <- subset(human_mouse, Common.Organism.Name == 'mouse, laboratory')

HOM_genes <- merge(human_genes, mouse_genes, by='DB.Class.Key')
tubulin_linked <- read.xlsx('./gene_lists/Tubilin linked genes (10-26-21).xlsx', sheetIndex = 1)
colnames(tubulin_linked)[1] <- 'HUMAN_GENES'
tubulin_linked <- merge(tubulin_linked, ortho, by='HUMAN_GENES')
tubulin_linked <- na.omit(tubulin_linked[tubulin_linked$RAT_GENES %in% diff_genes,])

tubulin_linked <- fpkm_clean[fpkm_clean$GeneName %in% tubulin_linked$MOUSE_SYMBOL,]
tubulin_linked <- tubulin_linked[!duplicated(tubulin_linked$GeneName),]
rownames(tubulin_linked) <- tubulin_linked$GeneName
tubulin_linked$GeneName <- NULL
x <- t(scale(t(tubulin_linked)))
x <- na.omit(x)
draw_heatmap(x, 'tubulin_linked_genes', res=300, height = 3500, width=2500)
