library(dplyr)
library(tidyr)
library(DESeq2)
library(tidyverse)

# Set wd
setwd('~/gitHub/panacea-rnaseq/cold_slurry_drg/')
dir.create('./output/DESEQ', showWarnings = F)
source('./scripts/functions.R')

# Read SIF
sif <- read.csv('./counts/sif.tsv', sep='\t', row.names = 1)
sif <- sif %>% mutate('group' = paste(type, side, time, sep = '_'))

# Read counts file
counts_df <- read.csv('./counts/metaReadCount.csv')
counts_df <- counts_df %>% column_to_rownames(var='GeneID')
gene_names <- data.frame('gene_id'=rownames(counts_df),
                         'gene_name'=counts_df$GeneName,
                         'description'=counts_df$Description)
counts_df <- counts_df[,rownames(sif)]
all_samples <- colnames(counts_df)

# get the dds object
dds = DESeqDataSetFromMatrix(as.matrix(counts_df), 
                             colData=sif, 
                             design=~group)

dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)
vsd <- assay(vst)
plot_pca('all_samples.png', vst, 'group')

# remove CS14d_2_ipsi as it is a suspected outlier
counts_df <- counts_df[, -which(colnames(counts_df) == 'CS14d_2_DRG_Ipsi')]
sif <- sif[-which(rownames(sif) == 'CS14d_2_DRG_Ipsi'), ]
dds = DESeqDataSetFromMatrix(as.matrix(counts_df), 
                             colData=sif, 
                             design=~group)

dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)
vsd <- assay(vst)
plot_pca('CS14D_2_remove.png', vst, 'group')

# CS ipsi 7D vs CTRL ipsi
CS_ipsi_7D_vs_ctrl_ipsi <- make_comparisons(dds, 
                                            c("group", "CS_ipsilateral_7D", "CTRL_ipsilateral_0H"),
                                            gene_names)
write.csv(CS_ipsi_7D_vs_ctrl_ipsi, './output/DESEQ/CS_ipsi_7D_vs_ctrl_ipsi.csv')
plot_volcano(CS_ipsi_7D_vs_ctrl_ipsi, 0.05, 2, 'CS_ipsi_7D_vs_ctrl_ipsi_padj_0.05_logfc_2')
plot_volcano(CS_ipsi_7D_vs_ctrl_ipsi, 0.05, 1.5, 'CS_ipsi_7D_vs_ctrl_ipsi_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsi_7D_vs_ctrl_ipsi, 0.05, 1, 'CS_ipsi_7D_vs_ctrl_ipsi_padj_0.05_logfc_1')

# CS ipsi 14D vs CTRL ipsi
CS_ipsi_14D_vs_ctrl_ipsi <- make_comparisons(dds, 
                                            c("group", "CS_ipsilateral_14D", "CTRL_ipsilateral_0H"),
                                            gene_names)
write.csv(CS_ipsi_14D_vs_ctrl_ipsi, './output/DESEQ/CS_ipsi_14D_vs_ctrl_ipsi.csv')
plot_volcano(CS_ipsi_14D_vs_ctrl_ipsi, 0.05, 2, 'CS_ipsi_14D_vs_ctrl_ipsi_padj_0.05_logfc_2')
plot_volcano(CS_ipsi_14D_vs_ctrl_ipsi, 0.05, 1.5, 'CS_ipsi_14D_vs_ctrl_ipsi_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsi_14D_vs_ctrl_ipsi, 0.05, 1, 'CS_ipsi_14D_vs_ctrl_ipsi_padj_0.05_logfc_1')


# RT 6H ipsi vs CS ipsi 7D
RT_ipsi_6H_vs_CS_ipsi_7 <- make_comparisons(dds, 
                                            c("group", "RT_ipsilateral_6H", "CS_ipsilateral_7D"),
                                            gene_names)
write.csv(RT_ipsi_6H_vs_CS_ipsi_7, './output/DESEQ/RT_ipsi_6H_vs_CS_ipsi_7.csv')
plot_volcano(RT_ipsi_6H_vs_CS_ipsi_7, 0.05, 2, 'RT_ipsi_6H_vs_CS_ipsi_7_padj_0.05_logfc_2')
plot_volcano(RT_ipsi_6H_vs_CS_ipsi_7, 0.05, 1.5, 'RT_ipsi_6H_vs_CS_ipsi_7_padj_0.05_logfc_1.5')
plot_volcano(RT_ipsi_6H_vs_CS_ipsi_7, 0.05, 1, 'RT_ipsi_6H_vs_CS_ipsi_7_padj_0.05_logfc_1')

# RT 6H ipsi vs CS ipsi 14D
RT_ipsi_6H_vs_CS_ipsi_14D <- make_comparisons(dds, 
                                            c("group", "RT_ipsilateral_6H", "CS_ipsilateral_14D"),
                                            gene_names)
write.csv(RT_ipsi_6H_vs_CS_ipsi_14D, './output/DESEQ/RT_ipsi_6H_vs_CS_ipsi_14D.csv')
plot_volcano(RT_ipsi_6H_vs_CS_ipsi_14D, 0.05, 2, 'RT_ipsi_6H_vs_CS_ipsi_14D_padj_0.05_logfc_2')
plot_volcano(RT_ipsi_6H_vs_CS_ipsi_14D, 0.05, 1.5, 'RT_ipsi_6H_vs_CS_ipsi_14D_padj_0.05_logfc_1.5')
plot_volcano(RT_ipsi_6H_vs_CS_ipsi_14D, 0.05, 1, 'RT_ipsi_6H_vs_CS_ipsi_14D_padj_0.05_logfc_1')

#14 day CS ipsi vs 14 day CS contra
CS_ipsilateral_14D_vs_CS_contralateral_14D <- make_comparisons(dds, 
                                              c("group", "CS_ipsilateral_14D", "CS_contralateral_14D"),
                                              gene_names)
write.csv(CS_ipsilateral_14D_vs_CS_contralateral_14D , './output/DESEQ/CS_ipsilateral_14D_vs_CS_contralateral_14D.csv')
plot_volcano(CS_ipsilateral_14D_vs_CS_contralateral_14D , 0.05, 2, 'CS_ipsilateral_14D_vs_CS_contralateral_14D_padj_0.05_logfc_2')
plot_volcano(CS_ipsilateral_14D_vs_CS_contralateral_14D , 0.05, 1.5, 'CS_ipsilateral_14D_vs_CS_contralateral_14D_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsilateral_14D_vs_CS_contralateral_14D , 0.05, 1, 'CS_ipsilateral_14D_vs_CS_contralateral_14D_padj_0.05_logfc_1')


