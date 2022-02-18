library(dplyr)
library(DESeq2)
library(tidyverse)
library(org.Rn.eg.db)
library(clusterProfiler)

# Set wd
setwd('~/gitHub/panacea-rnaseq/cold_slurry_drg/')
dir.create('./output/DESEQ', showWarnings = F)
source('./scripts/functions.R')

# Read Human Rat ortho
human_rat <- read.csv('./gene_lists/human_rat_ortho.tsv', sep='\t')
colnames(human_rat) <- c('rat_stable_id', 'rat_gene_name', 'human_stable_id', 'human_gene_name')

# Read SIF
sif <- read.csv('./counts/sample_information.sif', sep='\t', row.names = 1)
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

# Remove CS14d_2_ipsi as it is a suspected outlier
counts_df <- counts_df[, -which(colnames(counts_df) == 'CS14d_2_DRG_Ipsi')]
sif <- sif[-which(rownames(sif) == 'CS14d_2_DRG_Ipsi'), ]
dds = DESeqDataSetFromMatrix(as.matrix(counts_df), 
                             colData=sif, 
                             design=~group)

dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)
vsd <- assay(vst)
plot_pca('CS14D_2_remove.png', vst, 'group')

fun_res <- list()


# CS ipsi vs CS contra at 3 days and 7 days
sif_mod <- subset(sif, time=='3D' | time=='7D')
counts_df_mod <- counts_df[,rownames(sif_mod)]
all_samples_mod <- colnames(counts_df_mod)

dds_mod = DESeqDataSetFromMatrix(as.matrix(counts_df_mod), 
                                 colData=sif_mod, 
                                 design=~side)

dds_mod <- DESeq(dds_mod)
vst_mod <- varianceStabilizingTransformation(dds_mod)
vsd_mod <- assay(vst_mod)

cs_3d_7d_ipsi_vs_3d_7d_contra <- make_comparisons(dds_mod, 
                                                  c("side", "ipsilateral", 
                                                    "contralateral"),
                                                  gene_names)

write.csv(cs_3d_7d_ipsi_vs_3d_7d_contra, './output/DESEQ/CS_ipsi_3D_7D_vs_CS_contra_3D_7D.csv')
diff_genes <- get_diff_genes_subset(cs_3d_7d_ipsi_vs_3d_7d_contra, 0.05, 2, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsi_3D_7D_vs_CS_contra_3D_7D.csv')

##
entrez_ids <- convert_to_entrez(rownames(diff_genes), OrgDb = org.Rn.eg.db)$ENTREZID
fun_res[['cs_3d_7d_ipsi_vs_3d_7d_contra']] <- run_functional_analysis('cs_3d_7d_ipsi_vs_3d_7d_contra', entrez_ids, 'rno', org.Rn.eg.db)

##
plot_volcano(cs_3d_7d_ipsi_vs_3d_7d_contra, 0.05, 2, 'CS_ipsi_3D_7D_vs_CS_contra_3D_7D_padj_0.05_logfc_2')


# CS ipsi 7D vs CTRL ipsi
CS_ipsi_7D_vs_ctrl_ipsi <- make_comparisons(dds, 
                                            c("group", "CS_ipsilateral_7D", "CTRL_ipsilateral_0H"),
                                            gene_names)
write.csv(CS_ipsi_7D_vs_ctrl_ipsi, './output/DESEQ/CS_ipsi_7D_vs_ctrl_ipsi.csv')
diff_genes <- get_diff_genes_subset(CS_ipsi_7D_vs_ctrl_ipsi, 0.05, 2, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsi_7D_vs_ctrl_ipsi.csv')
##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsi_7D_vs_ctrl_ipsi']] <- run_functional_analysis('CS_ipsi_7D_vs_ctrl_ipsi', entrez_ids, 'rno', org.Rn.eg.db)
##
plot_volcano(CS_ipsi_7D_vs_ctrl_ipsi, 0.05, 2, 'CS_ipsi_7D_vs_ctrl_ipsi_padj_0.05_logfc_2')
plot_volcano(CS_ipsi_7D_vs_ctrl_ipsi, 0.05, 1.5, 'CS_ipsi_7D_vs_ctrl_ipsi_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsi_7D_vs_ctrl_ipsi, 0.05, 1, 'CS_ipsi_7D_vs_ctrl_ipsi_padj_0.05_logfc_1')



# CS ipsi 14D vs CTRL ipsi
CS_ipsi_14D_vs_ctrl_ipsi <- make_comparisons(dds, 
                                            c("group", "CS_ipsilateral_14D", "CTRL_ipsilateral_0H"),
                                            gene_names)
write.csv(CS_ipsi_14D_vs_ctrl_ipsi, './output/DESEQ/CS_ipsi_14D_vs_ctrl_ipsi.csv')
diff_genes <- get_diff_genes_subset(CS_ipsi_14D_vs_ctrl_ipsi, 0.05, 2, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsi_14D_vs_ctrl_ipsi.csv')
##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsi_14D_vs_ctrl_ipsi']] <- run_functional_analysis('CS_ipsi_14D_vs_ctrl_ipsi', entrez_ids, 'rno', org.Rn.eg.db)


plot_volcano(CS_ipsi_14D_vs_ctrl_ipsi, 0.05, 2, 'CS_ipsi_14D_vs_ctrl_ipsi_padj_0.05_logfc_2')
plot_volcano(CS_ipsi_14D_vs_ctrl_ipsi, 0.05, 1.5, 'CS_ipsi_14D_vs_ctrl_ipsi_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsi_14D_vs_ctrl_ipsi, 0.05, 1, 'CS_ipsi_14D_vs_ctrl_ipsi_padj_0.05_logfc_1')


# CS ipsi 7D vs RT 6H ipsi
CS_ipsi_7D_vs_RT_ipsi_6H <- make_comparisons(dds, 
                                            c("group", "CS_ipsilateral_7D", "RT_ipsilateral_6H"),
                                            gene_names)
write.csv(CS_ipsi_7D_vs_RT_ipsi_6H, './output/DESEQ/CS_ipsi_7D_vs_RT_ipsi_6H.csv')
diff_genes <- get_diff_genes_subset(CS_ipsi_7D_vs_RT_ipsi_6H, 0.05, 1, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsi_7D_vs_RT_ipsi_6H.csv')

##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsi_7D_vs_RT_ipsi_6']] <- run_functional_analysis('CS_ipsi_7D_vs_RT_ipsi_6', entrez_ids, 'rno', org.Rn.eg.db)

plot_volcano(CS_ipsi_7D_vs_RT_ipsi_6H, 0.05, 2, 'CS_ipsi_7D_vs_RT_ipsi_6H_padj_0.05_logfc_2')
plot_volcano(CS_ipsi_7D_vs_RT_ipsi_6H, 0.05, 1.5, 'CS_ipsi_7D_vs_RT_ipsi_6H_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsi_7D_vs_RT_ipsi_6H, 0.05, 1, 'CS_ipsi_7D_vs_RT_ipsi_6H_padj_0.05_logfc_1')


# CS ipsi 14D vs RT 6H ipsi
CS_ipsi_14D_vs_RT_ipsi_6H <- make_comparisons(dds, 
                                            c("group", "CS_ipsilateral_14D", "RT_ipsilateral_6H"),
                                            gene_names)
write.csv(CS_ipsi_14D_vs_RT_ipsi_6H, './output/DESEQ/CS_ipsi_14D_vs_RT_ipsi_6H.csv')
diff_genes <- get_diff_genes_subset(CS_ipsi_14D_vs_RT_ipsi_6H, 0.05, 1, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsi_14D_vs_RT_ipsi_6H.csv')
##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsi_14D_vs_RT_ipsi_6H']] <- run_functional_analysis('CS_ipsi_14D_vs_RT_ipsi_6H', entrez_ids, 'rno', org.Rn.eg.db)

plot_volcano(CS_ipsi_14D_vs_RT_ipsi_6H, 0.05, 2, 'CS_ipsi_14D_vs_RT_ipsi_6H_padj_0.05_logfc_2')
plot_volcano(CS_ipsi_14D_vs_RT_ipsi_6H, 0.05, 1.5, 'CS_ipsi_14D_vs_RT_ipsi_6H_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsi_14D_vs_RT_ipsi_6H, 0.05, 1, 'CS_ipsi_14D_vs_RT_ipsi_6H_padj_0.05_logfc_1')


#14 day CS ipsi vs 14 day CS contra
CS_ipsilateral_14D_vs_CS_contralateral_14D <- make_comparisons(dds, 
                                              c("group", "CS_ipsilateral_14D", "CS_contralateral_14D"),
                                              gene_names)
write.csv(CS_ipsilateral_14D_vs_CS_contralateral_14D , './output/DESEQ/CS_ipsilateral_14D_vs_CS_contralateral_14D.csv')
diff_genes <- get_diff_genes_subset(CS_ipsilateral_14D_vs_CS_contralateral_14D, 0.05, 1, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsilateral_14D_vs_CS_contralateral_14D.csv')
##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsilateral_14D_vs_CS_contralateral_14D']] <- run_functional_analysis('CS_ipsilateral_14D_vs_CS_contralateral_14D', entrez_ids, 'rno', org.Rn.eg.db)

plot_volcano(CS_ipsilateral_14D_vs_CS_contralateral_14D , 0.05, 2, 'CS_ipsilateral_14D_vs_CS_contralateral_14D_padj_0.05_logfc_2')
plot_volcano(CS_ipsilateral_14D_vs_CS_contralateral_14D , 0.05, 1.5, 'CS_ipsilateral_14D_vs_CS_contralateral_14D_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsilateral_14D_vs_CS_contralateral_14D , 0.05, 1, 'CS_ipsilateral_14D_vs_CS_contralateral_14D_padj_0.05_logfc_1')


#7 day CS ipsi vs 7 day CS contra
CS_ipsilateral_7D_vs_CS_contralateral_7D <- make_comparisons(dds, 
                                                               c("group", "CS_ipsilateral_7D", "CS_contralateral_7D"),
                                                               gene_names)
write.csv(CS_ipsilateral_7D_vs_CS_contralateral_7D , './output/DESEQ/CS_ipsilateral_7D_vs_CS_contralateral_7D.csv')
diff_genes <- get_diff_genes_subset(CS_ipsilateral_7D_vs_CS_contralateral_7D, 0.05, 1, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsilateral_7D_vs_CS_contralateral_7D.csv')
##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsilateral_7D_vs_CS_contralateral_7D']] <- run_functional_analysis('CS_ipsilateral_7D_vs_CS_contralateral_7D', entrez_ids, 'rno', org.Rn.eg.db)

plot_volcano(CS_ipsilateral_7D_vs_CS_contralateral_7D , 0.05, 2, 'CS_ipsilateral_7D_vs_CS_contralateral_7D_padj_0.05_logfc_2')
plot_volcano(CS_ipsilateral_7D_vs_CS_contralateral_7D , 0.05, 1.5, 'CS_ipsilateral_7D_vs_CS_contralateral_7D_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsilateral_7D_vs_CS_contralateral_7D , 0.05, 1, 'CS_ipsilateral_7D_vs_CS_contralateral_7D_padj_0.05_logfc_1')



# 1D IPSI vs 6H RT IPSI
CS_ipsilateral_1D_vs_RT_ipsilateral_6H <- make_comparisons(dds, 
                                                           c("group", "CS_ipsilateral_1D", "RT_ipsilateral_6H"),
                                                           gene_names)
write.csv(CS_ipsilateral_1D_vs_RT_ipsilateral_6H  , './output/DESEQ/CS_ipsilateral_1D_vs_RT_ipsilateral_6H .csv')
diff_genes <- get_diff_genes_subset(CS_ipsilateral_1D_vs_RT_ipsilateral_6H, 0.05, 2, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsilateral_1D_vs_RT_ipsilateral_6H.csv')
##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsilateral_1D_vs_RT_ipsilateral_6H']] <- run_functional_analysis('CS_ipsilateral_1D_vs_RT_ipsilateral_6H', 
                                                                               entrez_ids, 
                                                                               'rno', 
                                                                               org.Rn.eg.db)

plot_volcano(CS_ipsilateral_1D_vs_RT_ipsilateral_6H , 0.05, 2, 'CS_ipsilateral_1D_vs_RT_ipsilateral_6H_padj_0.05_logfc_2')
plot_volcano(CS_ipsilateral_1D_vs_RT_ipsilateral_6H , 0.05, 1.5, 'CS_ipsilateral_1D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsilateral_1D_vs_RT_ipsilateral_6H , 0.05, 1, 'CS_ipsilateral_1D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1')


# 3D IPSI vs 6H RT IPSI
CS_ipsilateral_3D_vs_RT_ipsilateral_6H <- make_comparisons(dds, 
                                                           c("group", 
                                                             "CS_ipsilateral_3D", 
                                                             "RT_ipsilateral_6H"),
                                                           gene_names)
write.csv(CS_ipsilateral_3D_vs_RT_ipsilateral_6H  , './output/DESEQ/CS_ipsilateral_3D_vs_RT_ipsilateral_6H .csv')
diff_genes <- get_diff_genes_subset(CS_ipsilateral_3D_vs_RT_ipsilateral_6H, 0.05, 2, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsilateral_3D_vs_RT_ipsilateral_6H.csv')

##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsilateral_1D_vs_RT_ipsilateral_6H']] <- run_functional_analysis('CS_ipsilateral_1D_vs_RT_ipsilateral_6H', 
                                                                               entrez_ids, 
                                                                               'rno', 
                                                                               org.Rn.eg.db)

plot_volcano(CS_ipsilateral_3D_vs_RT_ipsilateral_6H , 0.05, 2, 'CS_ipsilateral_3D_vs_RT_ipsilateral_6H_padj_0.05_logfc_2')
plot_volcano(CS_ipsilateral_3D_vs_RT_ipsilateral_6H , 0.05, 1.5, 'CS_ipsilateral_3D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsilateral_3D_vs_RT_ipsilateral_6H , 0.05, 1, 'CS_ipsilateral_3D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1')


# 7D IPSI vs 6H RT IPSI
CS_ipsilateral_7D_vs_RT_ipsilateral_6H <- make_comparisons(dds, 
                                                           c("group", "CS_ipsilateral_7D", "RT_ipsilateral_6H"),
                                                           gene_names)
write.csv(CS_ipsilateral_7D_vs_RT_ipsilateral_6H  , './output/DESEQ/CS_ipsilateral_7D_vs_RT_ipsilateral_6H.csv')
diff_genes <- get_diff_genes_subset(CS_ipsilateral_7D_vs_RT_ipsilateral_6H, 0.05, 2, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsilateral_7D_vs_RT_ipsilateral_6H.csv')

##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsilateral_1D_vs_RT_ipsilateral_6H']] <- run_functional_analysis('CS_ipsilateral_1D_vs_RT_ipsilateral_6H', 
                                                                               entrez_ids, 
                                                                               'rno', 
                                                                               org.Rn.eg.db)

plot_volcano(CS_ipsilateral_7D_vs_RT_ipsilateral_6H , 0.05, 2, 'CS_ipsilateral_7D_vs_RT_ipsilateral_6H_padj_0.05_logfc_2')
plot_volcano(CS_ipsilateral_7D_vs_RT_ipsilateral_6H , 0.05, 1.5, 'CS_ipsilateral_7D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsilateral_7D_vs_RT_ipsilateral_6H , 0.05, 1, 'CS_ipsilateral_7D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1')



# 14D IPSI vs 6H RT IPSI
CS_ipsilateral_14D_vs_RT_ipsilateral_6H <- make_comparisons(dds, 
                                                            c("group", "CS_ipsilateral_14D", 
                                                              "RT_ipsilateral_6H"),
                                                            gene_names)
write.csv(CS_ipsilateral_14D_vs_RT_ipsilateral_6H, 
          './output/DESEQ/CS_ipsilateral_14D_vs_RT_ipsilateral_6H.csv')
diff_genes <- get_diff_genes_subset(CS_ipsilateral_14D_vs_RT_ipsilateral_6H, 0.05, 2, human_rat)
write.csv(diff_genes , './output/DESEQ/diff_CS_ipsilateral_14D_vs_RT_ipsilateral_6H.csv')

##
entrez_ids <- convert_to_entrez(rownames(diff_genes), org.Rn.eg.db)
fun_res[['CS_ipsilateral_1D_vs_RT_ipsilateral_6H']] <- run_functional_analysis('CS_ipsilateral_1D_vs_RT_ipsilateral_6H', 
                                                                               entrez_ids, 
                                                                               'rno', 
                                                                               org.Rn.eg.db)

plot_volcano(CS_ipsilateral_14D_vs_RT_ipsilateral_6H , 0.05, 2, 'CS_ipsilateral_14D_vs_RT_ipsilateral_6H_padj_0.05_logfc_2')
plot_volcano(CS_ipsilateral_14D_vs_RT_ipsilateral_6H , 0.05, 1.5, 'CS_ipsilateral_14D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1.5')
plot_volcano(CS_ipsilateral_14D_vs_RT_ipsilateral_6H , 0.05, 1, 'CS_ipsilateral_14D_vs_RT_ipsilateral_6H_padj_0.05_logfc_1')




# Make dotplots
all_samples <- list()
all_samples[['all']] <- fun_res
make_functional_analysis_dotplot(all_samples)

