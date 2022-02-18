library(org.Rn.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

setwd('~/gitHub/panacea-rnaseq/cold_slurry_drg/')
source('./scripts/functions.R')

de_genes <- read.csv('./output/DESEQ/diff_CS_ipsi_3D_7D_vs_CS_contra_3D_7D.csv')
de_genes <- na.omit(de_genes$human_gene_name)
ids <- convert_to_entrez(de_genes, from = 'SYMBOL', to = 'ENTREZID', org.Hs.eg.db)
writeLines(ids$ENTREZID, './output/DESEQ/entrez_ids_CS_ipsi_3D_7D_vs_CS_contra_3D_7D.txt', sep = ',')
