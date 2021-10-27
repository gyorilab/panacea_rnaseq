library(xlsx)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library("org.Mm.eg.db")



# set wd
setwd('./gitHub/panacea-rnaseq/cs_drg/')

# functions
get_colors <- function(){
  my_colors = brewer.pal(n = 10, name = "RdBu")
  my_colors = colorRampPalette(my_colors)(50)
  my_colors = rev(my_colors)
  return(my_colors)
}


draw_heatmap <- function(x, fname, height=1800,
                         width=1500){
  my_colors <- get_colors()
  h <- pheatmap(x,
                #scale = "row",
                color = my_colors,
                border_color = NA,
                fontsize_row = 4, 
                cluster_cols = F)
  png(paste0('./output/', fname, '.png'), width = width,
      height = height, res = 150)
    print(h)
  dev.off()
}

# Read the FPKM counts
fpkm <- read.csv('./counts/FPKM.csv')

# clean FPKM table
fpkm_clean = fpkm %>%
  select(GeneID, GeneName, ends_with("FPKM")) %>%
  as.data.frame()
rownames(fpkm_clean) <- fpkm_clean$GeneID
fpkm_clean$GeneID <-  NULL
colnames(fpkm_clean)[2:ncol(fpkm_clean)] <- gsub(".FPKM", "", colnames(fpkm_clean)[2:ncol(fpkm_clean)])
colnames(fpkm_clean)

# boxplot of the fpkm counts
boxplot(fpkm_clean, las = 2)

# boxplot of log2 fpkm counts
fpkm_clean[c(2:ncol(fpkm_clean))] <- log2(fpkm_clean[,c(2:ncol(fpkm_clean))] + 1)
boxplot(fpkm_clean, las = 2)

# read the list of input files
# common up regulated RAGs
commonupregulatedRAGs_in_DRG <- read.xlsx('./gene_lists/commonupregulatedRAGs in DRG.xlsx',
                                          sheetIndex = 1)
common_up_reg_rags <- fpkm_clean[fpkm_clean$GeneName %in% commonupregulatedRAGs_in_DRG$ID,]
common_up_reg_rags <- common_up_reg_rags[!duplicated(common_up_reg_rags$GeneName),]
rownames(common_up_reg_rags) <- common_up_reg_rags$GeneName
common_up_reg_rags$GeneName <- NULL
x <- t(scale(t(common_up_reg_rags)))
x <- na.omit(x)
draw_heatmap(x, 'common_up_rags')


# read gprotein  coupled receptors
g_protein <- read.xlsx('./gene_lists/G-protein coupled receptor.xlsx', sheetIndex = 1)

symbols <- mapIds(org.Mm.eg.db, keys = g_protein$ID, keytype = "ENSEMBL", column="SYMBOL")
g_protein$Symbol <- symbols
g_protein <- fpkm_clean[fpkm_clean$GeneName %in% g_protein$Symbol,]
g_protein <- g_protein[!duplicated(g_protein$GeneName),]
rownames(g_protein) <- g_protein$GeneName
g_protein$GeneName <- NULL
x <- t(scale(t(g_protein)))
x <- na.omit(x)
draw_heatmap(x, 'g_protein_coupled_receptor')

# read ion channel
ion_channel <- read.xlsx('./gene_lists/ion channel.xlsx', sheetIndex = 1)

symbols <- mapIds(org.Mm.eg.db, keys = ion_channel$ID, keytype = "ENSEMBL", column="SYMBOL")
ion_channel$Symbol <- symbols
ion_channel <- fpkm_clean[fpkm_clean$GeneName %in% ion_channel$Symbol,]
ion_channel <- ion_channel[!duplicated(ion_channel$GeneName),]
rownames(ion_channel) <- ion_channel$GeneName
ion_channel$GeneName <- NULL
x <- t(scale(t(ion_channel)))
x <- na.omit(x)
draw_heatmap(x, 'ion_channels')

# read kinases
kinases <- read.xlsx('./gene_lists/kinase.xlsx', sheetIndex = 1)
symbols <- mapIds(org.Mm.eg.db, keys = kinases$ID, keytype = "ENSEMBL", column="SYMBOL")
kinases$Symbol <- symbols
kinases <- fpkm_clean[fpkm_clean$GeneName %in% kinases$Symbol,]
kinases <- kinases[!duplicated(kinases$GeneName),]
rownames(kinases) <- kinases$GeneName
kinases$GeneName <- NULL
x <- t(scale(t(kinases)))
x <- na.omit(x)
draw_heatmap(x, 'kinases')

# read transcription regulators
treg <- read.xlsx('./gene_lists/transcription regulator.xlsx', sheetIndex = 1)
symbols <- mapIds(org.Mm.eg.db, keys = treg$ID, keytype = "ENSEMBL", column="SYMBOL")
treg$Symbol <- symbols
treg <- fpkm_clean[fpkm_clean$GeneName %in% treg$Symbol,]
treg <- treg[!duplicated(treg$GeneName),]
rownames(treg) <- treg$GeneName
treg$GeneName <- NULL
x <- t(scale(t(treg)))
x <- na.omit(x)
draw_heatmap(x, 'tregs')

# read transmembrane receptor
trans_rec <- read.xlsx('./gene_lists/transmembrane receptor.xlsx', sheetIndex = 1)
symbols <- mapIds(org.Mm.eg.db, keys = trans_rec$ID, keytype = "ENSEMBL", column="SYMBOL")
trans_rec$Symbol <- symbols
trans_rec <- fpkm_clean[fpkm_clean$GeneName %in% trans_rec$Symbol,]
trans_rec <- trans_rec[!duplicated(trans_rec$GeneName),]
rownames(trans_rec) <- trans_rec$GeneName
trans_rec$GeneName <- NULL
x <- t(scale(t(trans_rec)))
x <- na.omit(x)
draw_heatmap(x, 'trans_membrane_receptor')


