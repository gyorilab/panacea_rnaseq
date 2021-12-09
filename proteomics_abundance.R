library(xlsx)
library(stringr)
library(RColorBrewer)

wd <- paste0('~/gitHub/panacea-rnaseq/h_ipsc_paper/')
setwd(wd)
source('./functions.R')
# Read Human mouse ortho table
ortho <- read.csv('~/.biomart/human_mouse_ortho.tsv', sep='\t')

# Read ion channels
ion_c <- read.xlsx('./gene_lists/ion channel.xlsx', sheetIndex = 1)

# Read abundance data
# h DRG
pro_abundance <- read.xlsx('./gene_lists/Abundance score calculation (hDRG).xlsx',
                           sheetIndex = 1)

# Filter pro abundance to ion channels
pro_abundance <- pro_abundance[pro_abundance$Gene.Symbol %in% 
                                 ion_c$Symbol, ]

ion_channels <- list('Nav1.3' = 'SCN3A',
                     'Nav1.7' = 'SCN9A',
                     'Nav1.9' = 'SCN11A')

# Calculate abundance relative to each ion channel
# hdrg
for(i in 1:length(ion_channels)){
  ion = ion_channels[[i]]
  col_name <- paste0('relative_abundance_to_', names(ion_channels)[i])
  pro_abundance[, col_name] <- NA
  ion_row <- which(pro_abundance$Gene.Symbol == ion)
  for (r in 1:nrow(pro_abundance)) {
    g_abd <- pro_abundance[r, 'Abundance.score']
    pro_abundance[r, col_name] <- g_abd/pro_abundance[ion_row, 'Abundance.score']
  }
}

write.csv(pro_abundance, './output/hDRG_ion_channel_protein_abundance.csv')


# mdrg
mdrg_pro_abundance <- read.xlsx('./gene_lists/Abundance score calculation (mDRG).xlsx',
                                sheetIndex = 1)
colnames(mdrg_pro_abundance)[4] <- 'MOUSE_SYMBOL'
mdrg_pro_abundance <- merge(mdrg_pro_abundance, ortho, by='MOUSE_SYMBOL')

mdrg_pro_abundance <- mdrg_pro_abundance[mdrg_pro_abundance$HUMAN_SYMBOL %in% 
                                           ion_c$Symbol, ]
mdrg_pro_abundance$MOUSE_SYMBOL <- NULL
mdrg_pro_abundance <- mdrg_pro_abundance[, c(7, 2:ncol(mdrg_pro_abundance)-1)]
colnames(mdrg_pro_abundance)[1] <- 'Gene.name'

mdrg_pro_abundance <- get_relative_abundance(mdrg_pro_abundance)
write.csv(mdrg_pro_abundance, './output/mDRG_ion_channel_protein_abundance.csv')


# Nociceptor 4 Weeks
noc4W_pro_abundance <- read.xlsx('./gene_lists/Panacea_April2021_iPSC_proteinquant_working (hDRG 4 w).xlsx',
                                sheetIndex = 1)
colnames(noc4W_pro_abundance)[2] <- 'Gene.name'

noc4W_pro_abundance <- noc4W_pro_abundance[noc4W_pro_abundance$Gene.name %in% 
                                           ion_c$Symbol,]
noc4W_pro_abundance$Abundance.score <- rowMeans(noc4W_pro_abundance[ ,c(8,9)])
noc4W_pro_abundance <- get_relative_abundance(noc4W_pro_abundance)
write.csv(noc4W_pro_abundance, './output/noc4W_ion_channel_protein_abundance.csv')


# Nociceptor 8 weeks
noc8W_pro_abundance <- read.xlsx('./gene_lists/Panacea_April2021_iPSC_proteinquant_working (hDRG 8w).xlsx',
                                 sheetIndex = 1)
colnames(noc8W_pro_abundance)[2] <- 'Gene.name'

noc8W_pro_abundance <- noc8W_pro_abundance[noc8W_pro_abundance$Gene.name %in% 
                                             ion_c$Symbol,]
noc8W_pro_abundance$Abundance.score <- rowMeans(noc8W_pro_abundance[ ,c(8,9)])
noc8W_pro_abundance <- get_relative_abundance(noc8W_pro_abundance)
write.csv(noc8W_pro_abundance, './output/noc8W_ion_channel_protein_abundance.csv')

pro_abundance$Abundance.relative.to.Nav1.7 <- NULL
colnames(pro_abundance)[4] <- 'Gene.name'

pro_abundance_list <- list('human_DRG' = pro_abundance,
                           'mouse_DRG' = mdrg_pro_abundance,
                           'Nociceptor_4Weeks' = noc4W_pro_abundance,
                           'Nociceptor_8Weeks' = noc8W_pro_abundance
                           )

# read transcriptomics all ions
transcription_all_ions_c <- read.csv('./output/merged_df_all_ion_channels.csv', row.names = 1)

# subset to nav 1.7
transcirption_nav1.7 <- transcription_all_ions_c[, na.omit(str_extract(colnames(transcription_all_ions_c),
                                                                       '.*Nav1.7.*'))]
trans_vs_pro <- as.data.frame(matrix(ncol = 4, nrow = nrow(transcirption_nav1.7)))
colnames(trans_vs_pro) <- names(pro_abundance_list)
rownames(trans_vs_pro) <- rownames(transcirption_nav1.7)
  
for(cols in 1:ncol(transcirption_nav1.7)){
  pro_df <- pro_abundance_list[[cols]]
  pro_df <- pro_df[pro_df$Gene.name %in% rownames(transcirption_nav1.7), ]
  no_genes <- rownames(transcirption_nav1.7)[!rownames(transcirption_nav1.7) %in% pro_df$Gene.name]
  pro_df <- pro_df[!duplicated(pro_df$Gene.name), ]
  rownames(pro_df) <- pro_df$Gene.name
  pro_df[no_genes, ] <- 0
  pro_df$Gene.name <- rownames(pro_df)
  col_no <- which(colnames(pro_df) == na.omit(str_extract(colnames(pro_df), '.*1.7.*')))
  trans_vs_pro[, names(pro_abundance_list)[cols]] <- (pro_df[,col_no]+1)/(transcirption_nav1.7[, cols]+1)

}
trans_vs_pro[is.na(trans_vs_pro)] <- 0
trans_vs_pro[sapply(trans_vs_pro, is.infinite)] <- 0
x <- t(scale(t(trans_vs_pro)))
x <- na.omit(x)


breaksList <- seq(0, 10, by = 2)
colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaksList))
p <- pheatmap(trans_vs_pro,
              color = colors,
              border_color = NA,
              fontsize_row = 4, 
              cluster_cols = F,
              breaks = breaksList)

png(paste0('./output/', 'nav1.7_transcript_proteome', '.png'), height=4800,
    width=2000, res=250)
print(p)
dev.off()
