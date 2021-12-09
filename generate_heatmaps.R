library(xlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(org.Mm.eg.db)

wd <- paste0('~/gitHub/panacea-rnaseq/h_ipsc_paper/')
setwd(wd)
# Source functions
source(paste0(wd, '/functions.R'))

## load datasets
# Human DRG
h_drg <- read.xlsx(paste0(wd, '/data/human_mouse_gene_foldchange.xls'), 
                   sheetIndex = 1)
h_drg <- h_drg %>% dplyr::select(columns = Human_Symbol, Ave_TPM, Ave_TPM.1)
colnames(h_drg) <- c('Human_symbol', 'human_DRG_avg_tpm', 'mouse_DRG_avg_tpm')

# Nociceptors
n_4_weeks <- read.csv('./data/Nociceptor_4_weeks.csv', row.names = 1)
colnames(n_4_weeks)[1:5] <- c('Human_symbol', 'Nociceptor_4Weeks_1', 
                         'Nociceptor_4Weeks_2', 'Nociceptor_4Weeks_3', 'Nociceptor_4Weeks_avg_tpm')
n_4_weeks <- n_4_weeks %>% dplyr::select(Human_symbol, Nociceptor_4Weeks_avg_tpm)

n_8_weeks <- read.csv('./data/Nociceptor_8weeks.csv', row.names = 1)
colnames(n_8_weeks)[1:5] <- c('Human_symbol', 'Nociceptor_8Weeks_1', 
                              'Nociceptor_8Weeks_2', 'Nociceptor_8Weeks_3', 'Nociceptor_8Weeks_avg_tpm')
n_8_weeks <- n_8_weeks %>% dplyr::select(Human_symbol, Nociceptor_8Weeks_avg_tpm)


# merge nociceptors and human drgs
merged_df <- merge(h_drg, c(n_4_weeks), by='Human_symbol')
merged_df <- merge(merged_df, c(n_8_weeks), by='Human_symbol')


# Read ION channels
ion_channel <- read.xlsx(paste0(wd, '/gene_lists/ion channel.xlsx'), sheetIndex = 1)
ion_channel <- ion_channel %>% dplyr::select(ID,Symbol)

# GPCR
g_protein <- read.xlsx(paste0(wd, 'gene_lists/G-protein coupled receptor.xlsx'), sheetIndex = 1)
g_protein <- g_protein %>% dplyr::select(ID, Symbol)
colnames(g_protein) <- c('ID', 'Human_symbol')

# Kinases
kinases <- read.xlsx(paste0(wd, 'gene_lists/kinase.xlsx'), sheetIndex = 1)
kinases <- kinases %>% dplyr::select(ID, Symbol)
colnames(kinases) <- c('ID', 'Human_symbol')

# Transcription Factors
treg <- read.xlsx(paste0(wd, 'gene_lists/transcription regulator.xlsx'), sheetIndex = 1)
treg <- treg %>% dplyr::select(ID, Symbol)
colnames(treg) <- c('ID', 'Human_symbol')

# Immune related genes
immu_genes <- read.xlsx(paste0(wd, 'gene_lists/immuerelatedgenes (12-06-21).xlsx'), sheetIndex = 1)
immu_genes <- immu_genes %>% dplyr::select(ID, Symbol)
colnames(immu_genes) <- c('ID', 'Human_symbol')

# Tublin related genes
tublin <- read.xlsx(paste0(wd, 'gene_lists/Tubilin linked genes (10-26-21).xlsx'), sheetIndex = 1)
tublin <- tublin %>% dplyr::select(ID, Symbol)
colnames(tublin) <- c('ID', 'Human_symbol')

# common up-reg rags
regs <- read.xlsx(paste0(wd, 'gene_lists/commonupregulatedRAGs in DRG.xlsx'), sheetIndex = 1)
regs <- regs %>% dplyr::select(ID, Symbol)
colnames(regs) <- c('ID', 'Human_symbol')

# Drg exc
drg_exc  <- read.xlsx(paste0(wd, 'gene_lists/EXtargets2021-07-17 2.0.xlsx'), sheetIndex = 1)
drg_exc <- drg_exc %>% dplyr::select(Symbol, Entrez.Gene.Name)
colnames(drg_exc)[1] <- c('Human_symbol')

# Top 20 TF
top_tf <- read.xlsx('./gene_lists/TFs from Ginty paper.xlsx', sheetIndex = 1)
top_tf <- top_tf %>% dplyr::select('ID', 'Symbol')
colnames(top_tf) <- c('ID', 'Human_symbol')


#ipsc protein quant
ipsc_proquant <- read.xlsx('./gene_lists/Panacea_April2021_iPSC_proteinquant_working_forSam.xlsx', sheetIndex = 1)
#motor neuron
ipsc_proquant$motor_neurons_avg <- rowMeans(ipsc_proquant[, c(8,9,10)])
#cortical neurons
ipsc_proquant$cortical_neurons_avg <- rowMeans(ipsc_proquant[ ,c(11,12,13)])
#ipsc
ipsc_proquant$ipsc_avg <- rowMeans(ipsc_proquant[ ,c(14,15,16)])
# 2weeks
ipsc_proquant$nociceptor_2W_avg <- rowMeans(ipsc_proquant[ ,c(17,18,19)])
#4weeks
ipsc_proquant$nociceptor_4W_avg <- rowMeans(ipsc_proquant[ ,c(20,21)])
#8weeks
ipsc_proquant$nociceptor_8W_avg <- rowMeans(ipsc_proquant[ ,c(22,23)])

cnums <- which(colnames(ipsc_proquant) %in% na.omit(str_extract(colnames(ipsc_proquant), '.*avg')))
ipsc_proquant <- ipsc_proquant[,c(2,cnums)]
colnames(ipsc_proquant)[1] <- c('Human_symbol')


# Plot ipsc
ipsc_proquant_df <- ipsc_proquant[ipsc_proquant$Human_symbol %in% top_tf$Human_symbol, ]
ion_ipsc_proquant_df <- ipsc_proquant[ipsc_proquant$Human_symbol %in% ion_channel$Symbol, ]
ion_ipsc_proquant_df <- ion_ipsc_proquant_df[!duplicated(ion_ipsc_proquant_df$Human_symbol), ]
rownames(ion_ipsc_proquant_df) <- ion_ipsc_proquant_df$Human_symbol
ion_ipsc_proquant_df$Human_symbol <- NULL

draw_heatmap(ion_ipsc_proquant_df, fname = 'ipsc_ion_no_scale', res=200, height=1500, width = 500)
draw_heatmap(t(scale(t(ion_ipsc_proquant_df))), fname = 'ipsc_ion_scale', res=200, height=1500, width = 500)


rownames(ipsc_proquant_df) <- ipsc_proquant_df$Human_symbol
ipsc_proquant_df$Human_symbol <- NULL
draw_heatmap(ipsc_proquant_df, fname = 'ipsc_tf', res=200, height=1000, width = 500)
write.csv(merged_df_tf, './output/merged_df_transcription_factors.csv')


# Plot tf
merged_df_tf <- merged_df[merged_df$Human_symbol %in% top_tf$Human_symbol, ]
rownames(merged_df_tf) <- merged_df_tf$Human_symbol
merged_df_tf$Human_symbol <- NULL
x <- t(scale(t(merged_df_tf)))
draw_heatmap(x, fname = 'transcription_factors', res=300, height=1500, width = 1000)
write.csv(merged_df_tf, './output/merged_df_transcription_factors.csv')

# Plot DRG exc
merged_df_drg <- merged_df[merged_df$Human_symbol %in% drg_exc$Human_symbol, ]
rownames(merged_df_drg) <- merged_df_drg$Human_symbol
merged_df_drg$Human_symbol <- NULL
x <- t(scale(t(merged_df_drg)))
draw_heatmap(x, fname = 'drg_exc', res=300, height=4000)
write.csv(merged_df_drg, './output/merged_df_drg_exc.csv')


# Plot up-reg rags
merged_df_regs <- merged_df[merged_df$Human_symbol %in% regs$Human_symbol, ]
rownames(merged_df_regs) <- merged_df_regs$Human_symbol
merged_df_regs$Human_symbol <- NULL
x <- t(scale(t(merged_df_regs)))
draw_heatmap(x, fname = 'up-reg_rags', res=300, height=8000)
write.csv(merged_df_regs, './output/merged_df_regs.csv')


# Plot tubulin related genes
merged_df_tubulin <- merged_df[merged_df$Human_symbol %in% tublin$Human_symbol, ]
rownames(merged_df_tubulin) <- merged_df_tubulin$Human_symbol
merged_df_tubulin$Human_symbol <- NULL
x <- t(scale(t(merged_df_tubulin)))
draw_heatmap(x, fname = 'tubulin_linked_genes', res=300, width = 1000, height=2000)
write.csv(merged_df_tubulin, './output/merged_df_tubulin.csv')


# Plot immune related genes
merged_df_immu <- merged_df[merged_df$Human_symbol %in% immu_genes$Human_symbol, ]
rownames(merged_df_immu) <- merged_df_immu$Human_symbol
merged_df_immu$Human_symbol <- NULL
x <- t(scale(t(merged_df_immu)))
draw_heatmap(x, fname = 'immu_related_genes', res=300, width = 1000, height=2000)
write.csv(merged_df_immu, './output/merged_df_immune_related.csv')

# Plot Ion channels
merged_df_ion_c <- merged_df[merged_df$Human_symbol %in% ion_channel$Symbol, ]
rownames(merged_df_ion_c) <- merged_df_ion_c$Human_symbol
merged_df_ion_c$Human_symbol <- NULL
x <- t(scale(t(merged_df_ion_c)))
draw_heatmap(x, fname = 'ion_channels', res=300, height = 6000)
write.csv(merged_df_ion_c, './output/merged_df_ion_channels.csv')


# Plot GPCR
merged_df_gpcr <- merged_df[merged_df$Human_symbol %in% g_protein$Human_symbol, ]
rownames(merged_df_gpcr) <- merged_df_gpcr$Human_symbol
merged_df_gpcr$Human_symbol <- NULL
x <- t(scale(t(merged_df_gpcr)))
draw_heatmap(x, fname = 'gpcrs')
write.csv(merged_df_gpcr, './output/merged_df_gpcr.csv')


# Plot Kinases
merged_df_kinases <- merged_df[merged_df$Human_symbol %in% kinases$Human_symbol, ]
rownames(merged_df_kinases) <- merged_df_kinases$Human_symbol
merged_df_kinases$Human_symbol <- NULL
x <- t(scale(t(merged_df_kinases)))
draw_heatmap(x, fname = 'kinases', res = 300, height = 10000)
write.csv(merged_df_kinases, './output/merged_df_kinases.csv')


merged_df_treg <- merged_df[merged_df$Human_symbol %in% treg$Human_symbol, ]
merged_df_treg <- merged_df_treg[!duplicated(merged_df_treg$Human_symbol), ]
rownames(merged_df_treg) <- merged_df_treg$Human_symbol
merged_df_treg$Human_symbol <- NULL
x <- t(scale(t(merged_df_treg)))
draw_heatmap(x, fname = 'treg', res = 300, height = 22000)
write.csv(merged_df_treg, './output/merged_df_treg.csv')

ion_channels <- list('Nav1.3' = 'SCN3A',
                     'Nav1.7' = 'SCN9A',
                     'Nav1.9' = 'SCN11A')


ion_c_cols <- colnames(merged_df_ion_c)

for(C in 1:length(ion_c_cols)){
  cname <- str_replace(colnames(merged_df_ion_c)[C], '_avg_tpm', '')
  for (I in names(ion_channels)) {
    merged_df_ion_c[paste0(cname, '_', I)] <- NA
  }
  for(I in names(ion_channels)){
    ion_value <- merged_df_ion_c[ion_channels[[I]], 
                                 colnames(merged_df_ion_c)[C]]
    for (r in 1:nrow(merged_df_ion_c)) {
      diff = merged_df_ion_c[r, C]/ion_value
      merged_df_ion_c[r,  paste0(cname, '_', I)] <- diff
      }
  }
}
write.csv(merged_df_ion_c, './output/merged_df_all_ion_channels.csv')

# Make heatmaps for each channel
# Nav1.3
nav1.3_cols <- na.omit(str_extract(colnames(merged_df_ion_c), '.*1.3'))
nav1.3_df <- merged_df_ion_c[, nav1.3_cols]
x <- t(scale(t(nav1.3_df)))
x <- na.omit(x)
draw_heatmap(x, fname = 'nav1.3', res = 300, height = 9000)

# Nav 1.7
nav1.7_cols <- na.omit(str_extract(colnames(merged_df_ion_c), '.*1.7'))
nav1.7_df <- merged_df_ion_c[, nav1.7_cols]
x <- t(scale(t(nav1.7_df)))
x <- na.omit(x)
draw_heatmap(x, fname = 'nav1.7', res = 300, height = 9000)


# Nav 1.9
nav1.9_cols <- na.omit(str_extract(colnames(merged_df_ion_c), '.*1.9'))
nav1.9_df <- merged_df_ion_c[, nav1.9_cols]
x <- t(scale(t(nav1.9_df)))
x <- na.omit(x)
draw_heatmap(x, fname = 'nav1.9', res = 300, height = 9000)
