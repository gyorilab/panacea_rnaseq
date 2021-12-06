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



# Plot Ion channels
merged_df_ion_c <- merged_df[merged_df$Human_symbol %in% ion_channel$Symbol, ]
rownames(merged_df_ion_c) <- merged_df_ion_c$Human_symbol
merged_df_ion_c$Human_symbol <- NULL
x <- t(scale(t(merged_df_ion_c)))
draw_heatmap(x, fname = 'ion_channels', res=300, height = 6000)

# Plot GPCR
merged_df_gpcr <- merged_df[merged_df$Human_symbol %in% g_protein$Human_symbol, ]
rownames(merged_df_gpcr) <- merged_df_gpcr$Human_symbol
merged_df_gpcr$Human_symbol <- NULL
x <- t(scale(t(merged_df_gpcr)))
draw_heatmap(x, fname = 'gpcrs')

# Plot Kinases
merged_df_kinases <- merged_df[merged_df$Human_symbol %in% kinases$Human_symbol, ]
rownames(merged_df_kinases) <- merged_df_kinases$Human_symbol
merged_df_kinases$Human_symbol <- NULL
x <- t(scale(t(merged_df_kinases)))
draw_heatmap(x, fname = 'kinases', res = 300, height = 10000)


merged_df_treg <- merged_df[merged_df$Human_symbol %in% treg$Human_symbol, ]
merged_df_treg <- merged_df_treg[!duplicated(merged_df_treg$Human_symbol), ]
rownames(merged_df_treg) <- merged_df_treg$Human_symbol
merged_df_treg$Human_symbol <- NULL
x <- t(scale(t(merged_df_treg)))
draw_heatmap(x, fname = 'treg', res = 300, height = 22000)



