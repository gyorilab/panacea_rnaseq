library(ggplot2)

wd <- paste0('~/gitHub/panacea-rnaseq/h_ipsc_paper/')
setwd(wd)

# Read ION channels
ion_c <- read.csv('./output/merged_df_ion_channels.csv', row.names = 1)

hdrg_noc_8W <- ion_c %>% dplyr::select(human_DRG_avg_tpm, Nociceptor_8Weeks_avg_tpm)
hdrg_noc_8W$ratio <- (hdrg_noc_8W$human_DRG_avg_tpm+1)/(hdrg_noc_8W$Nociceptor_8Weeks_avg_tpm+1)
print(mean(hdrg_noc_8W$ratio))
print(sd(hdrg_noc_8W$ratio))
hdrg_noc_8W$similarity_score <- abs(hdrg_noc_8W$human_DRG_avg_tpm - 
                                      hdrg_noc_8W$Nociceptor_8Weeks_avg_tpm)/nrow(hdrg_noc_8W)
print(mean(hdrg_noc_8W$similarity_score))
print(sd(hdrg_noc_8W$similarity_score))

write.csv(hdrg_noc_8W, './output/hdrg_noc_8W_similarity_ratio.csv')
hdrg_noc_8W$Genes <- 'Genes'
ggplot(hdrg_noc_8W, aes(Genes, ratio)) +
  geom_violin(alpha = 0.5) + ylim(0,7)  +
  geom_dotplot(binaxis = "y",
               stackdir = 'center',
               dotsize = 0.01)
