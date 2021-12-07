library(xlsx)

wd <- paste0('~/gitHub/panacea-rnaseq/h_ipsc_paper/')
setwd(wd)

# Read ion channels
ion_c <- read.xlsx('./gene_lists/ion channel.xlsx', sheetIndex = 1)

# Read abundance data
pro_abundance <- read.xlsx('./gene_lists/Abundance score calculation (hDRG).xlsx',
                           sheetIndex = 1)

# Filter pro abundance to ion channels
pro_abundance <- pro_abundance[pro_abundance$Gene.Symbol %in% 
                                 ion_c$Symbol, ]

ion_channels <- list('Nav1.3' = 'SCN3A',
                     'Nav1.7' = 'SCN9A',
                     'Nav1.9' = 'SCN11A')

# Calculate abundance relative to each ion channel
for(i in 1:length(ion_channels)){
  ion = ion_channels[i]
  col_name <- paste0('relative_abundance_to_', names(ion_channels)[i])
  pro_abundance[, col_name] <- NA
  ion_row <- which(pro_abundance$Gene.Symbol == ion)
  for (r in 1:nrow(pro_abundance)) {
    g_abd <- pro_abundance[r, 'Abundance.score']
    pro_abundance[r, col_name] <- g_abd/pro_abundance[ion_row, 'Abundance.score']
  }
}

write.csv(pro_abundance, './output/hDRG_ion_channel_protein_abundance.csv')

