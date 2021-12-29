library(dplyr)
library(tidyr)
library(DESeq2)
library(tidyverse)

# Set wd
setwd('~/gitHub/panacea-rnaseq/cold_slurry_drg/')

# Read counts file
counts_df <- read.csv('./counts/FPKM.csv')
counts_df <- counts_df %>% column_to_rownames(var='GeneID')

View(counts_df)
