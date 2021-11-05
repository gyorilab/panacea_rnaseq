library(xlsx)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

# set wd
setwd('~/gitHub/panacea-rnaseq/woolf_sensory_ampliseq/')

# functions
get_colors <- function(){
  my_colors = brewer.pal(n = 10, name = "RdBu")
  my_colors = colorRampPalette(my_colors)(50)
  my_colors = rev(my_colors)
  return(my_colors)
}

draw_heatmap <- function(x, fname, height=2800,
                         width=2500, res=250, anno = NULL){
  my_colors <- get_colors()
  h <- pheatmap(x,
                #scale = "row",
                color = my_colors,
                border_color = NA,
                fontsize_row = 4, 
                cluster_cols = F,
                annotation_row = anno,
                cluster_rows = F)
  png(paste0('./output/', fname, '.png'), width = width,
      height = height, res = res)
  print(h)
  dev.off()
}

# Read XLSX file
ampli_df <- read.csv('./data/sensory_ampliseq.csv')
ampli_df <- ampli_df[, c(1,3,4,5,6)]
ampli_df <- ampli_df %>% column_to_rownames('Gene')

# Read gene list
#gene_list <- read.csv('./data/gene_list.csv')
gene_list <- read.xlsx('./data/Gene_list.xlsx', sheetIndex = 3)
  
png('./output/boxplot_raw.png', width=1500, height = 1500, res=150)
print(boxplot(ampli_df))
dev.off()

ampli_log2 = log2(ampli_df + 1)
png('./output/boxplot_log_2.png', width=1500, height = 1500, res=150)
print(boxplot(ampli_log2))
dev.off()


genes_sub <- ampli_log2[rownames(ampli_log2) %in% gene_list$Gene,]
genes_sub <- genes_sub[, c(4,2,3,1)]
x <- t(scale(t(genes_sub)))
x <- na.omit(x)
# Create annotation col
anno <- gene_list[gene_list$Gene %in% rownames(x), ]
rownames(anno) <- anno$Gene
anno$Gene <- NULL
x <- x[rownames(anno),]
draw_heatmap(x, 'all_samples', anno = anno, res=300, height = 1200, width=1500)

x <- t(scale(t(genes_sub[, c(1,2,3)])))
x <- na.omit(x)
x <- x[rownames(anno),]
draw_heatmap(x, 'no_X5i.28d', anno = anno, res=300, height = 1200, width=1500)

