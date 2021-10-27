
library(xlsx)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)



# set wd
setwd('./gitHub/panacea-rnaseq/cs_drg/')

# functions
get_colors <- function(){
  my_colors = brewer.pal(n = 10, name = "RdBu")
  my_colors = colorRampPalette(my_colors)(50)
  my_colors = rev(my_colors)
  return(my_colors)
}


draw_heatmap <- function(x, fname){
  my_colors <- get_colors()
  h <- pheatmap(x,
                #scale = "row",
                color = my_colors,
                border_color = NA,
                fontsize_row = 4)
  png(paste0('./output/', fname, '.png'), width = 1500,
      height = 1800, res = 150)
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
colnames(fpkm_clean)[2:ncol(fpkm_clean)] <- gsub(".FPKM", "", colnames(fpkm_clean))

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

