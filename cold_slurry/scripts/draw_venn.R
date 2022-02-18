library(ggvenn)

setwd('./scripts/')


# read RAGS
commonupregulatedRAGs_in_DRG <- xlsx::read.xlsx('../gene_lists/commonupregulatedRAGs in DRG.xlsx',
                                                sheetIndex = 1)

# Read comparisons
CS_ipsilateral_1D_vs_RT_ipsilateral_6H <- read.csv('../output/DESEQ/diff_CS_ipsilateral_1D_vs_RT_ipsilateral_6H.csv')
CS_ipsilateral_1D_vs_RT_ipsilateral_6H <- CS_ipsilateral_1D_vs_RT_ipsilateral_6H[CS_ipsilateral_1D_vs_RT_ipsilateral_6H$rat_gene_name %in% 
                                                                                   commonupregulatedRAGs_in_DRG$ID, ]$rat_gene_name

CS_ipsilateral_3D_vs_RT_ipsilateral_6H <- read.csv('../output/DESEQ/diff_CS_ipsilateral_3D_vs_RT_ipsilateral_6H.csv')
CS_ipsilateral_3D_vs_RT_ipsilateral_6H <- CS_ipsilateral_3D_vs_RT_ipsilateral_6H[CS_ipsilateral_3D_vs_RT_ipsilateral_6H$rat_gene_name %in% 
                                                                                   commonupregulatedRAGs_in_DRG$ID, ]$rat_gene_name

CS_ipsilateral_7D_vs_RT_ipsilateral_6H <- read.csv('../output/DESEQ/diff_CS_ipsilateral_7D_vs_RT_ipsilateral_6H.csv')
CS_ipsilateral_7D_vs_RT_ipsilateral_6H <- CS_ipsilateral_7D_vs_RT_ipsilateral_6H[CS_ipsilateral_7D_vs_RT_ipsilateral_6H$rat_gene_name %in% 
                                                                                   commonupregulatedRAGs_in_DRG$ID, ]$rat_gene_name


CS_ipsilateral_14D_vs_RT_ipsilateral_6H <- read.csv('../output/DESEQ/diff_CS_ipsilateral_14D_vs_RT_ipsilateral_6H.csv')
CS_ipsilateral_14D_vs_RT_ipsilateral_6H <- CS_ipsilateral_14D_vs_RT_ipsilateral_6H[CS_ipsilateral_14D_vs_RT_ipsilateral_6H$rat_gene_name %in% 
                                                                                   commonupregulatedRAGs_in_DRG$ID, ]$rat_gene_name

x <- list(
  '1D vs 6H' = CS_ipsilateral_1D_vs_RT_ipsilateral_6H, 
  '3D vs 6H' = CS_ipsilateral_3D_vs_RT_ipsilateral_6H, 
  '7D vs 6H' = CS_ipsilateral_7D_vs_RT_ipsilateral_6H,
  '14D vs 6H' = CS_ipsilateral_14D_vs_RT_ipsilateral_6H
)

png('../output/images/rag_diff_venn.png', width = 1500,
    height = 1500, res = 200)
p <- ggvenn(
      x, 
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
      stroke_size = 0.5, set_name_size = 4
    )
print(p)
dev.off()