# functions
get_colors <- function(){
  my_colors = brewer.pal(n = 10, name = "RdBu")
  my_colors = colorRampPalette(my_colors)(50)
  my_colors = rev(my_colors)
  return(my_colors)
}

draw_heatmap <- function(x, fname, height=4800,
                         width=2000, res=250){
  my_colors <- get_colors()
  h <- pheatmap(x,
                #scale = "row",
                color = my_colors,
                border_color = NA,
                fontsize_row = 4, 
                cluster_cols = F)
  dir.create('./output', showWarnings = F)
  png(paste0('./output/', fname, '.png'), width = width,
      height = height, res = res)
  print(h)
  dev.off()
}

get_relative_abundance <- function(pro_abundance_df){
  ion_channels <- list('Nav1.3' = 'SCN3A',
                       'Nav1.7' = 'SCN9A',
                       'Nav1.9' = 'SCN11A')
  for(i in 1:length(ion_channels)){
    ion = ion_channels[[i]]
    col_name <- paste0('relative_abundance_to_', names(ion_channels)[i])
    pro_abundance_df[, col_name] <- NA
    if(ion %in% pro_abundance_df$Gene.name == T){
      ion_row <- which(pro_abundance_df$Gene.name == ion)
      for (r in 1:nrow(pro_abundance_df)) {
        g_abd <- pro_abundance_df[r, 'Abundance.score']
        pro_abundance_df[r, col_name] <- g_abd/pro_abundance_df[ion_row, 'Abundance.score']
      }
    }
  }
  return(pro_abundance_df)
}
