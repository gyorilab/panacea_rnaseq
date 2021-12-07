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
