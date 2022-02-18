library(EnhancedVolcano)


draw_heatmap <- function(x, fname, height=2800,
                         width=2500, res=250, annot_df=F){
  my_colors <- get_colors()
  h <- pheatmap(x,
                annotation_col = annot_df,
                #scale = "row",
                color = my_colors,
                border_color = NA,
                fontsize_row = 6, 
                cluster_cols = F, cluster_rows = T, show_colnames = F,
  )
  png(paste0('./output/', fname, '.png'), width = width,
      height = height, res = res)
  print(h)
  dev.off()
}


# Generate PCA
plot_pca <- function(file, vsd, Intgroup){
  dir.create('./output/images/pca', showWarnings = F, recursive = T)
  pca_data <- plotPCA(vsd, Intgroup, returnData=T)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  png(paste0('./output/images/pca/', file), width = 780*4, height = 720*4, res=300)
  p <- ggplot(pca_data, aes(PC1, PC2, color=pca_data[, Intgroup], label=name)) + 
    geom_point(size=3) +
    geom_text(aes(label=name), hjust=0, vjust=1.5) +
    xlab(paste0("PC1: ", percent_var[1], "%variance")) +
    ylab(paste0("PC2: ", percent_var[2], "%variance")) +
    coord_fixed()
  print(p)
  dev.off()
}

# Generate heatmaps
pheatmap_genarate<- function(file, vsd){
  dir.create('images/heatmaps', showWarnings = F)
  png(file, width = 780*4, height = 720*4, res=300)
  dists_count <- cor(vsd)
  print(pheatmap(dists_count,legend = TRUE, main = "Clustering heatmap for normalized samples"))
  dev.off()
}


# make comparisons
make_comparisons <- function(dds, contrast, gene_names){
  res <- results(dds, contrast = contrast)
  res$gene_id <- rownames(res)
  res <- merge(as.data.frame(res), gene_names, by.x='gene_id')
  res <- res %>% column_to_rownames(var='gene_id')
  res <- res[, c(7,8, 1:(ncol(res)-2))]
  return(res)
}


# generate volcano plots
plot_volcano <- function(res, padj_cutoff, logFC, titlename){
  
  keyvals <- ifelse(
    res$log2FoldChange < -logFC & res$padj < padj_cutoff, 'blue3',
    ifelse(res$log2FoldChange > logFC & res$padj < padj_cutoff, 'red3',
           'black'))
  up_reg <- nrow(subset(res, padj < padj_cutoff & log2FoldChange > logFC))
  down_reg <- nrow(subset(res, padj < padj_cutoff & log2FoldChange < -logFC))
  ns <- nrow(subset(res, padj > padj_cutoff))
  
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red3'] <- paste0('Upregulated ', up_reg)
  names(keyvals)[keyvals == 'black'] <- paste0('NS ', nrow(res)-(up_reg+down_reg)) 
  names(keyvals)[keyvals == 'blue3'] <- paste0('Downregulated ', down_reg)
  dir.create('./output/images/volcano/', recursive = T, showWarnings = F)
  pngname <- paste("./output/images/volcano/volcano_",titlename,".png",sep="")
  
  #top_10 <- c(up_reg, down_reg)
  #imp_genes <- c('Atf3', 'Sox11', 'Sprr1a', 'Gal', 'Gpr151')
  #rownames(res)[which(res$gene_name %in% imp_genes)] <- imp_genes
  png(pngname, width = 980*4, height = 650*4, res=300)
  print(EnhancedVolcano(res,
                        lab = res$gene_name,
                        x = 'log2FoldChange',
                        y = 'padj',
                        selectLab = c('Atf3', 'Sox11', 'Sprr1a', 'Gal', 'Gpr151'),
                        ylab = bquote(~-Log[10]~Adjusted~italic(P)),
                        FCcutoff = FALSE,
                        pCutoff = FALSE,
                        colCustom = keyvals,
                        drawConnectors = T,
                        widthConnectors = 0.3,
                        colAlpha = 4/5,
                        title = titlename,
                        pointSize = 3.0,
                        labSize = 5.0,
                        xlim = c(-10, 10),
                        colConnectors = 'black',
                        #boxedLabels = T,
                        typeConnectors = 'closed', ylim = c(-4, 130) 
                        
                        ))
  dev.off() 
}


get_diff_genes_subset <- function(res, padj_cutoff, logFC, annot){
  # map rat genes to human
  diff_genes <- subset(res, padj<padj_cutoff & (log2FoldChange >= logFC | log2FoldChange <= -logFC))
  annot <- annot[, -c(2,3)]
  diff_genes$rat_stable_id <- rownames(diff_genes)
  diff_genes <- diff_genes[, c(ncol(diff_genes), 1:ncol(diff_genes)-1)]
  diff_merged <- merge(diff_genes, annot, by='rat_stable_id', all.x=T)
  diff_merged <- diff_merged[!duplicated(diff_merged$rat_stable_id),]
  diff_merged <- diff_merged[, c(ncol(diff_merged), 1:ncol(diff_merged)-1)]
  rownames(diff_merged) <- NULL
  diff_merged <- diff_merged %>% column_to_rownames(var='rat_stable_id')
  colnames(diff_merged)[2] <- 'rat_gene_name'
  return(diff_merged)
}


convert_to_entrez <- function(ensembl, from='ENSEMBL', 
                              to='ENTREZID', OrgDb){
  DE_gene = bitr(ensembl, fromType=from,
                 toType=to, OrgDb=OrgDb)
  return(DE_gene)
}


run_functional_analysis <- function(comp_name, entrez_ids, organism, OrgDb){
  dir.create('./output/functional_analysis/KEGG', showWarnings = F, recursive = T)
  dir.create('./output/functional_analysis/GO_CC', showWarnings = F, recursive = T)
  dir.create('./output/functional_analysis/GO_BP', showWarnings = F, recursive = T)
  dir.create('./output/functional_analysis/GO_MF', showWarnings = F, recursive = T)
  
  # create files
  keggFile <- paste0('./output/functional_analysis/KEGG/', comp_name, '.csv')
  ccFile <- paste0('./output/functional_analysis/GO_CC/', comp_name, '.csv')
  bpFile <- paste0('./output/functional_analysis/GO_BP/', comp_name, '.csv')
  mfFile <- paste0('./output/functional_analysis/GO_MF/', comp_name, '.csv')
  
  kegg_result <- enrichKEGG(gene = entrez_ids, organism = organism, pvalueCutoff = 0.05)
  go_cc_all <- enrichGO(gene = entrez_ids, OrgDb = OrgDb, ont = "CC", pvalueCutoff = 0.05)
  go_bp_all <- enrichGO(gene = entrez_ids, OrgDb = OrgDb, ont = "BP", pvalueCutoff = 0.05)
  go_mf_all <- enrichGO(gene = entrez_ids, OrgDb = OrgDb, ont = "MF", pvalueCutoff = 0.05)
  
  # write csvs
  write.csv(as.data.frame(kegg_result), file=paste(keggFile, sep="")) 
  write.csv(as.data.frame(go_cc_all), file=paste(ccFile, sep="")) 
  write.csv(as.data.frame(go_bp_all), file=paste(bpFile, sep="")) 
  write.csv(as.data.frame(go_mf_all), file=paste(mfFile, sep=""))
  
  # barplots
  dir.create('./output/images/functional_analysis/KEGG', showWarnings = F, recursive = T)
  dir.create('./output/images/functional_analysis/GO_CC', showWarnings = F, recursive = T)
  dir.create('./output/images/functional_analysis/GO_BP', showWarnings = F, recursive = T)
  dir.create('./output/images/functional_analysis/GO_MF', showWarnings = F, recursive = T)
  image_outfile <- './output/images/functional_analysis/'
  
  if (length(kegg_result$Description)!=0) {
    max_desc <- max(sapply(kegg_result$Description,nchar))
    KEGG_PlotFile <- paste0(image_outfile, 'KEGG/', 'kegg_', comp_name,'.png')
    if (max_desc <50) {
      png(KEGG_PlotFile, width = 480*4, height = 420*4, res=300)
    } else if(max_desc > 50 & max_desc <= 100){
      png(KEGG_PlotFile, width = 520*4, height = 420*4, res=300)
    } else if(max_desc > 100){
      png(KEGG_PlotFile, width = 700*4, height = 420*4, res=300)
    }
    print(barplot(kegg_result, showCategory=20,title = "KEGG"))
    dev.off()
  }
  
  if (length(go_cc_all$Description)!=0) {
    max_desc <- max(sapply(go_cc_all$Description, nchar))
    go_cc_PlotFile <- paste0(image_outfile, 'GO_CC/', 'GO_CC_', comp_name,'.png')
    if (max_desc <50) {
      png(go_cc_PlotFile, width = 480*4, height = 420*4, res=300)
    } else if(max_desc > 50 & max_desc <= 100){
      png(go_cc_PlotFile, width = 520*4, height = 420*4, res=300)
    } else if(max_desc > 100){
      png(go_cc_PlotFile, width = 700*4, height = 420*4, res=300)
    }
    print(barplot(go_cc_all, showCategory=20,title = "GO_CC"))
    dev.off()
  }
  
  if (length(go_bp_all$Description)!=0) {
    max_desc <- max(sapply(go_bp_all$Description, nchar))
    go_bp_PlotFile <- paste0(image_outfile, 'GO_BP/', 'GO_BP_', comp_name,'.png')
    if (max_desc <50) {
      png(go_bp_PlotFile, width = 480*4, height = 420*4, res=300)
    } else if(max_desc > 50 & max_desc <= 100){
      png(go_bp_PlotFile, width = 520*4, height = 420*4, res=300)
    } else if(max_desc > 100){
      png(go_bp_PlotFile, width = 1000*4, height = 420*4, res=300)
    }
    print(barplot(go_bp_all, showCategory=20,title = "GO_BP"))
    dev.off()
  }
  
  
  if (length(go_mf_all$Description)!=0) {
    max_desc <- max(sapply(go_mf_all$Description, nchar))
    go_mf_PlotFile <- paste0(image_outfile, 'GO_MF/', 'GO_MF_', comp_name,'.png')
    if (max_desc <50) {
      png(go_mf_PlotFile, width = 480*4, height = 420*4, res=300)
    } else if(max_desc > 50 & max_desc <= 100){
      png(go_mf_PlotFile, width = 520*4, height = 420*4, res=300)
    } else if(max_desc > 100){
      png(go_mf_PlotFile, width = 700*4, height = 420*4, res=300)
    }
    print(barplot(go_mf_all, showCategory=20,title = "GO_mf"))
    dev.off()
  }
  
  return(list(kegg_result, go_cc_all, go_bp_all, go_mf_all))
  
}


make_functional_analysis_dotplot <- function(organized_enrichments){
  index=1
  for (comparison_set in organized_enrichments) {
    comparison_set_name <- 'all_samples'
    names_comparison_set <- names(comparison_set)
    transpose_comparison_set <- transpose(comparison_set)
    
    keggs <- transpose_comparison_set[[1]]
    kegg_comp = merge_result(keggs)
    if(any(kegg_comp@compareClusterResult$pvalue < 0.05)){
      kegg_comp@compareClusterResult <- kegg_comp@compareClusterResult %>% 
        group_by(Cluster)  %>% dplyr::slice(1:20)
      
      
      kegg_comp_file <- paste("./output/images/functional_analysis/KEGG/kegs_comparison_all",".png",sep="")
      
      png(kegg_comp_file, width = 2500, height = 3500, res = 300)

      p_kegg_comp = dotplot(kegg_comp, showCategory=20)
      print(p_kegg_comp + theme(axis.text.x = element_text (angle=90, hjust=1)))
      dev.off()
    }
    
    
    go_ccs <- transpose_comparison_set[[2]]
    go_cc_comp = merge_result(go_ccs)
    if(any(go_cc_comp@compareClusterResult$pvalue < 0.05)){
      go_cc_comp@compareClusterResult <- go_cc_comp@compareClusterResult %>% 
        group_by(Cluster)  %>% dplyr::slice(1:20)
      
      go_cc_comp_file <- paste("./output/images/functional_analysis/GO_CC/go_cc_comparison_all",".png",sep="")
      png(go_cc_comp_file, width = 2300, height = 3500, res = 300)
      p_go_cc_comp = dotplot(go_cc_comp, showCategory=20)
      print(p_go_cc_comp + theme(axis.text.x = element_text (angle=90, hjust=1)))
      dev.off()
    }
    
    go_bps <- transpose_comparison_set[[3]]
    go_bp_comp = merge_result(go_bps)
    if(any(go_bp_comp@compareClusterResult$pvalue < 0.05)){
      go_bp_comp@compareClusterResult <- go_bp_comp@compareClusterResult %>% 
        group_by(Cluster)  %>% dplyr::slice(1:20)
      
      
      go_bp_comp_file <- paste("./output/images/functional_analysis/GO_BP/go_bp_comparison_all",".png",sep="")
      png(go_bp_comp_file, width = 4200, height = 5500, res = 300)
      p_go_bp_comp = dotplot(go_bp_comp,showCategory=20)
      print(p_go_bp_comp + theme(axis.text.x = element_text (angle=90, hjust=1)))
      dev.off()
    }
    
    
    go_mfs <- transpose_comparison_set[[4]]
    go_mf_comp = merge_result(go_mfs)
    
    if(any(go_bp_comp@compareClusterResult$pvalue < 0.05)){
      go_mf_comp@compareClusterResult <- go_mf_comp@compareClusterResult %>% 
        group_by(Cluster)  %>% dplyr::slice(1:20)
      
      

      go_mf_comp_file <- paste("./output/images/functional_analysis/GO_MF/go_mf_comparison_all",".png",sep="")
      png(go_mf_comp_file, width = 3500, height = 4200, res = 300)
      p_go_mf_comp = dotplot(go_mf_comp, showCategory=20)
      print(p_go_mf_comp + theme(axis.text.x = element_text (angle=90, hjust=1)))
      dev.off()
    }
    
    index=index+1
    
  }
  
}
