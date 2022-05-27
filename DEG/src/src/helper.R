
volcanoPlot = function(test = dataset,p.threshold = 0.05, logfc.threshold = 1.5,plot.name = "Differential Expression", gene.list = NULL){
  
  vp <- ggplot(test) + 
    geom_point(aes(x = log2FoldChange, y = -log10(padj), color = (-log10(padj) >= -log10(p.threshold) & abs(log2FoldChange) >= logfc.threshold )), size = 1) + 
    geom_text_repel(data = test %>% 
                      filter(symbol %in% all_of(gene.list)), 
                    aes(label = symbol, x = log2FoldChange, y = -log10(padj)), box.padding = unit(.7, "lines"),hjust= 0.30) +
    geom_vline(xintercept = c(logfc.threshold,-1*logfc.threshold), linetype = 'dashed', color = 'red') +
    geom_hline(yintercept = c(-log10(p.threshold)), linetype = 'dashed', color = 'red') +
    ggtitle(plot.name) +
    xlab("log2 Fold Change") + 
    ylab("-log10 PVal") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))
  
  return(vp)
  
}




enrichGO <-  function(gene_vector,category = "BP",thr = 1){
  gene_set <- as.data.frame(gene_vector) %>% rownames_to_column("SYMBOL")
  colnames(gene_set) <- c("SYMBOL","PADJ")
  gene_names <- bitr(names(gene_vector), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db) 
  gene_set <- inner_join(gene_set,gene_names,by="SYMBOL") %>% dplyr::select(c(ENTREZID,PADJ)) 
  genes <- gene_set$ENTREZID
  gene_list <-gene_set$PADJ
  names(gene_list) <- genes
  enrichment <- gseGO(geneList = sort(gene_list, decreasing = TRUE),
                      OrgDb        = org.Hs.eg.db,
                      ont          = category,
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      pvalueCutoff = thr,
                      pAdjustMethod = "bonferroni",
                      verbose      = FALSE)
}

enrichKP <-  function(gene_vector,thr = 1){
  gene_set <- as.data.frame(gene_vector) %>% rownames_to_column("SYMBOL")
  colnames(gene_set) <- c("SYMBOL","PADJ")
  gene_names <- bitr(names(gene_vector), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db) 
  gene_set <- inner_join(gene_set,gene_names,by="SYMBOL") %>% dplyr::select(c(ENTREZID,PADJ)) 
  genes <- gene_set$ENTREZID
  gene_list <-gene_set$PADJ
  names(gene_list) <- genes
  enrichment <- enrichKEGG(gene         = genes,
                           organism     = 'hsa',
                           pvalueCutoff = thr, 
                           pAdjustMethod = "bonferroni")
}


enrichGS <-  function(gene_vector,thr = 1){
  gene_set <- as.data.frame(gene_vector) %>% rownames_to_column("SYMBOL")
  colnames(gene_set) <- c("SYMBOL","PADJ")
  gene_names <- bitr(names(gene_vector), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Mm.eg.db) 
  gene_set <- inner_join(gene_set,gene_names,by="SYMBOL") %>% select(c(ENTREZID,PADJ)) 
  genes <- gene_set$ENTREZID
  gene_list <-gene_set$PADJ
  names(gene_list) <- genes
  enrichment <- gseKEGG(geneList     = sort(gene_list, decreasing = TRUE),
                        organism     = 'mmu',
                        minGSSize    = 120,
                        pvalueCutoff = thr,
                        nPermSimple = 10000,
                        verbose      = FALSE)
}
