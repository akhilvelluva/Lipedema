
# Current and Final Used for Proposal -------------------------------------
# All Samples are involved except outlier ---------------------------------
# Paired and unpaired samples using voom


# Load packages -----------------------------------------------------------
source("src/helper.R")
source("src/packages.R")

# Start from here ---------------------------------------------------------
coldata <- combined.metadata
countdata <- combined.samples


# Remove known outliers ---------------------------------------------------
coldata <- coldata %>% filter(ID != 33631)
countdata <- countdata[,rownames(coldata)]
all(rownames(coldata) == colnames(countdata))

coldata <- coldata %>% filter(Tissue == "BELLY")
countdata <- countdata[,rownames(coldata)]

all(colnames(countdata) == rownames(coldata))

# Mutate data -------------------------------------------------------------

# coldata <- coldata %>% mutate(State = 
#                                 case_when(
#                                   Tissue == "BELLY" ~ "CONTROLS",
#                                   TRUE ~ State
#                                 ))

coldata <- coldata %>% mutate(Tissue = factor(Tissue),
                              State = factor(State))



# Analysis ----------------------------------------------------------------

# Visulaiztion and PCA ----------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(countdata),
                              colData = coldata,
                              design = ~ 1)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# Used a different design, I consider both the state and tissue as source of variation
# with levels from state -> tissue
vsd <- vst(dds, blind=FALSE) 
batch1 <- factor(vsd$State)
batch2 <- factor(vsd$Tissue)
design <- model.matrix(~batch1+batch2)


assay(vsd) <- removeBatchEffect(assay(vsd),
                                batch = vsd$ID,
                                design = design)

pcaData <- plotPCA(vsd, intgroup=c("Tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
res.pca <- PCA(pcaData[,c(1:2)], graph = FALSE)

fviz_pca_biplot(res.pca,obs.scale = 1,label = "none",
                var.scale = 1, alpha=0, col.var="steelblue",
                title = "Paired and Unpaired Samples") + 
  geom_point(aes(shape = factor(vsd$State), colour = factor(vsd$Tissue))) +
  guides(shape = guide_legend(title = "Samples"),
         colour = guide_legend(title = "Depots")) +
  labs(x = paste0("PC1: ",percentVar[1],"% variance"), 
       y = paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel( aes( label=vsd$ID),direction="both",size=3) +
  theme (plot.title = element_text (hjust = 0.5))

ggsave("results/PCA_belly.png", )

# Differential expression -------------------------------------------------
# Analysis ----------------------------------------------------------------
library(edgeR)
library(limma)
library(statmod)


coldata <- coldata %>% dplyr::select(Samples,ID,Tissue,State)


dge <- DGEList(counts=countdata, samples = coldata, group = coldata$State)
keep<- rowSums(cpm(dge)>1) >= 6
dge <- dge[keep,]
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
dge$counts

design <- model.matrix(~coldata$State)

# Model fitting -----------------------------------------------------------

# x2 iteration
# this will estimate the weight of observed levels/factors
de_belly <- NULL
de_belly <- voom(dge, design, plot=TRUE)
model_belly <- duplicateCorrelation(de_belly, design, block = coldata$ID)
# this will include the block variable (ID or pairings) as well as corr estimates
de_belly <- voom(dge,design,plot=TRUE,block=coldata$ID,
                correlation=model_belly$consensus)
model_belly <- duplicateCorrelation(de_belly, design, block = coldata$ID)


fit_belly <- lmFit(de_belly, design, block=coldata$ID, 
                  correlation=model_belly$consensus)
fit_belly <- eBayes(fit_belly, robust = TRUE)

de_results_belly <- topTable(fit_belly, coef=2, sort.by = "p", n =Inf)

tmp <- gsub("\\..*","",row.names(de_results_belly))
de_results_belly$symbol <- mapIds(org.Hs.eg.db, keys=tmp, column="SYMBOL", keytype="ENSEMBL", multiVals="first") 
de_results_belly$Entrez <- mapIds(org.Hs.eg.db, keys=tmp, column="ENTREZID", keytype="ENSEMBL", multiVals="first") 
de_results_belly <- na.omit(de_results_belly)
de_results_belly <- de_results_belly[order(de_results_belly$adj.P.Val),]
de_results_belly$Abs.logFC <- abs(de_results_belly$logFC)
de_sig_belly <- de_results_belly %>% filter(adj.P.Val <= 0.1)

write.csv(de_results_belly,"results/de_belly.csv")


# de_genes <- de_results_belly %>% dplyr::select(symbol,logFC,adj.P.Val)
# colnames(de_genes) <- c("symbol","log2FoldChange","padj")
# volcanoPlot(de_genes,0.1,0.5,"Paired and Unpaired Samples")


gene_vector <- sig_genes_new %>% 
  select(symbol, logFC)
gene_names <- gene_vector$symbol
gene_vector <- gene_vector$logFC
names(gene_vector) <- gene_names

bp_lipedema <- enrichGO(gene_vector)
kegg_lipedema <- enrichKP(gene_vector)
bp_lipedema <- setReadable(bp_lipedema,org.Hs.eg.db,"ENTREZID")
kegg_lipedema <- setReadable(kegg_lipedema,org.Hs.eg.db,"ENTREZID")

kegg_df <- kegg_lipedema@result
bp_df <- bp_lipedema@result

msigdb <- NULL
msigdb$HALL <- gmtPathways("data/h.all.v7.5.1.symbols.gmt")
msigdb$BP <- gmtPathways("data/c5.go.bp.v7.5.1.symbols.gmt")


hallmark_lipedema <- fgsea(pathways=msigdb$HALL, 
                           stats=sort(gene_vector, decreasing = T), nPermSimple = 10000)
msigbp_lipedema <- fgsea(pathways=msigdb$BP, 
                         stats=sort(gene_vector, decreasing = T), nPermSimple = 10000)


