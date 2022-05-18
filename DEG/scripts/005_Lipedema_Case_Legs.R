
# Current and Final Used for Proposal -------------------------------------
# All Samples are involved except outlier ---------------------------------
# Paired and unpaired samples using voom
# Only Lipedema Legs where Belly is Controls ------------------------------




# Load packages -----------------------------------------------------------
source("src/helper.R")
source("src/packages.R")



# Start from here ---------------------------------------------------------
coldata <- combined.metadata
countdata <- combined.samples


# Remove known outliers ---------------------------------------------------

coldata <- coldata %>% filter(Tissue %in% c("LEGS","BELLY"),State == "LIPEDEMA")
countdata <- countdata[,rownames(coldata)]

all(colnames(countdata) == rownames(coldata))

# Mutate data -------------------------------------------------------------

coldata <- coldata %>% mutate(State =
                                case_when(
                                  Tissue == "BELLY" ~ "CONTROLS",
                                  TRUE ~ State
                                ))

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

ggsave("results/PCA_case_legs.png", )

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
de_case_legs <- NULL
de_case_legs <- voom(dge, design, plot=TRUE)
model_case_legs <- duplicateCorrelation(de_case_legs, design, block = coldata$ID)
# this will include the block variable (ID or pairings) as well as corr estimates
de_case_legs <- voom(dge,design,plot=TRUE,block=coldata$ID,
                     correlation=model_case_legs$consensus)
model_case_legs <- duplicateCorrelation(de_case_legs, design, block = coldata$ID)


fit_case_legs <- lmFit(de_case_legs, design, block=coldata$ID, 
                       correlation=model_case_legs$consensus)
fit_case_legs <- eBayes(fit_case_legs, robust = TRUE)

de_results_case_legs <- topTable(fit_case_legs, coef=2, sort.by = "p", n =Inf)

tmp <- gsub("\\..*","",row.names(de_results_case_legs))
de_results_case_legs$symbol <- mapIds(org.Hs.eg.db, keys=tmp, column="SYMBOL", keytype="ENSEMBL", multiVals="first") 
de_results_case_legs$Entrez <- mapIds(org.Hs.eg.db, keys=tmp, column="ENTREZID", keytype="ENSEMBL", multiVals="first") 
de_results_case_legs <- na.omit(de_results_case_legs)
de_results_case_legs <- de_results_case_legs[order(de_results_case_legs$adj.P.Val),]
de_results_case_legs$Abs.logFC <- abs(de_results_case_legs$logFC)
de_sig_case_legs <- de_results_case_legs %>% filter(adj.P.Val <= 0.1)

write.csv(de_results_case_legs,"results/de_case_legs.csv")


# de_genes <- de_results_case_legs %>% dplyr::select(symbol,logFC,adj.P.Val)
# colnames(de_genes) <- c("symbol","log2FoldChange","padj")
# volcanoPlot(de_genes,0.1,0.5,"Paired and Unpaired Samples")
