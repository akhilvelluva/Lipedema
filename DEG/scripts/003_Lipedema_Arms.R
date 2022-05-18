
# Current and Final Used for Proposal -------------------------------------
# All Samples are involved except outlier ---------------------------------
# Paired and unpaired samples using voom


# Load packages -----------------------------------------------------------
source("src/helper.R")
source("src/packages.R")

# preprocessing of paired -------------------------------------------------
paired.samples <- read.table("data/paired/feautureCounts_multimapper_NoDuplRem.summary",
                             header=TRUE, row.names=1)
paired.samples <- paired.samples[ ,6:ncol(paired.samples)]
colnames(paired.samples) <- gsub("*_starAligned.sortedByCoord.out.bam$", "", 
                                 colnames(paired.samples))
colnames(paired.samples) <- gsub("X*......$", "", colnames(paired.samples))


# preprocessing of control ------------------------------------------------

# load files from a path i.e ( Data )
# a list of dataframes will be enlisted in a variable ( dataset )

files_control <- list.files(path = "data/control",pattern = ".txt", full.names = TRUE)
files <- c(files_control)
dataset <- lapply(files, vroom, delim = "\t", col_names = FALSE) 
filenames_control <- list.files(path = "data/control",
                                pattern = ".txt") %>% gsub("*_STAR.grch38-aln_sorted.bam.txt", "", .)
filenames <- c(filenames_control)
names(dataset) <- filenames


# Prepare assay  ---------------------------------------

# Merge all dataframes enlisted in the variable
# This creates a NAMED count matrix from the dataframes
for(d in names(dataset)){
  colnames(dataset[[d]]) <- c("ENSEMBL",d)
  
}

unpaired.samples <- dataset %>% reduce(inner_join, by="ENSEMBL") %>% as.data.frame
ensembl_list <- unpaired.samples$ENSEMBL
unpaired.samples <- unpaired.samples %>% dplyr::select(-c(ENSEMBL)) 
rownames(unpaired.samples) <- ensembl_list




# Combine count data ------------------------------------------------------
combined.samples <- merge(unpaired.samples,paired.samples,by.x="row.names", 
                          by.y="row.names", all.x=F, all.y=F) %>% 
  column_to_rownames(.,"Row.names")

colnames(combined.samples) <- c("C1","C11","C15","C16","C21","C3","C9","10",
                                "11","1","12","13","14","15","16","17","18",
                                "19","20","21","2","22","23","24","25","27",
                                "28","29","30","31","3","32","33","34","35",
                                "36","37","38","39","40","4","42","43","44","45",
                                "46","47","48","5","6","7","8","9")


# Paired sample information -----------------------------------------------
paired.meta <- read.csv("data/paired//sampleInfo.csv",sep = "\t",header = TRUE,nrows = 0)
paired.meta$State <- "Lipedema"
paired.meta <- paired.meta %>% dplyr::select(Samples,ID,Tissue,Sick,State)
paired.meta <- paired.meta %>% filter(ID %in% c(16388,24396,28423))

# Unpaired sample information ---------------------------------------------
unpaired.meta <- vroom("data/sampleInfo.csv") %>% 
  filter(Observation == "CONTROL") %>% dplyr::select(-c("Observation"))
colnames(unpaired.meta) <- c("Samples","State","Tissue")
unpaired.meta$Sick <- "healthy"
unpaired.meta$ID <- paste0("C",sub(".*LipoControl_", "", unpaired.meta$Samples)) 
unpaired.meta$Samples <- unpaired.meta$ID


# combine metada ----------------------------------------------------------
combined.metadata <- plyr::rbind.fill(unpaired.meta,paired.meta)
rownames(combined.metadata) <- combined.metadata$Samples


combined.samples <- combined.samples[,rownames(combined.metadata)]

combined.metadata <- combined.metadata %>% mutate(Tissue = case_when(
  Tissue == "Arm" ~ "ARMS",
  Tissue == "Belly" ~ "BELLY",
  Tissue == "Femoral" ~ "LEGS",
  TRUE ~ Tissue
))

combined.metadata$State <- casefold(combined.metadata$State,upper = TRUE)

# Check if metadata matches count data ------------------------------------
all(rownames(combined.metadata) == colnames(combined.samples))


rownames(combined.metadata) <- factor(paste0(combined.metadata$Tissue,"_",combined.metadata$ID))
colnames(combined.samples) <- factor(paste0(combined.metadata$Tissue,"_",combined.metadata$ID))
combined.metadata$Group <- factor(paste0(combined.metadata$State,"_",combined.metadata$Tissue))




# Start from here ---------------------------------------------------------
coldata <- combined.metadata
countdata <- combined.samples


# Remove known outliers ---------------------------------------------------
coldata <- coldata %>% filter(ID != 33631)
countdata <- countdata[,rownames(coldata)]
all(rownames(coldata) == colnames(countdata))

coldata <- coldata %>% filter(Tissue == "ARMS")
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

ggsave("results/PCA_arms.png", )

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
de_arms <- NULL
de_arms <- voom(dge, design, plot=TRUE)
model_arms <- duplicateCorrelation(de_arms, design, block = coldata$ID)
# this will include the block variable (ID or pairings) as well as corr estimates
de_arms <- voom(dge,design,plot=TRUE,block=coldata$ID,
                correlation=model_arms$consensus)
model_arms <- duplicateCorrelation(de_arms, design, block = coldata$ID)


fit_arms <- lmFit(de_arms, design, block=coldata$ID, 
                  correlation=model_arms$consensus)
fit_arms <- eBayes(fit_arms, robust = TRUE)

de_results_arms <- topTable(fit_arms, coef=2, sort.by = "p", n =Inf)

tmp <- gsub("\\..*","",row.names(de_results_arms))
de_results_arms$symbol <- mapIds(org.Hs.eg.db, keys=tmp, column="SYMBOL", keytype="ENSEMBL", multiVals="first") 
de_results_arms$Entrez <- mapIds(org.Hs.eg.db, keys=tmp, column="ENTREZID", keytype="ENSEMBL", multiVals="first") 
de_results_arms <- na.omit(de_results_arms)
de_results_arms <- de_results_arms[order(de_results_arms$adj.P.Val),]
de_results_arms$Abs.logFC <- abs(de_results_arms$logFC)
de_sig_arms <- de_results_arms %>% filter(adj.P.Val <= 0.1)

write.csv(de_results_arms,"results/de_arms.csv")


# de_genes <- de_results_arms %>% dplyr::select(symbol,logFC,adj.P.Val)
# colnames(de_genes) <- c("symbol","log2FoldChange","padj")
# volcanoPlot(de_genes,0.1,0.5,"Paired and Unpaired Samples")
