
# Intersection ------------------------------------------------------------

data <- NULL
data[["LEGS_LEGS"]] <- de_sig_legs
data[["ARMS_ARMS"]] <- de_sig_arms
data[["ABD_LEGS"]] <- de_sig_case_legs
data[["ABD_ARMS"]] <- de_sig_case_arms



openxlsx::write.xlsx(data, "results/significant_genes.xlsx")


gene_list <- NULL
for (genes in names(data)) {
  print(genes)
  gene_list[[genes]] <- data[[genes]] %>% select(symbol) %>% unlist %>% unique
  
}

df <- data.frame(gene1 = "NULL",gene2 = "NULL", count = 0)

for(gene in names(gene_list)){
  print(gene)
  for(other in names(gene_list)){
    gene_count <- intersect(gene_list[[gene]],gene_list[[other]]) %>% length
    tmp <- data.frame(gene1 = gene,gene2 = other,count =gene_count)
    df <- rbind(df,tmp)
  }
}

df <- df %>% filter(gene1 != 'NULL')
g=graph.data.frame(df,directed=FALSE)
g <- get.adjacency(g,attr='count',spars=FALSE)

write.csv(g,"results/Lipedema_Audit.csv")


venn_data = list(
  "LEGS-LEGS" = gene_list$LEGS_LEGS,
  "ARMS" = gene_list$ARMS_ARMS,
  BL = gene_list$ABD_LEGS,
  BA = gene_list$ABD_ARMS
)


ggvenn(
  venn_data, 
  #fill_color = c("indianred", "darkturqoise","yellow"),
  stroke_size = 0.5, set_name_size = 4
)

intersect(de_sig_case_legs$symbol,de_sig_case_arms$symbol) %>% 
  intersect(.,de_sig_arms$symbol) %>%
  intersect(.,de_sig_legs$symbol) %>% 
  intersect(.,de_sig_subset$symbol)

ggsave("results/intersect_genes.png")
