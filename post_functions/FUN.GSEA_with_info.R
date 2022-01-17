#function GSEA implements ranked GSEA algorithm for ranked gene expression profile (gene_profile)
#and gene list (geneset) and returns similarity value based on GSEA (from -1 to 1).
#This version outputs genes with significant impact on GSEA score

GSEA_with_info <- function(gene_profile_sorted, geneset){
  geneset <- intersect(geneset,names(gene_profile_sorted))
  gene_profile_list = names(gene_profile_sorted)

  #obtain indeces corresponding to positions of geneset in gene_profile
  indeces = which(gene_profile_list %in% geneset)
  
  #sort genes in geneset based on their positions in gene_profile
  indeces_sorted = sort(indeces,decreasing=FALSE)
  gene_names <- gene_profile_list[indeces_sorted]
  n = length(geneset)
  N = length(gene_profile_list)
  
  #calculate Kolmogorov statistic (ES)
  p_hit = sapply(1:N,function(x){sum(indeces_sorted<=x)/n})
  p_miss = sapply(1:N,function(x){(x-sum(indeces_sorted<=x))/(N-n)})
  #p_hit = sapply(c(indeces_sorted,indeces_sorted-1,indeces_sorted+1),function(x){sum(indeces_sorted<=x)/n})
  #p_miss = sapply(c(indeces_sorted,indeces_sorted-1,indeces_sorted+1),function(x){(x-sum(indeces_sorted<=x))/(N-n)})
  
  es = p_hit-p_miss
  similarity_value = es[which.max(abs(es))]
  #print(similarity_value)

  #index_max <- which.max(abs(es))
  #index_max <- c(indeces_sorted,indeces_sorted-1,indeces_sorted+1)[which.max(abs(es))]
  index_max <- c(1:N)[which.max(abs(es))]
  #print(plot(es[order(c(indeces_sorted,indeces_sorted-1,indeces_sorted+1),decreasing = F)]))
  #print(plot(es,type="l"))

  if (similarity_value>0){
    gene_table <- data.frame(Gene_name=gene_names,Index=indeces_sorted,
                             Significant=ifelse(indeces_sorted<=index_max,"Y","N"))
    direction="positively"
    percent_top <- sum(gene_table$Significant=="Y")/n*100
  }else{
    gene_table <- data.frame(Gene_name=gene_names,Index=indeces_sorted,
                             Significant=ifelse(indeces_sorted>=index_max,"Y","N"))
    direction="negatively"
    percent_top <- sum(gene_table$Significant=="Y")/n*100
  }
  
  return(list(Similarity=similarity_value, Direction=direction, Top=percent_top, Genes=gene_table, ES = es))
}