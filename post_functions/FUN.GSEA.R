#function GSEA implements ranked GSEA algorithm for ranked gene expression profile (gene_profile)
#and gene list (geneset) and returns similarity value based on GSEA (from -1 to 1).

GSEA <- function(gene_profile_sorted, geneset){
        gene_profile_list = names(gene_profile_sorted)
        
        #obtain indices corresponding to positions of geneset(part of set) in gene_profile(full set)
        indeces = which(gene_profile_list %in% geneset)

        #sort genes in geneset based on their positions in gene_profile
        indeces_sorted = sort(indeces, decreasing=FALSE)
        n = length(geneset)
        N = length(gene_profile_list)
                
        #calculate Kolmogorov statistic (ES)
        p_hit = sapply(c(indeces_sorted, indeces_sorted-1, indeces_sorted+1), function(x){sum(indeces_sorted<=x)/n})
        p_miss = sapply(c(indeces_sorted, indeces_sorted-1, indeces_sorted+1), function(x){(x-sum(indeces_sorted<=x))/(N-n)})
        es = p_hit - p_miss
        similarity_value = es[which.max(abs(es))]

        return(similarity_value)
}