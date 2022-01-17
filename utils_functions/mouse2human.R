### create mouse 2 human entrez ID mapping

mouse2human <- function(human_ids){
    mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

    # human / mouse
    human2mouse_df <-getLDS( attributes=c("entrezgene_id"), 
                            attributesL=c("entrezgene_id"),
                            filters="entrezgene_id", values=human_ids, 
                            mart=mart1, martL=mart2
                            )

    #one to one mapping filter
    hcol <- human2mouse_df[, 1]
    mcol <- human2mouse_df[, 2]
    ht <- as.data.frame(table(hcol))
    mt <- as.data.frame(table(mcol))
    one2one <- human2mouse_df[(hcol %in% ht[ht$Freq == 1,]$hcol) & (mcol %in% mt[mt$Freq == 1,]$mcol),]
    colnames(one2one) <- c('human', 'mouse')
    return(one2one)
}