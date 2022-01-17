####################################################################
#This function maps Entrez IDs from one species to the orthologs of the other
#(such that only mutually exclusive orthologs are left)

Map_to_orthologs <- function(gene_list,species_from,species_to,type="Entrez"){
  if (species_from %in% c("mmusculus","Mus musculus","Mouse","Mice","Mus")){
    species_from <- "mmusculus"
  }else if(species_from %in% c("hsapiens","Homo sapiens","Human","Humans","Homo")){
    species_from <- "hsapiens"
  }

  if (species_to %in% c("mmusculus","Mus musculus","Mouse","Mice","Mus")){
    species_to <- "mmusculus"
  }else if(species_to %in% c("hsapiens","Homo sapiens","Human","Humans","Homo")){
    species_to <- "hsapiens"
  }
  
  #A) Find orthologs
  library(biomaRt)
  ensembl <- useMart("ensembl")
  species_1 <- paste0(species_from,"_gene_ensembl")
  species_2 <- paste0(species_to,"_gene_ensembl")
  dataset_from = useDataset(species_1, mart=ensembl)
  dataset_to = useDataset(species_2, mart=ensembl)
  if (type=="Entrez"){
    From_to_to_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                   values=gene_list,
                                   mart=dataset_from,attributesL=c("entrezgene_id"), martL=dataset_to)
  }else if (type=="Symbol" &  species_from=="hsapiens"){
    From_to_to_orthologs <- getLDS(attributes=c("hgnc_symbol","entrezgene_id"), filters="hgnc_symbol",
                                   values=gene_list,
                                   mart=dataset_from,attributesL=c("entrezgene_id"), martL=dataset_to)[,c(2:3)]
  }else if (type=="Symbol" & species_from=="mmusculus"){
    From_to_to_orthologs <- getLDS(attributes=c("mgi_symbol","entrezgene_id"), filters="mgi_symbol",
                                   values=gene_list,
                                   mart=dataset_from,attributesL=c("entrezgene_id"), martL=dataset_to)[,c(2:3)]
  }else{
    stop("Invalid type")
  }
  colnames(From_to_to_orthologs) <- c("From_Entrez","To_Entrez")

  #B) Leave only genes that uniquely map to orthologs in the second species
  #and which are uniquely mapped to orthologs in the first species
  From_to_unique <- From_to_to_orthologs
  dupl_from_Entrez <- unique(From_to_to_orthologs$From_Entrez[which(duplicated(From_to_to_orthologs$From_Entrez))])
  unique_from_Entrez <- unique(From_to_to_orthologs$From_Entrez[!(From_to_to_orthologs$From_Entrez %in% dupl_from_Entrez)])
  print(paste0("Genes with non-unique orthologs: ",length(dupl_from_Entrez)))
  print(paste0("Genes with unique orthologs: ",length(unique_from_Entrez)))

  dupl_to_Entrez <- unique(From_to_to_orthologs$To_Entrez[which(duplicated(From_to_to_orthologs$To_Entrez))])
  unique_to_Entrez <- unique(From_to_to_orthologs$To_Entrez[!(From_to_to_orthologs$To_Entrez %in% dupl_to_Entrez)])

  From_to_unique <- From_to_to_orthologs[From_to_to_orthologs$From_Entrez %in% unique_from_Entrez & From_to_to_orthologs$To_Entrez %in% unique_to_Entrez,]
  print(paste0("Mutually exclusive orthologs: ",nrow(From_to_unique)))
  rownames(From_to_unique) <- From_to_unique$From_Entrez
  return(From_to_unique)
}