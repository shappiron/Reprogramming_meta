#This function performs functional GSEA with KEGG, REACTOME and GO BP ontology

Run_functional_GSEA <- function(ranked_list,species="Mus musculus",type="Symbol",
                                n_permutations=1000,min_size=15,max_size=500,
                                with_MF=F, ontologies_list=NaN, ontologies_dict=NaN){

  ranked_list <- ranked_list[order(ranked_list[,2],decreasing = T),]
  ranked_vector <- ranked_list[,2]
  names(ranked_vector) <- as.character(ranked_list[,1])
  
  library(fgsea)
  library(data.table)
  library(ggplot2)
  library(msigdbr)

  m_df_h = msigdbr(species = species, category = "H")
  m_df_bp = msigdbr(species = species, category = "C5",subcategory = "BP")
  m_df_kegg = msigdbr(species = species, category = "C2",subcategory = "CP:KEGG")
  m_df_reactome = msigdbr(species = species, category = "C2",subcategory = "CP:REACTOME")
  m_df <- rbind(rbind(m_df_h, m_df_bp,m_df_kegg),m_df_reactome)
  
  if (with_MF==T){
    m_df_mf <- msigdbr(species = species, category = "C5",subcategory = "MF")
    m_df <- rbind(m_df,m_df_mf)
  }

  if (typeof(ontologies_dict) == 'logical'){
    if (ontologies_list==F){
      functions_list <- levels(factor(m_df$gs_name))
    }else{
      functions_list <- levels(factor(ontologies_list))
    }

    functions_annotation <- list()
    for (i in functions_list){
      temp_data <- m_df[m_df$gs_name==i,]
      if (type=="Symbol"){
        temp_data <- toupper(as.character(temp_data$gene_symbol))
      }else if (type=="Entrez"){
        temp_data <- as.character(temp_data$entrez_gene)      
      }
      functions_annotation[[i]] <- temp_data
    }
  }else{
    functions_annotation <- ontologies_dict
  }

  fgseaRes <- fgsea(pathways = functions_annotation, 
                    stats    = ranked_vector,
                    nperm = n_permutations,
                    minSize  = min_size,
                    maxSize  = max_size)

  return(fgseaRes)
}

