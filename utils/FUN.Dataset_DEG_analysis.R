#################################################################################################
#This script performs preprocessing of dataset 
#and analysis of differentially expressed genes (using limma or edgeR)

Dataset_DEG_analysis <- function(sample_eset,phenodata,
                                 group_col="Group",control_group,
                                 batch_col="",
                                 Entrez=T,
                                 type="microarray",
                                 function_path="D:/Documents/Science/PhD/Functions/"){
  library(matrixStats)
  library(genefilter)
  dataset_results <- list()
  
  if(type=="microarray"){
    if (Entrez==T){
    #1. Transform to Entrez
    data_Entrez <- sample_eset[fData(sample_eset)[,grepl("entrez",tolower(colnames(fData(sample_eset))))]!="" & !is.na(fData(sample_eset)[,grepl("entrez",tolower(colnames(fData(sample_eset))))]),]
    print(paste0("There were ", nrow(sample_eset)," probes"))
    print(paste0("There are ", nrow(data_Entrez)," probes with Entrez IDs"))
    
    source(paste0(function_path,"FUN.Microarray_to_Entrez.R"))
    source(paste0(function_path,"FUN.Microarray_dictionary_create.R"))
    
    dictionary_output <- Microarray_dictionary_create(data_Entrez)
    data_Entrez <- Microarray_to_Entrez(dictionary_output$dataset,dictionary_output$dataset_Entreztoprobe,dictionary_output$dataset_probetoEntrez)
    }else{
      data_Entrez <- sample_eset
    }


    #2. Calculate DEG using limma
    source(paste0(function_path,"FUN.Run_limma_for_dataset.R"))

    output <- Run_limma_for_dataset(data_Entrez,phenodata = phenodata,
                                    control_group=control_group,
                                    batch_col=batch_col)

  }else if(type=="RNAseq"){
      #0. RLE normalization
    source(paste0(function_path,"FUN.RLE_normalization.R"))
    data_RLE <- RLE_normalization(sample_eset)
    #plot(density(log2(exprs(data_RLE)+1)[,1]))
    #for (i in 2:ncol(data_RLE)){
    #  lines(density(log2(exprs(data_RLE)+1)[,i]),col=i)
    #}
    
    if (Entrez==T){
      #1. Transform to Entrez
      source(paste0(function_path,"FUN.Ensembl_mouse_dictionary_create.R"))
      source(paste0(function_path,"FUN.Ensembl_to_entrez.R"))
      EntreztoID_rnaseq <- Ensembl_mouse_dictionary_create(sample_eset)
      data_Entrez <- Ensembl_to_entrez(sample_eset,EntreztoID_rnaseq)
      data_Entrez <- new("ExpressionSet", exprs=data_Entrez,phenoData=phenoData(sample_eset))
      #data_RLE_Entrez <- Ensembl_to_entrez(data_RLE,EntreztoID_rnaseq)
      #data_RLE_Entrez <- new("ExpressionSet", exprs=data_RLE_Entrez,phenoData=phenoData(sample_eset))
      
      print(paste0("There were ", nrow(sample_eset)," probes"))
      print(paste0("There are ", nrow(data_Entrez)," probes with Entrez IDs"))
      
    }else{
      data_Entrez <- sample_eset
      data_RLE_Entrez <- data_RLE
    }

    gene_list <- rownames(data_Entrez)

    #2. Calculate DEG using edgeR
    source(paste0(function_path,"FUN.Run_edgeR_for_dataset.R"))
    source(paste0(function_path,"FUN.Run_edgeR.R"))
    library(edgeR)
    library(statmod)

    output <- Run_edgeR_for_dataset(data_Entrez,phenodata = phenodata,
                                    control_group=control_group,
                                    batch_col=batch_col)
    
  }else{
    print("Unknown type of data provided")
    stop()
  }
  
  return(output)
}