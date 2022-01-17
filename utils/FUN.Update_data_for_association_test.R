##########################################
#This script performs analysis of individual datasets for association study

Update_data_for_association_test <- function(data_list,data_norm,
                                             phenodata,
                                             group_col="Group",
                                             control_group,
                                             batch_col="",
                                             current_gse=F,
                                             function_path=function_path,
                                             type="microarray",
                                             Entrez=T){
  
    source(paste0(function_path,"FUN.Dataset_DEG_analysis.R"))

    phenodata$Group <- phenodata[,group_col]
    temp_results <- Dataset_DEG_analysis(sample_eset = data_norm,phenodata=phenodata,
                                         group_col="Group",control_group=control_group,
                                         batch_col=batch_col,Entrez=Entrez,
                                         function_path = function_path,type = type)
    names(temp_results)

    for (b in seq_along(names(temp_results))){
      if (current_gse!=F){
        names(temp_results)[b] <- paste(current_gse,names(temp_results)[b],sep="_")
      }
      names(temp_results)[b] <- gsub("_-","",names(temp_results)[b])
      names(temp_results)[b] <- gsub(" ","",names(temp_results)[b])
    }

    if (length(data_list)==0){
      data_list <- temp_results
    }else if(class(data_list)=="list"){
      data_list <- c(data_list,temp_results)
    }else{
      data_list <- temp_results
    }
    
    return(data_list)
}  