#this function is for gathering metafor results to a dataframe
gather <- function(data, col, name, subcol=NA){
    cont <- data.frame(row.names=as.character(names(data)))
    for (g in names(data)){
        if (!is.na(subcol)){
            if (col == "fit.stats"){
                #cont[g, name] <- flatten(data[[g]][[col]])[[subcol]]
                cont[g, name] <- unlist(data[[g]][[col]])[[subcol]]
            }else{
                vec <- data[[g]][[col]]
                if (subcol <= length(vec)){
                    cont[g, name] <- vec[[subcol]]
                }
            }
        }else{
            cont[g, name] <- data[[g]][[col]]
        }
        
    }  
    return(cont)
}