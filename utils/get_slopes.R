## This function aggregate time-course datasets into corresponding slopes vector by limma lmFit
## the result consist of two objects: list of datasets with statistics ["logFC", "SE", "P.Value", "FDR"]
## and dataframe with slopes only for full intersection of genes within datasets

#stoplist example:
#c("GSE116309$OKMS", "GSE10871$OSKM", 
#  "GSE116309$OK+9MS", "GSE38509$GFP", 
#  "GSE38509$OSK", "GSE46321$C/EBPα- OSKM")
#genelist example: #merged_clusters #treat_vs_control_intersection

get_slopes <- function(ultradf, stoplist=c(), genelist=c()){
    slopes <- data.frame(row.names=rownames(ultradf[[1]]$data), tmp=rownames(ultradf[[1]]$data))
    fdrs <- data.frame(row.names=rownames(ultradf[[1]]$data), tmp=rownames(ultradf[[1]]$data))
    datasets_list <- c()
    for (key in names(ultradf)){
        if (key %in% stoplist){
            next
        }
        
        #extract data
        if (length(genelist) != 0){
            df <- ultradf[[key]]$data[genelist,]
            df <- df[complete.cases(df),]
        }else{
            df <- ultradf[[key]]$data
        }
        pheno <- ultradf[[key]]$pheno
        time <-ultradf[[key]]$time
        
        #create design matrix
        design <- data.frame(row.names=rownames(pheno))
        design['Intercept'] <- rep(1, length(rownames(pheno)))
        design['Day'] <- time
        
        #analyze
        fit <- lmFit(df, design=design)
        fit <- eBayes(fit)
        result <- topTable(fit, adjust.method = "BH", sort.by="none", coef="Day", n=Inf, confint=T)
        result[key] <- result["logFC"]
        result$SE <- (result[,3] - result[,2]) / 3.92
        result$FDR <- result$adj.P.Val
        
        datasets_list[[key]] <- result[c("logFC", "SE", "P.Value", "FDR")]
        
        cat(key, "Pheno:", nrow(pheno), "Total:", nrow(result), "Passed:", nrow(result[result['adj.P.Val'] < 0.05,]), '\n')
        
        #merge
        slopes <- merge(slopes, result[key], by="row.names", sort=FALSE)
        row.names(slopes) <- slopes$Row.names
        slopes <- slopes[, -1]
    }
    slopes <- slopes[, -1]
    return(list("datasets"=datasets_list, "slopes"=slopes))

}