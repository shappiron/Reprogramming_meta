# #################################
# This functions performs metafor statistical analysis with specified
# data matrix, descriptions and random terms
# For each gene the function construct INTERCEPTION-ONLY model and calculate
# corresponding statistic

Metafor_analysis <- function(current_data,
                             current_std=F,
                             current_desc,
                             random_terms,
                             LOO=FALSE,
                             Corr1=FALSE){
        library(metafor)
  
        #1) Manage inputs
        data_for_test = as.matrix(current_data)
        if (current_std!=F){
                stds_for_test = as.matrix(current_std)
        }
        description_for_test = current_desc


        #2) Prepare outputs
        output <- data.frame(list(ID = rownames(data_for_test),
                                  coefficient=c(rep(0, nrow(data_for_test))),
                                  se=(rep(0, nrow(data_for_test))),
                                  pvalue=c(rep(1, nrow(data_for_test))),
                                  FDR = c(rep(1, nrow(data_for_test))),
                                  intervention_LOO = c(rep(0, nrow(data_for_test))),
                                  coefficient_LOO = c(rep(0, nrow(data_for_test))),
                                  pvalue_LOO=c(rep(1, nrow(data_for_test))),
                                  FDR_LOO = c(rep(1, nrow(data_for_test))),
                                  intervention_robust = c(rep(0, nrow(data_for_test))),
                                  coefficient_robust = c(rep(0, nrow(data_for_test))),
                                  pvalue_robust=c(rep(1, nrow(data_for_test))),
                                  FDR_robust = c(rep(1, nrow(data_for_test)))),
                                  row.names = rownames(data_for_test))

        #3) Run analysis
        main_errors = 0
        aaa <- Sys.time()
        
        for (g in 1:nrow(data_for_test)){
                #i) Get subset of data
                temp_id <- rownames(data_for_test)[g]
                where_expressed = !is.na(data_for_test[temp_id,])
                current_data = data_for_test[temp_id, where_expressed]
                if (current_std!=F){
                    current_std = stds_for_test[temp_id, where_expressed]
                }
                #get only those datasets which express given gene
                current_description = description_for_test %>% filter(as.vector(where_expressed)) 
                
                if(sum(where_expressed)>2){
                        #ii) run analysis
                        #Specify random terms
                        number_random <- length(random_terms)
                        random_list <- list()
                        for (i in seq_along(random_terms)){
                                #if number of unique factors > 1 - use it
                                if (length(levels(factor(current_description[,random_terms[i]])))>1){
                                        random_list[[i]] <- factor(current_description[,random_terms[i]])
                                }
                        }
                        #random_term_list <- list(~1|random_list[[1]],~1|random_list[[2]],~1|random_list[[3]],
                        #                         ~1|random_list[[4]],~1|random_list[[5]],~1|random_list[[6]])
                        
                        #HARDCODEEEE
                        random_term_list <- list(~1|random_list[[1]], ~1|random_list[[2]])
                        random_term_list <- random_term_list[1:length(random_list)] #why this row needed?
                        
                        if (all(Corr1==F)){
                                if (current_std!=F){
                                        tryCatch({
                                                result = rma.mv(yi=current_data, V= current_std^2,
                                                                random = random_term_list,
                                                                method="REML")
                                                output[temp_id,"coefficient"] <- result$b[1]
                                                output[temp_id,"pvalue"] <- result$pval[1]
                                                output[temp_id,"se"] <- result$se[1]
                                        }, error=function(e){
                                                main_errors=main_errors+1
                                                print("Main\n")
                                                cat("ERROR :",conditionMessage(e), "\n")})
                                        
                                }else{
                                        tryCatch({
                                                result = rma.mv(yi=current_data, V=1,
                                                                random = random_term_list,
                                                                method="REML")
                                                output[temp_id,"coefficient"] <- result$b[1]
                                                output[temp_id,"pvalue"] <- result$pval[1]
                                                output[temp_id,"se"] <- result$se[1]
                                        }, error=function(e){
                                                main_errors=main_errors+1
                                                print("Main\n")
                                                cat("ERROR :",conditionMessage(e), "\n")})
                                }        
                        }else{
                                temp_Corr1 <- Corr1[levels(factor(current_description[,random_terms[1]])),
                                               levels(factor(current_description[,random_terms[1]]))]

                                if (current_std!=F){
                                        tryCatch({
                                                result = rma.mv(yi=current_data, V= current_std^2,
                                                                random = random_term_list,
                                                                R=list('random_list[[1]]'=temp_Corr1),
                                                                Rscale=0,
                                                                method="REML")
                                                output[temp_id,"coefficient"] <- result$b[1]
                                                output[temp_id,"pvalue"] <- result$pval[1]
                                                output[temp_id,"se"] <- result$se[1]
                                        }, error=function(e){
                                                main_errors=main_errors+1
                                                print("Main\n")
                                                cat("ERROR :",conditionMessage(e), "\n")})
                                        
                                }else{
                                        tryCatch({
                                                result = rma.mv(yi=current_data, V=1,
                                                                random = random_term_list,
                                                                R=list('random_list[[1]]'=temp_Corr1),
                                                                Rscale=0,
                                                                method="REML")
                                                output[temp_id,"coefficient"] <- result$b[1]
                                                output[temp_id,"pvalue"] <- result$pval[1]
                                                output[temp_id,"se"] <- result$se[1]
                                        }, error=function(e){
                                                main_errors=main_errors+1
                                                print("Main\n")
                                                cat("ERROR :",conditionMessage(e), "\n")})        
                                }
                                
                        }
                               

                        #iii) Run robust and LOO test
                        if (LOO!=F){
                                pvals_robust <- c()
                                coeffs_robust <- c()
                                ints_robust <- c()
                                
                                LOO_factor <- factor(current_description[,LOO])
                                for (s in levels(LOO_factor)){
                                        remaining_samples <- rownames(current_description[current_description[,LOO]!=s,])
                                        temp_data <- current_data[remaining_samples]
                                        if (current_std!=F){
                                                temp_std <- current_std[remaining_samples]
                                        }
                                        temp_description <- current_description[remaining_samples,]

                                        temp_random_list <- list()
                                        test_corr=0
                                        for (i in seq_along(random_terms)){
                                                if (length(levels(factor(temp_description[,random_terms[i]])))>1){
                                                        temp_random_list[[i]] <- factor(temp_description[,random_terms[i]])
                                                }else{
                                                        if(i==1){
                                                                test_corr=1
                                                        }
                                                }
                                        }

                                        #temp_random_term_list <- list(~1|temp_random_list[[1]],~1|temp_random_list[[2]],~1|temp_random_list[[3]],
                                        #                         ~1|temp_random_list[[4]],~1|temp_random_list[[5]],~1|temp_random_list[[6]])
                                        temp_random_term_list <- list(~1|temp_random_list[[1]],~1|temp_random_list[[2]])
                                        temp_random_term_list <- temp_random_term_list[1:length(temp_random_list)]
                                        if (all(Corr1==F) || test_corr==1){
                                                if (current_std!=F){
                                                        tryCatch({
                                                                result = rma.mv(yi=temp_data, V= temp_std^2,
                                                                                random = temp_random_term_list,
                                                                                method="REML")
                                                                coeffs_robust <- c(coeffs_robust,result$b[1])
                                                                pvals_robust <- c(pvals_robust,result$pval[1])
                                                                ints_robust <- c(ints_robust,s)
                                                        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                                                }else{
                                                        tryCatch({
                                                                result = rma.mv(yi=temp_data, V=1,
                                                                                random = temp_random_term_list,
                                                                                method="REML")
                                                                coeffs_robust <- c(coeffs_robust,result$b[1])
                                                                pvals_robust <- c(pvals_robust,result$pval[1])
                                                                ints_robust <- c(ints_robust,s)
                                                        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                                                }
                                        }else{
                                                temp_Corr1 <- Corr1[levels(factor(temp_description[,random_terms[1]])),
                                                               levels(factor(temp_description[,random_terms[1]]))]
                                                
                                                if (current_std!=F){
                                                        tryCatch({
                                                                result = rma.mv(yi=temp_data, V= temp_std^2,
                                                                                random = temp_random_term_list,
                                                                                R=list('temp_random_list[[1]]'=temp_Corr1),
                                                                                Rscale=0,
                                                                                method="REML")
                                                                coeffs_robust <- c(coeffs_robust,result$b[1])
                                                                pvals_robust <- c(pvals_robust,result$pval[1])
                                                                ints_robust <- c(ints_robust,s)
                                                        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                                                }else{
                                                        tryCatch({
                                                                result = rma.mv(yi=temp_data, V=1,
                                                                                random = temp_random_term_list,
                                                                                R=list('temp_random_list[[1]]'=temp_Corr1),
                                                                                Rscale=0,
                                                                                method="REML")
                                                                coeffs_robust <- c(coeffs_robust,result$b[1])
                                                                pvals_robust <- c(pvals_robust,result$pval[1])
                                                                ints_robust <- c(ints_robust,s)
                                                        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                                                }
                                        }

                                }
                                 
                                if (length(pvals_robust)!=0){
                                        output[temp_id,"pvalue_robust"] <- max(pvals_robust)
                                        output[temp_id,"intervention_robust"] <- ints_robust[which.max(pvals_robust)]
                                        output[temp_id,"coefficient_robust"] <- coeffs_robust[which.max(pvals_robust)]
                                        output[temp_id,"pvalue_LOO"] <- min(pvals_robust)
                                        output[temp_id,"intervention_LOO"] <- ints_robust[which.min(pvals_robust)]
                                        output[temp_id,"coefficient_LOO"] <- coeffs_robust[which.min(pvals_robust)]
                                }
                        }
                }
        }
        bbb <- Sys.time()
        print(bbb-aaa)
        print(main_errors)
        
        output$FDR <- p.adjust(output$pvalue,method="BH")
        if (LOO!=F){
                output$FDR_robust <- p.adjust(output$pvalue_robust,method="BH")
                output$FDR_LOO <- p.adjust(output$pvalue_LOO,method="BH")
        }
        return(output)

}



