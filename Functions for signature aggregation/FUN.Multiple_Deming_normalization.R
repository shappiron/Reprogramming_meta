# This script performs multiple Deming minimization and provides Deming output,
# best coefficients along with quality figures

Multiple_Deming_normalization <- function(logFC_list,
                                          number_top=0,
                                          coeff_thresholds=c(0.1),
                                          se_table=F,
                                          FDR_threshold=0.05,
                                          number_iter=1, 
                                          max_value=F,
                                          function_path="D:/Documents/Science/PhD/Functions/"){
  
  output <- list()

  #a) Plot densities
  show_density <- function(temp, ylim=c(0, 3), title=""){
    plot(density(temp[[1]]$logFC), col=1, xlim=c(-3,3), ylim=ylim, main=title)
    for (i in 1:length(temp)){
      lines(density(temp[[i]]$logFC), col=i)
    }
  }
  show_density(logFC_list,c(0,8),title="Distribution at the beginning")

  
  #b) Normalize by sd
  logFC_list_sd <- logFC_list
  for (i in names(logFC_list_sd)){
    if (se_table==T){
      logFC_list_sd[[i]]$SE <- logFC_list_sd[[i]]$SE / sd(logFC_list_sd[[i]]$logFC) #????&&&&&&
    }
    logFC_list_sd[[i]]$logFC <- logFC_list_sd[[i]]$logFC / sd(logFC_list_sd[[i]]$logFC)
  }
  show_density(logFC_list_sd,c(0,5),title="Distribution after scaling")

    
  #c) Run Deming regression
  coeff_threshs <- coeff_thresholds
  deminglist <- list()
  minimums <- c()
  for (a in 1:length(coeff_threshs)){
    print("Starting with coefficient threshold:")
    print(a)
    for (i in 1:number_iter){
      print(paste0("Iteration: ", i))
      deminglist[[(a-1) * number_iter + i]] <- Deming_minimizer(logFC_list = logFC_list_sd,
                                                  number_top = number_top,
                                                  coeff_threshold = coeff_threshs[a], #slopes correlation coeff threshold
                                                  FDR_threshold = FDR_threshold,
                                                  function_path=function_path)
      minimums = c(minimums, deminglist[[(a-1)*number_iter+i]]$minimum)
    }
  }
  
  output$Deming_list <- deminglist
  output$Minimums <- minimums
  

  #d) Visualize coefficient consistency in boxplot format
  source(paste(function_path,"FUN.Calculate_paired_correlation.R",sep=""))
  
  max_abs_score <- max(output$Minimums)
  min_abs_score <- min(output$Minimums)
  
  #table with normalization coefficients for each dataset
  a = as.data.frame(deminglist[[1]]$coefs)
  for (i in 2:length(deminglist)){
    a = cbind(a, deminglist[[i]]$coefs)
  }
  colnames(a) = 1:length(deminglist)
  rownames(a) <- names(logFC_list_sd)
    
  a$group <- row.names(a)
  a.m <- melt(a, id.vars = "group")
  b = data.frame()
  for (i in 1:length(deminglist)){
    b = rbind(b, deminglist[[i]]$minimum)
  }
  b = cbind(b, 1:length(deminglist))
  b = cbind(b, rep(coeff_thresholds, each=number_iter))
  colnames(b) = c("minimum", "variable", "Threshold")
  b$variable = as.factor(b$variable)
  library(dplyr)
  a.m$minimum = left_join(a.m, b, by = "variable")
  a.m$minimum$group <- factor(a.m$group,levels=names(logFC_list_sd))
  
  if (max_value!=F){
    figure <- ggplot(a.m$minimum, aes(group, log10(value)))+
      geom_boxplot()+
      geom_hline(yintercept = 0,lwd=1)+
      #geom_hline(yintercept = 1,lwd=1,lty=2)+
      coord_cartesian(ylim=c(-max_value,max_value))+
      xlab("Signature")+
      ylab("Normalization coefficient (log10)")+
      geom_jitter(aes(color=minimum,shape = as.factor(Threshold)),size=3,alpha=0.7)+
      theme_bw()+
      scale_shape_discrete(name="Rho threshold",
                           labels=coeff_thresholds)+
      scale_color_gradient2(high="blue3",low="red3",mid = "red3",
                            guide = "colorbar",
                            midpoint = min_abs_score,
                            limits=c(min_abs_score,max_abs_score))+
      theme(strip.background = element_rect(fill="white"),
            axis.text.x = element_text(size=12,angle=45,hjust=1,colour="black"),
            axis.text.y = element_text(size=12,colour="black"),
            axis.title=element_text(size=16),
            title=element_text(size = 18),
            #legend.position = c(0.91,0.8),
            #legend.justification = c(0.73,0.15),
            #legend.title=element_text(size=18),
            #legend.text = element_text(size=14),
            plot.title = element_text(hjust=0.5))
    
  }else{
    figure <- ggplot(a.m$minimum, aes(group, log10(value)))+
      geom_boxplot()+
      geom_hline(yintercept = 0,lwd=1)+
      #geom_hline(yintercept = 1,lwd=1,lty=2)+
      coord_cartesian()+
      xlab("Signature")+
      ylab("Normalization coefficient (log10)")+
      geom_jitter(aes(color=minimum,shape = as.factor(Threshold)),size=3,alpha=0.7)+
      theme_bw()+
      scale_shape_discrete(name="Rho threshold",
                           labels=coeff_thresholds)+
      scale_color_gradient2(high="blue3",low="red3",mid = "red3",
                            guide = "colorbar",
                            midpoint = min_abs_score,
                            limits=c(min_abs_score,max_abs_score))+
      theme(strip.background = element_rect(fill="white"),
            axis.text.x = element_text(size=12,angle=45,hjust=1,colour="black"),
            axis.text.y = element_text(size=12,colour="black"),
            axis.title=element_text(size=16),
            title=element_text(size = 18),
            #legend.position = c(0.91,0.8),
            #legend.justification = c(0.73,0.15),
            #legend.title=element_text(size=18),
            #legend.text = element_text(size=14),
            plot.title = element_text(hjust=0.5))
    
  }

  print(figure)
  output$Quality <- figure
  
  #e) Normalize data
  library(deming)
  best_coefs <- which.min(minimums)
  final_deming <- deminglist[[best_coefs]]
  final_deming_coefs <- final_deming$coefs
  names(final_deming_coefs) <- names(logFC_list_sd)

  output$Final_deming_coefs <- final_deming_coefs

  
  #f) Normalize by deming coefficients
  logFC_list_sd_norm <- logFC_list_sd
  for (a in names(logFC_list_sd_norm)){
    if (se_table==T){
      logFC_list_sd_norm[[a]]$SE <- logFC_list_sd_norm[[a]]$SE / final_deming_coefs[a]
    }
    logFC_list_sd_norm[[a]]$logFC <- logFC_list_sd_norm[[a]]$logFC / final_deming_coefs[a]
  }
  show_density(logFC_list_sd_norm,c(0,5), title="Distribution after Deming normalization")
  
 
  output$Normalized_data <- logFC_list_sd_norm   

    
  return(output)
}
#output of structure:
# Normalized_data - list of DataFrames of structure (logFC, SE, Pval, FDR)
# Deming_list - list of deming coefficients and objective minimums obtained at different iterations
# Quality - barplot with obtained coefficients distribution
# Minimums - vector of objective minimums obtained at different iterations
# Final_deming_coefs - best deming coefficients chosen based on the best objective values
