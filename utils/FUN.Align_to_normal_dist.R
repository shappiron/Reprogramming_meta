#This function aligns values of the sampling to the normal distribution

Align_to_normal_dist <- function(data,by="cols"){
  if (by=="rows"){
    new_data <- t(data)
  }else{
    new_data <- data
  }
  
  for (i in 1:ncol(new_data)){
    temp_col <- new_data[,i]
    number_values <- sum(!is.na(temp_col))
    normal_dist <- rnorm(number_values)
    normal_dist <- sort(normal_dist)
    temp_col_nona <- temp_col[!is.na(temp_col)]
    temp_col_nona[order(temp_col_nona,decreasing = F)] <- normal_dist
    new_data[which(!is.na(new_data[,i])),i] <- temp_col_nona
  }

  if (by=="rows"){
    new_data <- t(new_data)
  }
  
  return(new_data)
}