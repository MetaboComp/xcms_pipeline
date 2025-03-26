#### scaleAwayFromZero
# Scales parameters up above < 0 

scaleAwayFromZero <- function(PTBefore, PTAfter){
  
  PTFinal <- data.frame(matrix(nrow=dim(PTAfter)[1], ncol=dim(PTAfter)[2]))
  colnames(PTFinal) <- colnames(PTAfter)
  rownames(PTFinal) <- rownames(PTAfter)
  
  #Checking if features lost during drift correction
  if(dim(PTBefore)[2] == dim(PTAfter)[2]){
    for(i in 1:ncol(PTBefore)){
      dataStandardized <- (PTAfter[,i] - min(PTAfter[,i])) / (max(PTAfter[,i]) - min(PTAfter[,i]))
      PTFinal[,i] <- dataStandardized * (max(PTBefore[,i]) - min(PTBefore[,i])) + min(PTBefore[,i])
    }
  } else {
    for(i in 1:ncol(PTAfter)){
      tempFeature <- PTBefore[,which(colnames(PTBefore) %in% colnames(PTAfter))[i]]
      
      dataStandardized <- (PTAfter[,i] - min(PTAfter[,i])) / (max(PTAfter[,i]) - min(PTAfter[,i]))
      PTFinal[,i] <- dataStandardized * (max(PTBefore[,which(colnames(PTBefore) %in% colnames(PTAfter))[i]]) - min(PTBefore[,which(colnames(PTBefore) %in% colnames(PTAfter))[i]])) + min(PTBefore[,which(colnames(PTBefore) %in% colnames(PTAfter))[i]])
    }
  }
  
  return(PTFinal)
}