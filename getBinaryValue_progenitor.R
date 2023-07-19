#Rsansale 3/23/20
#getBinaryValue
#inputs: values, cutoff, outputVal
#outputs: binVal

getBinaryValue <- function(values, cutoff, outputVal){
  
  if(length(cutoff) ==1){
    cutoff = matrix(rep(cutoff,times=nrow(values)),nrow =nrow(values))
  }
binVal <- matrix(rep(0,times = nrow(values)),nrow = nrow(values), ncol = 1)
for(i in 1:nrow(values)){
  binVal[i,] = (values[i,] > cutoff[i])*1
  
  binVal[i,] <- ifelse(binVal[i,] == 1, outputVal[2], outputVal[1])
  
}
return(binVal)
}
