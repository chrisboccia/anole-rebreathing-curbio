matchstick <- function(dataset, metadata, matchable, match_id, as_char=0) {
  add_on <- data.frame(matrix(nrow = length(matchable), ncol=(length(metadata[1,])-1)))
  matchable <- trimws(matchable)
  colnames(add_on) <- colnames(metadata[,-match_id])
  rownames(metadata) <- metadata[,match_id]
  if(length(as_char) > 0){
    for(g in 1: length(as_char)){
      metadata[,as_char[g]] <- as.character(metadata[,as_char[g]])
    }
    
  }
  for (e in 1:length(matchable)){
    add_on[e,] <- metadata[as.character(matchable[e]),-match_id]
    
  }
  
  matched <- cbind(dataset,add_on)
  return(matched)
}##description below
####function for matching up data across multiple datasets####
####dataset = the dataset you want to integrate data into
###metadata = dataframe of data to add onto data set
###both these data sets should be linked by IDs, matchable 
###(should be the column you want to match in dataset); should be in same order as dataset
###script will pull data from metadata to fill out dataset
###match_id = column in metadata to match with
###as_char--pass a vector of columns from metaadata to be interpreted as character (otherwise they will be read as factors)
###run function below to be able to use it
#####

###define var ending function used below####
var_endings <- function(vars, fcat){
  
  pasted <- c(rep("0", length(vars)))
  for(j in 1:length(vars)){
    pasted[j] <- paste(vars[j], fcat)
    
  }
  
  return(pasted)
}
####