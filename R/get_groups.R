#### get_groups
# Returning a vector of group belongings based on a filename vector

get_groups <- function(filenames) {
  groups <- rep(NA,length(filenames))
  groups[which(grepl("sQC", filenames))]  <- "sQC"
  groups[which(grepl("ltQC", filenames))] <- "ltQC"
  groups[is.na(groups)] <- "sample"
  return(groups)
}