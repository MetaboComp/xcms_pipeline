#### add_group_column
# Adds a group column to LCMS_data object with sample, sQC and ltQC and blank

add_group_column <- function(LCMS_data) { 
  LCMS_data$group <- rep(NA,nrow(LCMS_data))
  LCMS_data$group[LCMS_data$sample %>% grep("sQC",ignore.case = T,.)] <- "sQC"
  LCMS_data$group[LCMS_data$sample %>% grep("ltQC",ignore.case=T,.)] <- "ltQC"
  LCMS_data$group[is.na(LCMS_data$group)] <- "sample"
  return(LCMS_data)
}