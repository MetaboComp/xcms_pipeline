#### retrieveBatch
# retrieve peakTable of the specified batch together with meta data from meta data table

retrieveBatch <- function(peakTable,
                          meta,
                          select){
  peakTable = peakTable[grepl(select, rownames(peakTable), ignore.case = T), ]
  if (is.null(dim(meta))) 
    meta = matrix(meta, ncol = 1)
  meta = meta[grepl(select, meta$fullname, ignore.case = T), ]
  return(list(peakTable = peakTable, meta = meta))
}