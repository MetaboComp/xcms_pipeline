#### get_batch_nos
# Returning a vector of batch numbers as numbers only

get_batch_nos <- function(filenames) {
  batch_codes <- stri_extract(filenames, regex="B[0-9]+W[0-9]+")
  bcode <- substring(batch_codes, 2, nchar(batch_codes))
  bn <- stri_extract(bcode, regex="^[0-9]+")
  return(bn)
}