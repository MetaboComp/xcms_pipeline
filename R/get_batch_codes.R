#### get_batch_codes
# Collect batch codes from filenames

get_batch_codes <- function(filenames) {
  batch_codes <- stri_extract(filenames, regex = "B[0-9]+W[0-9]+")
  return(batch_codes)
}