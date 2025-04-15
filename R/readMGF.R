###readMGF####
readMGF <- function(filename,
                    pdata = NULL,
                    centroided = TRUE,
                    smoothed = FALSE,
                    verbose = isMSnbaseVerbose(),
                    cache = 1,
                    delZero = FALSE) {
  if (verbose)
    cat("Scanning", filename, "...\n")
  
  mgf <- scan(file = filename, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)
  
  ## From http://www.matrixscience.com/help/data_file_help.html#GEN
  ## Comment lines beginning with one of the symbols #;!/ can be included,
  ## but only outside of the BEGIN IONS and END IONS statements that delimit an MS/MS dataset.
  cmts <- grep("^[#;!/]", mgf)
  if (length(cmts))
    mgf <- mgf[-cmts]
  
  begin <- grep("BEGIN IONS", mgf) + 1L
  end <- grep("END IONS", mgf) - 1L
  
  n <- length(begin)
  
  if (verbose) {
    cnt <- 1L
    pb <- txtProgressBar(min = 0L, max = n, style = 3L)
  }
  
  spectra <- vector("list", length = n)
  fdata <- vector("list", length = n)
  
  for (i in seq(along = spectra)) {
    if (verbose) {
      setTxtProgressBar(pb, cnt)
      cnt <- cnt + 1L
    }
    specInfo <- MSnbase:::extractMgfSpectrum2Info(mgf = mgf[begin[i]:end[i]],
                                                  centroided = centroided)
    if(any(grepl("RTINMINUTES", names(specInfo$fdata))) && specInfo$spectrum@rt==0){
      specInfo$spectrum@rt <- as.double(specInfo$fdata[grepl("RTINMINUTES", names(specInfo$fdata))])*60
    }
    
    spectra[[i]] <- specInfo$spectrum
    fdata[[i]] <- specInfo$fdata
  }
  
  whichZero<-""
  
  ##If user wants to remove all 0 intensity fragments this part will activate and save some space
  if(delZero==TRUE){
    for(i in 1:length(spectra)){
      zeros<-which(spectra[[i]]@intensity==0)
      if(length(zeros)>0){
        spectra[[i]]@intensity<-spectra[[i]]@intensity[-zeros]
        spectra[[i]]@mz<-spectra[[i]]@mz[-zeros]
      }
    }
  }
  
  return(spectra)
}