###filterMGF####
# filename = .mgf file from msdial, mse to process
# to_keep = the masses and rts in a dataframe that you're looking for
# dppm = delta parts per million mass difference between precursors and masses of interest
# drt = delta retention time is difference in rt between precursors and masses of interest
# int_filter = intensity under which fargments will be filtered out of MS2s
# centroided = pass-on parameter for MSnbase:::extractMgfSpectrum2Info
# verbose = If

filterMGF <- function(filename,
                      result_dir,
                      to_keep,
                      dppm=10,
                      drt=NULL,
                      int_filter=NULL,
                      centroided = TRUE){
  
  options(digits=10)
  
  mgf <- scan(file = filename, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)
  
  cmts <- grep("^[#;!/]", mgf)
  if (length(cmts))
    mgf <- mgf[-cmts]
  
  begin <- grep("BEGIN IONS", mgf) + 1L
  end <- grep("END IONS", mgf) - 1L
  
  n <- length(begin)
  
  # spectra <- vector("list", length = n)
  fdata <- vector("list", length = n)
  whichToKeep <- c()
  
  for (i in seq(along = c(1:n))) {
    suppressWarnings(specInfo <- MSnbase:::extractMgfSpectrum2Info(mgf = mgf[begin[i]:end[i]],
                                                                   centroided = centroided))
    
    if(length(strsplit(specInfo$fdata[grepl("PEPMASS", names(specInfo$fdata))], " ")[[1]]) > 1){
      curr_mz <- as.double(strsplit(specInfo$fdata[grepl("PEPMASS", names(specInfo$fdata))], " ")[[1]][1])
    } else {
      curr_mz <- as.double(specInfo$fdata[grepl("PEPMASS", names(specInfo$fdata))])
    }
    
    curr_rt <- as.double(specInfo$fdata[grepl("RT", names(specInfo$fdata))])
    if(!any(grepl("INSECONDS", names(specInfo$fdata)))){
      curr_rt <- curr_rt * 60
    }
    dmz <- curr_mz/10^6*dppm
    
    which_match <- which(to_keep[,1] > curr_mz - dmz &
                           to_keep[,1] < curr_mz + dmz)
    
    
    if(length(which_match)!=0){
      if(!is.null(drt)){
        which_match <- which_match[which(to_keep[which_match, 2] > curr_rt - drt &
                                           to_keep[which_match, 2] < curr_rt + drt)]
      }
      
      if(length(which_match) > 0){
        whichToKeep <- c(whichToKeep,
                         i)
      }
    }
  }
  
  
  
  
  
  if(is.null(whichToKeep)){
    message("Problem: No matching MS2 spectra. Cannot continue comparisons.")
    return(F)
  }
  
  file_name <- paste0(dirname(filename), result_dir,"/Filtered_", basename(filename))
  for(i in 1:length(whichToKeep)){
    spectra <- mgf[begin[whichToKeep[i]]:end[whichToKeep[i]]]
    
    if(any(grepl("Num peaks:", spectra, ignore.case=T))){
      which_num_peaks <- which(grepl("Num peaks:", spectra, ignore.case=T))
      spectra <- spectra[!grepl("Num peaks:", spectra, ignore.case=T)]
      spectra_start <- which_num_peaks
    } else if(any(grepl("CHARGE", spectra, ignore.case=T))){
      spectra_start <- which(grepl("CHARGE", spectra, ignore.case=T))+1
    } else {
      spectra_start <- which(grepl("PEPMASS", spectra, ignore.case=T))+1
    }
    
    if(!is.null(int_filter)){
      which_to_remove <- c()
      
      if(grepl("\\t", spectra[spectra_start], ignore.case=T)){
        split_pattern <- "\\t"
      } else if (grepl(" ", spectra[spectra_start], ignore.case=T)){
        split_pattern <- " "
      }
      
      for(j in spectra_start:length(spectra)){
        if(as.integer(strsplit(spectra[j], split_pattern)[[1]][2]) < int_filter){
          which_to_remove <- c(which_to_remove,
                               j)
        }
      }
      
      if(length(which_to_remove) > 0){
        spectra_start_n <- length(spectra)
        spectra <- spectra[-which_to_remove]
        if(length(spectra) == (spectra_start-1)){
          next
        }
      }
    }
    
    
    
    write("BEGIN IONS", file=file_name, append=T)
    write(spectra, file=file_name, append=T)
    write("END IONS", file=file_name, append=T)
  }
  
  return(T)
}