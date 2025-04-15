###cosineSim####
#' cosineSim - A function for comparing 2 ms2 spectra, one from a SQLite DB and the other from user inputted data
#'
#' @param combList A "combList" object to be used for the analysis
#' @param metaData MetaData formatted according to example file (see matchMSfunction() & getChemForm())
#' @param dbName Name / path to a SQLite of the same format as the example file. To make a new one see function "buildDB()"
#'
#' @return A list of mgf-objects resulting from MS1 vs Theoretical dPPM match + Precursor vs Theoretical dPPM match + MS1 vs Precursor dPPM*2 match + MS1 RT vs MS2 RT within rtDiff
#' @export
#'

#NIST database mzWeight = 3, intWeight = 0.6; MassBank mzWeight = 2, intWeight = 0.5.

cosineSim <- function(ms2SampSpec, ms2StdSpec, dPPM=10, mzWeight=2, intWeight=0.5, mzPrec=NULL, ms2StdMetaData=NULL, absThreshold=NULL, precFilter=TRUE){
  #mzWeight & intWeight defaults = MONA massbank settings
  
  #Creating the output data frame
  alignPlot<-list()
  
  #Cleaning up the spectra, removing anything below or equal to the precursor m/z and adding mzPrec to the alignPlot obj in order to plot the precursor
  if(precFilter==TRUE){
    mzPrecDiff<-(dPPM/10^6)*mzPrec
    ms2SampSpec[[1]]<-ms2SampSpec[[1]][which(ms2SampSpec[[1]]$mz<=(mzPrec+(mzPrecDiff*2))),]
    
    #Checking that any fragments left, otherwise returns dummy-DF with cosine=0
    if(nrow(ms2SampSpec[[1]][1])==0){
      alignPlot$simScore<-(-1)
      return(alignPlot)
    }
    
    ms2StdSpec<-ms2StdSpec[which(ms2StdSpec$mz<=(mzPrec+(mzPrecDiff*2))),]
    # alignPlot$sampMZPrec<-mzPrec
  }
  
  alignPlot$sampMZPrec<-mzPrec
  
  #If ms2StdMetaData provided then add that information for plotting
  if(!is.null(ms2StdMetaData)){
    alignPlot$ms2StdMetaData<-ms2StdMetaData
  } else {
    alignPlot$ms2StdMetaData<-""
  }
  
  #If user has specified absThreshold, anything below that fraction will be removed
  if(!is.null(absThreshold) && absThreshold > 0){
    # ms2SampSpec[[1]]$int<-ms2SampSpec[[1]]$int/max(ms2SampSpec[[1]]$int)
    # ms2StdSpec$int<-ms2StdSpec$int/max(ms2StdSpec$int)
    if(all(ms2SampSpec[[1]]$int<absThreshold)){
      alignPlot$simScore <- (-1)
      return(alignPlot)
    } else if(any(ms2SampSpec[[1]]$int<absThreshold)){
      ms2SampSpec[[1]]<-ms2SampSpec[[1]][-which(ms2SampSpec[[1]]$int<absThreshold),]
    }
    #### For LEVEL 1 the noise should have been removed prior to saving the spectrum ####
    #### Can't have standard spectrum intensity for LEVEL2 since we don't know the noise threshold values ####
    # if(any(ms2StdSpec$int < absThreshold)){
    #   ms2StdSpec <- ms2StdSpec[-which(ms2StdSpec$int < absThreshold), ]
    # }
    # if(any(ms2StdSpec$int<percentFilter)){
    #   ms2StdSpec<-ms2StdSpec[-which(ms2StdSpec$int<percentFilter),]
    # }
  }
  
  #Checking that any fragments left, otherwise returns dummy-DF with cosine=0
  if(nrow(ms2SampSpec[[1]][1])==0){
    alignPlot$simScore <- (-1)
    return(alignPlot)
  }
  
  #Loop to find which fragments can be matched based on dPPM
  for(i in 1:length(ms2SampSpec[[1]]$mz)){
    mzDiff<-(dPPM/10^6)*ms2SampSpec[[1]]$mz[i]
    if(any(ms2StdSpec$mz < ms2SampSpec[[1]]$mz[i]+mzDiff & ms2StdSpec$mz > ms2SampSpec[[1]]$mz[i]-mzDiff)){
      ms2StdSpec$mz[which(ms2StdSpec$mz < ms2SampSpec[[1]]$mz[i]+mzDiff & ms2StdSpec$mz > ms2SampSpec[[1]]$mz[i]-mzDiff)]<-ms2SampSpec[[1]]$mz[i]
      alignPlot$matched<-c(alignPlot$matched,as.double(ms2SampSpec[[1]]$mz[i]))
    }
  }
  
  #Merging the spectra based on matched up masses
  alignDF <- merge(ms2SampSpec[[1]], ms2StdSpec, by=1, all = TRUE)
  colnames(alignDF)<-c("mz","intSamp","intStd")
  alignDF[is.na(alignDF)]<-0
  
  
  #Weighting mz and int based on user input
  alignDF[,4]<-alignDF[,1]^mzWeight * alignDF[,2]^intWeight
  alignDF[,5]<-alignDF[,1]^mzWeight * alignDF[,3]^intWeight
  colnames(alignDF)[4:5]<-c("sampWt","stdWt")
  
  # #Testing removing non-matched masses
  # if(any(alignDF[,2]==0)){
  #   alignDF<-alignDF[-which(alignDF[,2]==0),]
  # }
  # if(any(alignDF[,3]==0)){
  #   alignDF<-alignDF[-which(alignDF[,3]==0),]
  # }
  
  #Stopping comparison if no matching peaks
  if(dim(alignDF)[1]==0){
    alignPlot$simScore<-0.00001
    return(alignPlot)
  }
  
  
  #Merging fragments with identical mass together and taking their mean value to symbolize their combined intensity
  #Mean value might be the wrong approach?
  # if(mergeFrags==TRUE){
  #   tempAlignDF<-matrix(data=NA,nrow=length(unique(alignDF[,1])),ncol=5)
  #   colnames(tempAlignDF)<-c("mz","intSamp","intStd","sampWt","stdWt")
  #   for(i in 1:length(unique(alignDF[,1]))){
  #     tempAlignDF[i,]=colMeans(alignDF[which(alignDF[,1]==unique(alignDF[,1])[i]),])
  #   }
  #   alignDF<-as.data.frame(tempAlignDF)
  # }
  
  #Reducing matched to the actual fragments used
  alignPlot$matched<-alignPlot$matched[which(alignPlot$matched %in% alignDF$mz)]
  
  #Performing similarity scoring using 
  simScore <- as.vector((alignDF$sampWt %*% alignDF$stdWt)/(sqrt(sum(alignDF$sampWt^2)) * sqrt(sum(alignDF$stdWt^2))))
  
  print(ms2StdMetaData)
  
  #Head-to-tail plot of spectra
  alignPlot$mzSamp<-ms2SampSpec[[1]]$mz
  alignPlot$mzStd<-ms2StdSpec$mz
  alignPlot$highestIntSamp<-max(ms2SampSpec[[1]]$int)
  alignPlot$highestIntStd<-max(ms2StdSpec$int)
  alignPlot$intSamp<-ms2SampSpec[[1]]$int/max(ms2SampSpec[[1]]$int)*100
  alignPlot$intStd<-ms2StdSpec$int/max(ms2StdSpec$int)*100
  alignPlot$simScore<-simScore
  alignPlot$rtSamp<-ms2SampSpec$rt
  alignPlot$rtStd<-ms2StdMetaData$rt
  
  return(alignPlot)
}