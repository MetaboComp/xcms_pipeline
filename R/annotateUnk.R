###annotateUnk####
#' AnnotateUnk - A function for going through all the matches for masses of interest in a pair of .mzML & .mgf, supplied as a "matchMSUnk" function output
#'
#' @param combListUnkObj A "combListUnkObj" object to be used for the analysis, containing MS2s which were matched using "matchMSUnk"
#' @param chromPol
#' @param dbName Name / path to a SQLite of the same format as the example file. To make a new one see function "buildDB()"
#' @param dPPM
#' @param rtWin
#' @param precFilter
#' @param cosineCutOff
#' @param dirName
#' @param mzWeight
#' @param intWeight
#' @param minPeakMatch
#'
#' @return A list of mgf-objects resulting from MS1 vs Theoretical dPPM match + Precursor vs Theoretical dPPM match + MS1 vs Precursor dPPM*2 match + MS1 RT vs MS2 RT within rtDiff
#' @export
#'


annotateUnk<-function(combListUnkObj,
                      chromPol,
                      dbName,
                      dPPM=10,
                      rtWin=60,
                      precFilter=TRUE,
                      percentFilter=NULL,
                      absThresh=NULL,
                      cosineCutOff=0.8,
                      dirName="matchMSReports",
                      mzWeight=0,
                      intWeight=1,
                      minPeakMatch=3){
  
  #Creating a list in which to collect all outputs from cosineSim()
  alignPlotSamplesList<-list()
  
  #Variable for counting matches and total number of spectra considered
  nMatches<-0
  totalSpectra<-0
  
  #Loop through every _file-pair_ which had MS2s for masses of interest
  for(i in 1:length(combListUnkObj$cleanL)){
    
    #Mainting list to match input file
    alignPlotSamplesList[[i]]<-list()
    
    if(length(combListUnkObj$cleanL[[i]])<1){
      next
    }
    
    #Loop through every _precursor mz_ which had a match in the current file-pair
    for(j in 1:length(combListUnkObj$cleanL[[i]])){
      
      #Mainting list to match input file
      alignPlotSamplesList[[i]][[j]]<-list()
      
      #Loop going through all the RT and MZ matches for each precursor to compare against DB
      for(k in 1:length(unique(combListUnkObj$ms1Info[[i]][[j]])[,1])){
        
        alignPlotSamplesList[[i]][[j]][[k]]<-list()
        
        
        #Fetching mzRT for precursor of current MS2match
        mzRT<-unique(combListUnkObj$ms1Info[[i]][[j]])[k,c(1:2)]
        ms2Std<-suppressMessages(fetchMS2(precMZ=as.double(mzRT[1]), chromPol=chromPol, precRT=as.double(mzRT[2]), dbName=dbName, dPPM=dPPM, rtWin=rtWin))
        
        if(is.null(ms2Std)){
          alignPlotSamplesList[[i]][[j]][[k]]<-"No match in DB"
          next
        }
        
        if(precFilter==TRUE){
        mzPrec<-as.double(mzRT[1])
        } else {
          mzPrec<-NULL
        }
        
        toCheck<-which(combListUnkObj$ms1Info[[i]][[j]]$tempRT1==mzRT[1,2] & combListUnkObj$ms1Info[[i]][[j]]$tempMS1mz==mzRT[1,1])
        
        #Loop through every MS2 for each mass of interest which had a match against the DB 
        for(l in 1:length(toCheck)){
          
          #Grabbing the sample MS2 metadata & spectrum
          ms2SampFullFormat<-combListUnkObj$cleanL[[i]][[j]][[toCheck[l]]]
          ms2Samp<-list(data.frame(mz=ms2SampFullFormat@mz,int=ms2SampFullFormat@intensity))
          ms2Samp$rt=ms2SampFullFormat@rt
          
          # print("combListUnkObj")
          # print(combListUnkObj)
          # 
          # print("which(combListUnkObj$ms1Info[[i]][[j]]$tempMS1mz==mzRT[1,1])")
          # print(which(combListUnkObj$ms1Info[[i]][[j]]$tempMS1mz==mzRT[1,1]))
          # 
          # print("combListUnkObj$ms1Info[[i]][[j]]$tempMS1mz")
          # print(combListUnkObj$ms1Info[[i]][[j]]$tempMS1mz)
          # 
          # print("mzRT[1,1]")
          # print(mzRT[1,1])
          # 
          # print("toCheck")
          # print(toCheck)
          
          
          switch(as.character(l%%3), #which(which(combListUnkObj$ms1Info[[i]][[j]]$tempMS1mz==mzRT[1,1]) %in% toCheck)
                 "1" = {ms2Samp$CE <- 10},
                 "2" = {ms2Samp$CE <- 20},
                 "0" = {ms2Samp$CE <- 40})
          
          alignPlotList<-list()
          
          #Checking if MS2 was found in local DB and if it was cosineSiming all matches
          if(!is.null(ms2Std)){
            #Comparing current MS2 spectra to all spectras found for the same precursor mass in DB
            for(m in 1:length(ms2Std$spectra)){
              
              #Checking if CE is matching with library and skipping if not
              if(!is.null(ms2Samp$CE)){
                if(ms2Samp$CE != ms2Std$peakInfo$CollE[m]){
                  next
                }
              }
              
              #Performing cosine similarity check
              alignPlotList[[m]]<-cosineSim(ms2Samp,ms2Std$spectra[[m]], dPPM=dPPM, mzWeight=mzWeight, intWeight=intWeight, mzPrec=mzPrec, ms2StdMetaData = ms2Std$peakInfo[m,], precFilter=precFilter)
              
              if(!is.null(ms2Samp$CE)){
                alignPlotList[[m]]$CE <- ms2Samp$CE
              } else {
                alignPlotList[[m]]$CE <- ""
              }
              
              alignPlotList[[m]]$adduct<-ms2Std$peakInfo$adduct[m]
              alignPlotList[[m]]$stdMZPrec <- ms2Std$peakInfo$mz[m]
              
              # print(alignPlotList[[m]])
              # print(length(alignPlotList[[m]]$matched))
              # print(minPeakMatch)
              
              if(alignPlotList[[m]]$simScore>cosineCutOff && length(alignPlotList[[m]]$matched) >= minPeakMatch){
                nMatches<-nMatches+1
              }
              
              totalSpectra<-totalSpectra+1
            }
          }
          
          
          #Using the same structure as input file to save plots but with 4th layer since several stds in DB might have been matched
          alignPlotSamplesList[[i]][[j]][[k]][[l]]<-alignPlotList
          
        }
      }
    }
  }
  
  # print(alignPlotSamplesList)
  
  # [[]] == MS1&MS2 pair
  # [[]][[]] == MoI matches
  # [[]][[]][[]] == Number of spectra matches / MoI
  # [[]][[]][[]][[]] == Matches in DB for each MS2
  
  
  ##############################################################
  #####Printing plots in a pdf for all the comparisons made#####
  ##############################################################
  
  pdf(file=paste0(dirName,"\\","CosineSims.pdf"), onefile=T)
  plot.new()
  plot.window(xlim = c(0, 20), ylim = c(-10, 10))
  text(10,9,paste("MS2 cosine similarity analysis of a", chromPol,"dataset."), cex=1.3)
  text(10,7,paste("DB used:",dbName))
  text(10,6,paste("dPPM between DB and samples:",dPPM))
  text(10,5,paste("RT window between DB and samples:",rtWin))
  text(10,4,paste("Precursor filter used:",precFilter))
  text(10,3,paste("Cosine sim. cutoff for report:", cosineCutOff))
  text(10,0,paste("A total of",nMatches,"/",totalSpectra, "comparisons survived the cutoff."))
  
  
  # paste(unique(combListUnkObj$whichIsWhere[,c(3:4)]))
  # paste(Sys.time())
  # 
  # dev.off()
  
  #Checking how many MS2s there are to track wiw
  whichWIW<-0
  
  suppressWarnings(
    #Looping through all the file-pairs which have been checked
    for(i in 1:length(alignPlotSamplesList)){
      
      #Looping through all the precursor masses which had a match in the DB
      for(l in 1:length(alignPlotSamplesList[[i]])){
        whichWIW<-whichWIW+1
        
        #Looping through all the MS2s which were matched against DB
        for(n in 1:length(alignPlotSamplesList[[i]][[l]])){
          
          #Looping through all std matches from DB
          for(k in 1:length(alignPlotSamplesList[[i]][[l]][[n]])){
            
            #Checking if there was a match, otherwise proceeding to next
            if(typeof(alignPlotSamplesList[[i]][[l]][[n]])=="character" || length(alignPlotSamplesList[[i]][[l]][[n]][[k]])<1){
              next
            }
            
            #Looping through all the stds which were matched against the MS2 at hand
            for(j in 1:length(alignPlotSamplesList[[i]][[l]][[n]][[k]])){
              if (alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]]$simScore<cosineCutOff || length(alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]]$matched)<minPeakMatch){
                next
              }
              print(plotSim(alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]],
                            ms2ID=alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]]$ms2StdMetaData$ms2ID,
                            wiw=whichWIW,
                            wiwList=combListUnkObj$whichIsWhere,
                            molName=alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]]$ms2StdMetaData$annotation,
                            adduct=alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]]$adduct,
                            CE=alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]]$CE)) #, plotAll=plotAll
              writeLines(paste("Match for:",i,l,n,k,j,"RT:",alignPlotSamplesList[[i]][[l]][[n]][[k]][[j]]$rt))
            }  
          }
        }
      }
    })
  
  dev.off()
  return(alignPlotSamplesList)
}