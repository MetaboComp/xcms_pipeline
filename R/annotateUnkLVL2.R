###annotateUnkMONA####
#' AnnotateUnkMONA - A function for going through all the matches for masses of interest in a pair of .mzML & .mgf, supplied as a "matchMSUnk" function output
#'
#' @param combListUnkObj A "combListUnkObj" object to be used for the analysis, containing MS2s which were matched using "matchMSUnk"
#' @param chromPol Chromatography and polarity combined: "RP", "RN", "HP", "HN"
#' @param dbMSDialPath DB file created using the "buildMSDialDB" and "readMSPFile" functions
#' @param dPPM dPPM between precursor and peak m/z
#' @param precFilter Precursor filter, removing any m/z peaks above the precursor mass
#' @param percentFilter Filter removing any peaks with intensity below a certain percentage of the highest peak in the spectra
#' @param cosineCutoff A cut-off filter for which matches to put in the report pdf-file based on how high the cosine score is
#' @param dirName Directory where reports will be saved
#' @param mzWeight Weighting for mz match in cosine dot similarity product comparison, traditionally 2 for online DB comparisons
#' @param intWeight Weighting for intensity matching in cosine dot similarity product comparison, traditionally 0.5 for online DB comparisons
#' @param minPeakMatch A cut-off filter for which matches to put in the report pdf-file based on how many mz peaks that matched
#'
#' @return A list of mgf-objects resulting from MS1 vs Theoretical dPPM match + Precursor vs Theoretical dPPM match + MS1 vs Precursor dPPM*2 match + MS1 RT vs MS2 RT within rtDiff
#' @export
#'

annotateUnkLVL2<-function(combListUnkObj,
                          chromPol,
                          dbMSDialPath,
                          dPPM,
                          precFilter=TRUE,
                          absThreshold=NULL,
                          cosineCutOff=0.8,
                          dirName="matchMSReports",
                          mzWeight=2,
                          intWeight=0.5,
                          minPeakMatch=3,
                          useCE=F,
                          useRT=F){
  
  #Creating a list in which to collect all outputs from cosineSim()
  alignPlotSamplesList<-list()
  
  #Checking which level of MONA identification the user is looking for
  MONAmatches<-list()
  mMatches<-0
  totalMONAspectra<-0
  
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
        mzPrec<-as.double(mzRT[1])

        tempMona<-queryLVL2DB(massOI=mzPrec, dPPM=dPPM, polarity=chromPol, dbName=dbMSDialPath)
        
        if(typeof(tempMona)=="character"){
          alignPlotSamplesList[[i]][[j]][[k]]<-"No match found in MONA"
          next
        }
        
        toCheck<-which(combListUnkObj$ms1Info[[i]][[j]]$tempRT1==mzRT[1,2] & combListUnkObj$ms1Info[[i]][[j]]$tempMS1mz==mzRT[1,1])
        
        MONAmatches<-list()
        
        for(l in 1:length(toCheck)){
          for(m in 1:length(tempMona$spectra)){
            #Grabbing the sample MS2 metadata & spectrum
            ms2SampFullFormat<-combListUnkObj$cleanL[[i]][[j]][[toCheck[l]]]
            ms2Samp<-list(data.frame(mz=ms2SampFullFormat@mz,int=ms2SampFullFormat@intensity))
            ms2Samp$rt=ms2SampFullFormat@rt
            
            #cosineSim between MONA spectrum and the current spectrum
            MONAmatches[[m]]<-cosineSim(ms2SampSpec=ms2Samp, ms2StdSpec=tempMona$spectra[[m]], dPPM=dPPM, mzWeight=mzWeight, intWeight=intWeight, mzPrec=mzPrec, ms2StdMetaData = tempMona[[2]][m,], precFilter=precFilter, absThreshold=absThreshold)
            
            if(!is.null(ms2Samp$CE)){
              MONAmatches[[m]]$CE <- tempMona$metaData$ce[m]
            } else {
              MONAmatches[[m]]$CE <- ""
            }
            
            MONAmatches[[m]]$adduct<-tempMona[[2]]$adduct[m]
            MONAmatches[[m]]$stdMZPrec <- tempMona[[2]]$precursormz[m]
            
            if(MONAmatches[[m]]$simScore>cosineCutOff && length(MONAmatches[[m]]$matched) >= minPeakMatch){
              mMatches<-mMatches+1
            }
            totalMONAspectra<-totalMONAspectra+1
          } 
          
          alignPlotSamplesList[[i]][[j]][[k]][[l]]<-MONAmatches
        }
        progressBar$inc(((1/length(combListUnkObj$cleanL[[i]]))*0.5))
      }
    }
  }
  
  # print(alignPlotSamplesList)
  
  #suppressWarnings(
  #####Printing plots in a pdf for all the comparisons made#####
  pdf(file=paste0(dirName,"\\","MONASims.pdf"), onefile=T)
  plot.new()
  plot.window(xlim = c(0, 20), ylim = c(-10, 10))
  text(10,9,paste("MS2 cosine similarity analysis of a", chromPol,"dataset."), cex=1.3)
  text(10,7,paste("DB used:",dbMSDialPath))
  text(10,6,paste("dPPM between DB and samples:",dPPM))
  #text(10,5,paste("RT window between DB and samples:",rtWin))
  text(10,4,paste("Precursor filter used:",precFilter))
  text(10,3,paste("Cosine sim. cutoff for report:", cosineCutOff))
  text(10,0,paste("A total of",mMatches,"/",totalMONAspectra, "comparisons survived the cutoff."))
  
  #Checking how many MS2s there are to track wiw
  whichWIW<-0
  
  
  #Looping through all the file-pairs which have been checked
  for(i in 1:length(alignPlotSamplesList)){
    
    #Looping through all the precursor masses which had a match in the DB
    for(j in 1:length(alignPlotSamplesList[[i]])){
      whichWIW<-whichWIW+1
      
      #Looping through all the MS2s which were matched against DB
      for(k in 1:length(alignPlotSamplesList[[i]][[j]])){
        
        #Looping through all std matches from DB
        for(l in 1:length(alignPlotSamplesList[[i]][[j]][[k]])){
          
          #Checking if there was a match, otherwise proceeding to next
          # print(typeof(alignPlotSamplesList[[i]][[j]][[k]][[l]]))
          if(typeof(alignPlotSamplesList[[i]][[j]][[k]][[l]])=="character"){
            next
          }
          
          #Looping through all the stds which were matched against the MS2 at hand
          for(m in 1:length(alignPlotSamplesList[[i]][[j]][[k]][[l]])){
            # print(alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$simScore)
            # print(cosineCutOff)
            # print(length(alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$matched))
            # print(minPeakMatch)
            
            if (alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$simScore<cosineCutOff || length(alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$matched)<minPeakMatch){
              next
            }
            
            # print("plotSim")
            print(plotSim(alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]],
                          ms2ID=alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$ms2StdMetaData$id,
                          wiw=whichWIW, wiwList=combListUnkObj$whichIsWhere,
                          molName=alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$ms2StdMetaData$annotation,
                          adduct=alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$ms2StdMetaData$adduct),
                          CE=alignPlotSamplesList[[i]][[j]][[k]][[l]]$CE) #, plotAll=plotAll
            writeLines(paste("Match for:",i,j,k,l,m,"RT:",alignPlotSamplesList[[i]][[j]][[k]][[l]][[m]]$rt))
          }  
        }
      }
    }
  }

  dev.off()
  return(alignPlotSamplesList)
}
