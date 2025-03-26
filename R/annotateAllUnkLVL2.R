###annotateAllUnkMONA####
#' AnnotateAllUnkMONA - A function for going through all the matches for masses of interest in a pair of .mzML & .mgf, supplied as a "matchMSUnk" function output
#'
#' @param combListUnkObj A "combListUnkObj" object to be used for the analysis, containing MS2s which were matched using "matchMSUnk"
#' @param metaData MetaData formatted according to example file (see matchMSfunction() & getChemForm())
#' @param dbName Name / path to a SQLite of the same format as the example file. To make a new one see function "buildDB()"
#' @param chromPol
#' @param dPPM
#' @param precFilter
#' @param cosineCutOff
#' @param MONA An integer defining whether no MONA check should be done (0), MONA should be done for compounds with no match at lvl 1 (1) or for all compounds (2)
#'
#' @return A list of mgf-objects resulting from MS1 vs Theoretical dPPM match + Precursor vs Theoretical dPPM match + MS1 vs Precursor dPPM*2 match + MS1 RT vs MS2 RT within rtDiff
#' @export
#'


annotateAllUnkLVL2<-function(metaData,
                             adduct,
                             MSFilePath,
                             chromPol,
                             dbMona,
                             dPPM=10,
                             dRT=NULL,
                             precFilter=TRUE,
                             absThreshold=NULL,
                             cosineCutOff=0.8,
                             dirName="matchMSReports",
                             mzWeight=2,
                             intWeight=0.5,
                             minPeakMatch=3,
                             useCE=F,
                             useRT=F){
  
  # massOI<-as.numeric(metaData[,which(colnames(metaData)==adduct)])
  MS2FileList<-metaData$MS2File
  # if(!dir.exists(dirName))
  # suppressWarnings(dir.create(dirName))
  
  if(any(is.na(MS2FileList))){
    MS2FileList<-MS2FileList[-which(is.na(MS2FileList))]
  }
  
  ###Loop for MS2 spectra to be added to DB; Peak in MS1? RT from MS1 matching RT from MS2? RT from 
  importList<-list()
  MS1list<-list()
  MS2File<-character()
  wiwList<-data.frame()
  
  # print(length(MS2FileList))
  for(i in 1:length(MS2FileList)){
    
    #Checking if still on same MS2 (.mgf) file, if not it's loading new one into memory
    MS2File<-paste0(MSFilePath,"/",MS2FileList[i])
    mgf<-readMGF(MS2File)
    
    if(length(mgf)==0){
      next
    }
    
    n<-1
    tempSpectra<-list()
    tempMS1mz<-vector()
    tempRT1<-vector()
    tempMS1INT<-vector()
    
    for(j in 1:length(mgf)){
      preMZ2<-mgf[[j]]@precursorMz
      preRT2<-mgf[[j]]@rt
      preINT2<-mgf[[j]]@precursorIntensity
      
      tempSpectra[[n]]<-mgf[[j]]
      tempMS1mz[n]<-preMZ2
      tempRT1[n]<-preRT2
      tempMS1INT[n]<-preINT2
      
      n<-n+1
    }
    
    importList[[i]]<-tempSpectra
    MS1list[[i]]<-data.frame(tempMS1mz,tempRT1, tempMS1INT)
    wiwList<-rbind(wiwList,cbind(NA,NA,metaData$MS2Files[i]))
    
  }
  
  combListUnkObj<-list("cleanL"=importList, "ms1Info"=MS1list, "adduct"=adduct, "whichIsWhere"=wiwList) #"chromPlots"=chromPlots
  
  progressBar$inc(1/3)
  progressBar$set(message="Matching against DB")
  
  ##############################
  ####Cosine sim of all MS2s####
  
  #Creating a list in which to collect all outputs from cosineSim()
  alignPlotSamplesList<-list()
  
  #Checking which level of MONA identification the user is looking for
  MONAmatches<-list()
  mMatches<-0
  totalMONAspectra<-0
  
  # print(combListUnkObj$cleanL)
  
  #Loop through every _file-pair_ which had MS2s for masses of interest
  for(i in 1:length(combListUnkObj$cleanL)){
    
    #Mainting list to match input file
    alignPlotSamplesList[[i]]<-list()
    
    if(length(combListUnkObj$cleanL[[i]])<1){
      progressBar$inc(((1/length(unique(combListUnkObj$ms1Info[[i]][,1]))) * (1/3)))
      print("NA here inside of if-statement?")
      print(((1/length(unique(combListUnkObj$ms1Info[[i]][,1]))) * (1/3)))
      next
    }
    
    #Loop through every _precursor mz_ which had a match in the current file-pair
    for(j in 1:length(unique(combListUnkObj$ms1Info[[i]][,1]))){
      
      #Mainting list to match input file
      alignPlotSamplesList[[i]][[j]]<-list()
      
      #Fetching mzRT for precursor of current MS2match
      mzPrec<-unique(combListUnkObj$ms1Info[[i]][,1])[j]
      
      tempMona<-queryLVL2DB(massOI=mzPrec, dPPM=dPPM, polarity=chromPol, dbName=dbMona)
      
      if(typeof(tempMona)=="character"){
        alignPlotSamplesList[[i]][[j]]<-"No match found in MONA"
        progressBar$inc((1/length(unique(combListUnkObj$ms1Info[[i]][,1]))) * (1/3))
        next
      }
      
      toCheck<-which(combListUnkObj$ms1Info[[i]]$tempMS1mz==mzPrec)
      
      MONAmatches<-list()
      
      for(k in 1:length(toCheck)){
        #Grabbing the sample MS2 metadata & spectrum
        ms2SampFullFormat<-combListUnkObj$cleanL[[i]][[toCheck[k]]]
        ms2Samp<-list(data.frame(mz=ms2SampFullFormat@mz,int=ms2SampFullFormat@intensity))
        ms2Samp$rt=ms2SampFullFormat@rt
        
        switch(as.character(k%%3),
               "1" = {ms2Samp$CE <- 10},
               "2" = {ms2Samp$CE <- 20},
               "0" = {ms2Samp$CE <- 40})
        
        for(m in 1:length(tempMona$spectra)){
          #Checking if CE is matching with library and skipping if not
          if(useCE && !is.null(ms2Samp$CE) && !is.null(tempMona$metaData$ce[m])){
            if(ms2Samp$CE != tempMona$metaData$ce[m]){
              next
            }
          } else if (useCE && (is.null(tempMona$metaData$ce[m]) || tempMona$metaData$ce[m]=="NULL")){
            next
          }
          
          #Checking if RT should be checked and if it is, non-within RTs will be discarded
          if(useRT && !is.null(ms2Samp$rt) && !is.null(tempMona$metaData$rt[m]) && tempMona$metaData$rt[m] != "NULL"){
            print(tempMona$metaData$rt[m])
            
            if(abs(ms2Samp$rt - as.double(tempMona$metaData$rt[m])) > 30){
              next
            }
          } else if(useRT && (is.null(tempMona$metaData$rt[m]) || tempMona$metaData$rt[m] == "NULL")){
            next
          }
          
          
          #cosineSim between MONA spectrum and the current spectrum
          MONAmatches[[m]]<-cosineSim(ms2Samp,
                                      tempMona$spectra[[m]],
                                      dPPM=dPPM,
                                      mzWeight=mzWeight,
                                      intWeight=intWeight,
                                      mzPrec=mzPrec,
                                      ms2StdMetaData = tempMona[[2]][m,],
                                      absThreshold=absThreshold,
                                      precFilter=precFilter)
          
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
        
        alignPlotSamplesList[[i]][[j]][[k]]<-MONAmatches
      }
    }
    
    progressBar$inc((1/length(unique(combListUnkObj$ms1Info[[i]][,1]))*(1/3)))
  }
  
  progressBar$set(message="Printing PDF")
  
  #####Printing plots in a pdf for all the comparisons made#####
  pdf(file=paste0(dirName,"\\","MONASims.pdf"), onefile=T)
  plot.new()
  plot.window(xlim = c(0, 20), ylim = c(-10, 10))
  text(10,9,paste("MS2 cosine similarity analysis of a", chromPol,"dataset."), cex=1.3)
  text(10,7,paste("DB used:",dbMona))
  text(10,6,paste("dPPM between DB and samples:",dPPM))
  #text(10,5,paste("RT window between DB and samples:",rtWin))
  text(10,4,paste("Precursor filter used:",precFilter))
  text(10,3,paste("Cosine sim. cutoff for report:", cosineCutOff))
  text(10,2,paste("Intensity filter cutoff:", absThreshold))
  text(10,0,paste("A total of",mMatches,"/",totalMONAspectra, "comparisons survived the cutoff."))
  
  
  # paste(unique(combListUnkObj$whichIsWhere[,c(3:4)]))
  # paste(Sys.time())
  # 
  # dev.off()
  
  #Checking how many MS2s there are to track wiw
  whichWIW<-0
  
  
  #Looping through all the file-pairs which have been checked
  for(i in 1:length(alignPlotSamplesList)){
    
    if(length(alignPlotSamplesList[[i]])==0){
      next
    }
    whichWIW<-whichWIW+1
    
    #Looping through all the precursor masses which had a match in the DB
    for(j in 1:length(alignPlotSamplesList[[i]])){
      
      #Looping through all the MS2s which were matched against DB
      for(k in 1:length(alignPlotSamplesList[[i]][[j]])){
        
        #Checking if there was a match, otherwise proceeding to next
        if(typeof(alignPlotSamplesList[[i]][[j]][[k]])=="character" || length(alignPlotSamplesList[[i]][[j]][[k]]) < 1){
          next
        }
        
        #Looping through all the stds which were matched against the MS2 at hand
        for(l in 1:length(alignPlotSamplesList[[i]][[j]][[k]])){
          if (alignPlotSamplesList[[i]][[j]][[k]][[l]]$simScore<cosineCutOff || length(alignPlotSamplesList[[i]][[j]][[k]][[l]]$matched)<minPeakMatch){
            next
          }
          print(plotSim(alignPlotSamplesList[[i]][[j]][[k]][[l]],
                        ms2ID=alignPlotSamplesList[[i]][[j]][[k]][[l]]$ms2StdMetaData$id,
                        wiw=whichWIW,
                        wiwList=combListUnkObj$whichIsWhere,
                        molName=alignPlotSamplesList[[i]][[j]][[k]][[l]]$ms2StdMetaData$annotation,
                        adduct=alignPlotSamplesList[[i]][[j]][[k]][[l]]$ms2StdMetaData$adduct,
                        CE=alignPlotSamplesList[[i]][[j]][[k]][[l]]$CE)) #, plotAll=plotAll
          # writeLines(paste("Match for:",i,j,k,l, "RT:",alignPlotSamplesList[[i]][[j]][[k]][[l]]$rt))
        }  
      }
    }
    progressBar$inc((1/length(alignPlotSamplesList)) * (1/3))
  }
  dev.off()
  
  return(list(combListUnkObj, alignPlotSamplesList))
}