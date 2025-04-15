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


annotateAllUnk<-function(metaData,
                         adduct,
                         MSFilePath,
                         chromPol,
                         dbName,
                         dPPM=10,
                         rtWin=30,
                         precFilter=TRUE,
                         absThreshold=NULL,
                         cosineCutOff=0.8,
                         dirName="matchMSReports",
                         mzWeight=0,
                         intWeight=1,
                         minPeakMatch=3){
  
  nMatches<-0
  totalSpectra<-0
  
  # massOI<-as.numeric(metaData[,which(colnames(metaData)==adduct)])
  MS2FileList<-metaData$MS2File
  # if(!dir.exists(dirName))
  # suppressWarnings(dir.create(dirName))
  
  # print(MS2FileList)
  
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
    MS1list[[i]]<-data.frame(round(tempMS1mz,5),round(tempRT1,1), round(tempMS1INT,0))
    colnames(MS1list[[i]]) <- c("tempMS1mz", "tempRT1", "tempMS1INT")
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
  
  # print(combListUnkObj$ms1Info)
  
  #Loop through every _file-pair_ which had MS2s for masses of interest
  for(i in 1:length(combListUnkObj$cleanL)){
    # print(i)
    
    #Mainting list to match input file
    alignPlotSamplesList[[i]]<-list()
    
    if(length(combListUnkObj$cleanL[[i]])<1){
      # progressBar$inc(((1/length(unique(combListUnkObj$ms1Info[[i]][,1]))) * (1/3)))
      # print("NA here inside of if-statement?")
      # print(((1/length(unique(combListUnkObj$ms1Info[[i]][,1]))) * (1/3)))
      next
    }
    
    #Loop through every _precursor mz_ which had a match in the current file-pair
    for(j in 1:length(combListUnkObj$cleanL[[i]])){
      # print(j)
      
      #Mainting list to match input file
      alignPlotSamplesList[[i]][[j]]<-list()
      
      #Loop going through all the RT and MZ matches for each precursor to compare against DB
      # for(k in 1:length(unique(combListUnkObj$ms1Info[[i]])[,1])){
      
      # alignPlotSamplesList[[i]][[j]][[k]]<-list()
      
      
      #Fetching mzRT for precursor of current MS2match
      # print(unique(combListUnkObj$ms1Info[[i]])[j,c(1:2)])
      mzRT<-unique(combListUnkObj$ms1Info[[i]])[j,c(1:2)]
      # print(chromPol)
      # print(mzRT)
      # print(dbName)
      # print(dPPM)
      # print(rtWin)
      # print(suppressMessages(fetchMS2(precMZ=as.double(mzRT[1]), chromPol=chromPol, precRT=as.double(mzRT[2]), dbName=dbName, dPPM=dPPM, rtWin=rtWin)))
      ms2Std<-suppressMessages(fetchMS2(precMZ=as.double(mzRT[1]), chromPol=chromPol, precRT=as.double(mzRT[2]), dbName=dbName, dPPM=dPPM, rtWin=rtWin))
      
      # print("After ms2fetch")
      # print(ms2Std)
      
      if(is.null(ms2Std)){
        alignPlotSamplesList[[i]][[j]]<-"No match in DB"
        next
      } else {
        # print(ms2Std$peakInfo)
      }
      # else {
      #   print("Returning")
      #   return()
      # }
      
      if(precFilter==TRUE){
        mzPrec<-as.double(mzRT[1])
      } else {
        mzPrec<-NULL
      }
      
      # print("Before tochck")
      
      toCheck<-which(combListUnkObj$ms1Info[[i]]$tempRT1==mzRT[1,2] & combListUnkObj$ms1Info[[i]]$tempMS1mz==mzRT[1,1])
      
      # print(toCheck)
      
      #Loop through every MS2 for each mass of interest which had a match against the DB 
      for(l in 1:length(toCheck)){
        # print(l)
        #Grabbing the sample MS2 metadata & spectrum
        ms2SampFullFormat<-combListUnkObj$cleanL[[i]][[toCheck[l]]]
        ms2Samp<-list(data.frame(mz=ms2SampFullFormat@mz,int=ms2SampFullFormat@intensity))
        ms2Samp$rt=ms2SampFullFormat@rt
        
        # print(which(which(combListUnkObj$ms1Info[[i]]$tempMS1mz==mzRT[1,1]) %in% toCheck)%%3)
        switch(as.character(which(which(combListUnkObj$ms1Info[[i]]$tempMS1mz==mzRT[1,1]) %in% toCheck)%%3),
               "1" = {ms2Samp$CE <- 10},
               "2" = {ms2Samp$CE <- 20},
               "0" = {ms2Samp$CE <- 40})
        
        # print(ms2Samp$CE)
        
        
        alignPlotList<-list()
        
        #Checking if MS2 was found in local DB and if it was cosineSiming all matches
        if(!is.null(ms2Std)){
          #Comparing current MS2 spectra to all spectras found for the same precursor mass in DB
          # print(ms2Std$peakInfo$CollE)
          # print(ms2Samp$CE)
          # whichCEToCheck <- which(ms2Std$peakInfo$CollE %in% ms2Samp$CE)
          # print(whichCEToCheck)
          
          # print(ms2Samp$CE)
          # print(which(ms2Std$peakInfo$CollE != ms2Samp$CE))
          # print(ms2Std$peakInfo)
          
          nListEntries <- 1
          for(m in (1:length(ms2Std$spectra))[-which(ms2Std$peakInfo$CollE != ms2Samp$CE)]){
            
            if(!is.null(ms2Samp$CE)){
              if(ms2Samp$CE != ms2Std$peakInfo$CollE[m]){
                next
              }
            }
            
            #Performing cosine similarity check
            alignPlotList[[nListEntries]]<-cosineSim(ms2Samp,
                                          ms2Std$spectra[[m]],
                                          dPPM=dPPM,
                                          mzWeight=mzWeight,
                                          intWeight=intWeight,
                                          mzPrec=mzPrec,
                                          ms2StdMetaData = ms2Std$peakInfo[m,],
                                          precFilter=precFilter,
                                          absThreshold=absThreshold)
            if(!is.null(ms2Samp$CE)){
              alignPlotList[[nListEntries]]$CE <- ms2Samp$CE
            } else {
              alignPlotList[[nListEntries]]$CE <- ""
            }
            
            alignPlotList[[nListEntries]]$adduct<-ms2Std$peakInfo$adduct[m]
            alignPlotList[[nListEntries]]$stdMZPrec <- ms2Std$peakInfo$mz[m]
            
            print(alignPlotList[[nListEntries]])
            
            # print(length(alignPlotList[[m]]$matched))
            # print(minPeakMatch)
            
            if(alignPlotList[[nListEntries]]$simScore>cosineCutOff && length(alignPlotList[[nListEntries]]$matched) >= minPeakMatch){
              nMatches<-nMatches+1
              
            }
            
            totalSpectra<-totalSpectra+1
            nListEntries <- nListEntries + 1
          }
        }
        
        
        #Using the same structure as input file to save plots but with 4th layer since several stds in DB might have been matched
        alignPlotSamplesList[[i]][[j]][[l]]<-alignPlotList
        
      }
      # }
    }
    
    progressBar$inc((1/length(unique(combListUnkObj$ms1Info[[i]][,1]))*(1/3)))
  }
  
  progressBar$set(message="Printing PDF")
  # close(pb)
  
  summaryDF <- data.frame(matrix(ncol=16))
  colnames(summaryDF) <- c("peakID", "injID", "mz", "RT", "int", "annotation", "adduct", "IDLVL", "peakID", "ms2ID", "precMZ", "precRT", "precINT", "CollE", "injID", "curator")
  
  #####Printing plots in a pdf for all the comparisons made#####
  pdf(file=paste0(dirName,"\\","CosineSims.pdf"), onefile=T)
  plot.new()
  plot.window(xlim = c(0, 20), ylim = c(-10, 10))
  text(10,9,paste("MS2 cosine similarity analysis of a", chromPol,"dataset."), cex=1.3)
  text(10,7,paste("DB used:",dbName))
  text(10,6,paste("dPPM between DB and samples:",dPPM))
  text(10,5,paste("RT window between DB and samples:",rtWin))
  text(10,4,paste("Precursor filter used:",precFilter))
  text(10,3,paste("Cosine sim. cutoff for report:", cosineCutOff))
  text(10,2,paste("Intensity filter cutoff:", absThreshold))
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
    whichWIW<-whichWIW+1
    
      #Looping through all the precursor masses which had a match in the DB
      for(l in 1:length(alignPlotSamplesList[[i]])){
        
        
        #Looping through all the MS2s which were matched against DB
        for(n in 1:length(alignPlotSamplesList[[i]][[l]])){
          
          # Looping through all std matches from DB
          # for(k in 1:length(alignPlotSamplesList[[i]][[l]][[n]])){
            
            #Checking if there was a match, otherwise proceeding to next
            if(typeof(alignPlotSamplesList[[i]][[l]][[n]])=="character" || length(alignPlotSamplesList[[i]][[l]][[n]])<1){
              next
            }
            
            #Looping through all the stds which were matched against the MS2 at hand
            for(j in 1:length(alignPlotSamplesList[[i]][[l]][[n]])){
              if (alignPlotSamplesList[[i]][[l]][[n]][[j]]$simScore<cosineCutOff || length(alignPlotSamplesList[[i]][[l]][[n]][[j]]$matched)<minPeakMatch){
                next
              }
              # message(l)
              # message(whichWIW)
              print(plotSim(alignPlotSamplesList[[i]][[l]][[n]][[j]],
                            ms2ID=alignPlotSamplesList[[i]][[l]][[n]][[j]]$ms2StdMetaData$ms2ID,
                            wiw=whichWIW,
                            wiwList=combListUnkObj$whichIsWhere,
                            molName=alignPlotSamplesList[[i]][[l]][[n]][[j]]$ms2StdMetaData$annotation,
                            adduct=alignPlotSamplesList[[i]][[l]][[n]][[j]]$adduct,
                            CE=alignPlotSamplesList[[i]][[l]][[n]][[j]]$CE)) #, plotAll=plotAll
              writeLines(paste("Match for:",i,l,n,j,"RT:",alignPlotSamplesList[[i]][[l]][[n]][[j]]$rt))
              summaryDF <- rbind(summaryDF,
                                 alignPlotSamplesList[[i]][[l]][[n]][[j]]$ms2StdMetaData)
            }  
          }
        # }
      }
    })
  
  dev.off()
  
  # saveRDS(combListUnkObj, paste0(dirName,"/combListUnkObj_allLVL1.rds"))
  
  return(list(combListUnkObj, alignPlotSamplesList))
  
}