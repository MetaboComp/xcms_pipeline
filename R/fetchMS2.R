###fetchMS2#####
#' fetchMS2 - A function to grab ms2spectra stored in DB based on precursor m/z, chromatography, polarity and optionally retention time
#'
#' @param precMZ An exact mass of interest for which ms2spectra should be retrieved from DB
#' @param chromPol String adhering to naming conventions for samples specifying chromatography and polarity of interest
#' @param precRT The retention time of the peak of interest. Optional.
#' @param projectID The project in which the authentic standard was analyzed. Optional.
#' @param dbName The database to look for ms2 spectra in. If not specified user will be prompted to input it. Specified as name of DB plus ".db" ("DBName.db").
#' @param dPPM The m/z window to search for precursor m/z in, specified as ppm. Default of 10 dppm.
#' @param rtWin Retention time window in which to search for the reference peak in seconds. Default of 30s.
#'
#' @return A list of two lists: precursor metadata (precResults) and the accompanying ms2 spectra (ms2spectra)
#' @export
#'

fetchMS2<-function(precMZ, chromPol, precRT="", projectID="", dbName, dPPM=10, rtWin=30){
  
  conn <- dbConnect(RSQLite::SQLite(), dbName)
  
  massDiff<-(dPPM/10^6)*precMZ
  
  # s1 <- sprintf("SELECT * FROM peaks INNER JOIN peakMS2link ON peaks.peakID=peakMS2link.peakID WHERE peaks.mz BETWEEN \"%s\" AND \"%s\"",
  #               round((precMZ-massDiff),5),
  #               round((precMZ+massDiff),5))
  s1 <- sprintf("SELECT * FROM peaks INNER JOIN peakMS2link ON peaks.peakID=peakMS2link.peakID WHERE peaks.mz BETWEEN %s AND %s",
                round((precMZ-massDiff),5),
                round((precMZ+massDiff),5))
  
  #Should be part of above statement once all the submission stuff has been properly implemented
  #INNER JOIN injections ON peaks.injID=injections.injID
  
  res<-dbSendQuery(conn, s1)
  precResults<-dbFetch(res)
  dbClearResult(res)
  precResults$RT<-as.double(precResults$RT)
  precResults$injID <- as.integer(precResults$injID)
  
  #If RT supplied, use it to filter results
  #Perhaps this should be done in the query itself to avoid filling up memory with shit
  if(precRT!=""){
    
    precResults<-precResults[which(precResults$RT<(precRT+(rtWin/2)) & precResults$RT>(ifelse(precRT-(rtWin/2)<0,0,precRT-(rtWin/2)))),]
    
    # print(precResults)
    if(dim(precResults)[1]==0){
      #writeLines(paste("No match in DB for",precMZ,"with RT:",precRT,"\n"))
      dbDisconnect(conn)
      return()
    }
  }
  
  #Code which should distinguish between entries in DB based on chrompol
  #results<-results[which(results$chromPol==chromPol),]
  
  #Fetching MS2spectra belonging to the precursor hits found
  if(dim(precResults)[1]==0){
    dbDisconnect(conn)
    return()
  } else {
    ms2spectra<-list()
    for(i in 1:length(precResults[,1])){
      s2 <- sprintf("SELECT mz, int FROM ms2spectra WHERE ms2ID=%s",
                    precResults$ms2ID[i])
      res<-dbSendQuery(conn, s2)
      ms2spectra[[i]]<-dbFetch(res)
      dbClearResult(res)
    }
    dbDisconnect(conn)
    
    ms2Obj<-list("spectra"=ms2spectra,"peakInfo"=precResults) #Add results here when chromPol added
    return(ms2Obj)
  }
}