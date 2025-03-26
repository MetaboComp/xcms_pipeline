# queryLVL2Db ###############
#' queryLVL2Db
#' @param massOI The mass of interest to use for querying the monaDB
#' @param dPPM The allowed delta PPM difference between the massOI and the entries in the MONA DB
#' @param neutralMass A boolean variable when true compares the massOI to the neutral mass entries in the DB and when false compares it to the ionized masses in the DB
#' @param chromPol The chromatography + polarity of the analysis: "RP" = Reversed phase positive, "HN" = HILIC negative, etc.
#' @param dbName The database to look for ms2 spectra in. If not specified user will be prompted to input it. Specified as name of DB plus ".db" ("DBName.db").
#' @param removeZero A boolean variable to indicate whether zero intensity fragments should be removed or not
#'
#' @export

queryLVL2DB <- function(massOI, dPPM, polarity, dbName, removeZero=T){
  
  massDiff<-(dPPM/10^6)*massOI
  
  conn <- dbConnect(RSQLite::SQLite(),dbName)
  
  #Checking if user wants to compare to neutral mass or ionized mass
  s1 <- sprintf("SELECT m2s.mz, m2s.int, m2s.ms2ID FROM ms2spectra m2s, compounds comp WHERE comp.compID=m2s.ms2ID AND comp.precursormz BETWEEN %s AND %s AND comp.ionmode = \"%s\"",
                round((massOI-massDiff),5),
                round((massOI+massDiff),5),
                polarity)
  s2 <- sprintf("SELECT comp.compID, comp.names, comp.precursortype, comp.rt,  comp.ce, comp.precursormz FROM compounds comp WHERE comp.precursormz BETWEEN %s AND %s AND comp.ionmode = \"%s\"",
                round((massOI-massDiff),5),
                round((massOI+massDiff),5),
                polarity)
  spectraTemp <- as.data.table(dbGetQuery(conn,s1))
  metaData <- as.data.frame(dbGetQuery(conn,s2))
    
  dbDisconnect(conn)
  
  #Checking that match was made in MONADB
  if(nrow(spectraTemp)<1 || nrow(metaData)<1){
    return("No match found in MONA")
  }
  
  if(removeZero==T && any(spectraTemp[,2]==0)){
    spectraTemp=spectraTemp[-which(spectraTemp[,2]==0),]
  }
  
  #Making spectras into a list of spectras
  colnames(spectraTemp)<-c("mz", "int", "id")
  spectraTemp$mz <- as.double(spectraTemp$mz)
  spectraTemp$int <- as.integer(spectraTemp$int)
  # spectraTemp$ms2ID <- as.integer(spectraTemp$ms2ID)
  spectraIDs<-metaData$compID
  spectra<-list(vector(length=length(spectraIDs)))
  for(i in 1:length(spectraIDs)){
    spectra[[i]]<-spectraTemp[which(spectraTemp$id==spectraIDs[i]),c(1:3)]
  }
  colnames(metaData)=c("id", "annotation","adduct","rt", "ce", "precursormz")
  
  data<-list("spectra"=spectra, "metaData"=metaData)
  
  return(data)
}