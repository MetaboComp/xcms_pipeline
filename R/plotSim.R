###plotSim####
#' plotSim - A function for comparing 2 ms2 spectra, one from a SQLite DB and the other from user inputted data
#'
#' @param alignPlot An "alignPlot" object to be plotted
#' @param ms2ID The ms2ID of the std in the DB to be plotted
#' @param wiw "What I want" ???
#' @param wiwList ???
#' @param dPPM delta PPM used for the similarity scoring
#' @param molName The name of the molecule in the DB which will be plotted
#' @param adduct The precursor ion type to be plotted
#'
#' @return A list of mgf-objects resulting from MS1 vs Theoretical dPPM match + Precursor vs Theoretical dPPM match + MS1 vs Precursor dPPM*2 match + MS1 RT vs MS2 RT within rtDiff
#' @export
#'

plotSim<-function(alignPlot, ms2ID="", wiw=NULL, wiwList=NULL, dPPM=10, molName="", adduct="", CE=""){
  
  plot.new()
  par(mar=c(3,3,3,3), oma=c(3,3,3,3))
  xlim<-c(0,(max(alignPlot$mzSamp,alignPlot$mzStd)+10))
  #Setting up the plot window in which to plot the bars
  plot.window(xlim = xlim, ylim = c(-125, 125))
  #Axis ticks
  yTicks <- c(-100, -50, 0, 50, 100)
  yTicksRight <- c(-alignPlot$highestIntStd,
                    -alignPlot$highestIntStd/2,
                    0,
                    alignPlot$highestIntSamp/2,
                   alignPlot$highestIntSamp)
  
  maxMZ<-max(alignPlot$mzSamp)
  if(max(alignPlot$mzStd)>maxMZ){
    maxMZ<-max(alignPlot$mzStd)
  }
  
  # if(plotAll==T){
  
  #Plotting a red line for each matched fragment and a blue for each unmatched
  for (i in 1:length(alignPlot$mzSamp)){
    # if(any(alignPlot$mzSamp[i]==alignPlot$matched)){
    #   lines(rep(alignPlot$mzSamp[i], 2), c(0, alignPlot$intSamp[i]), col = "red", lwd=0.8)
    # }else{
    lines(rep(alignPlot$mzSamp[i], 2), c(0, alignPlot$intSamp[i]), col = "blue", lwd=0.8)
    # }
  } 
  #Plotting a red line for each matched fragment and a blue for each unmatched
  for (i in 1:length(alignPlot$mzStd)){
    # if(any(alignPlot$mzStd[i]==alignPlot$matched)){
    #   lines(rep(alignPlot$mzStd[i], 2), c(0, -alignPlot$intStd[i]), col = "red", lwd=0.8)
    # } else {
    lines(rep(alignPlot$mzStd[i], 2), c(0, -alignPlot$intStd[i]), col = "blue", lwd=0.8)
    # }
  }
  
  #Plotting all red lines, making sure they end up on top of the blue ones
  for (i in 1:length(alignPlot$matched)){
    lines(rep(alignPlot$matched[i], 2), c(0, alignPlot$intSamp[max(which(alignPlot$mzSamp==alignPlot$matched[i]))]), col="red", lwd=0.8)
    lines(rep(alignPlot$matched[i], 2), c(0, -alignPlot$intStd[max(which(alignPlot$mzStd==alignPlot$matched[i]))]), col="red", lwd=0.8)
  }
  
  # } else {
  #   
  #   #Plotting a red line for each matched fragment and a blue for each unmatched
  #   for (i in 1:length(alignPlot$matched)){
  #       lines(rep(alignPlot$matched[i], 2), c(0, alignPlot$intSamp[which(alignPlot$mzSamp==alignPlot$matched[i])]), col = "red", lwd=0.8)
  #   } 
  #   #Plotting a red line for each matched fragment and a blue for each unmatched
  #   for (i in 1:length(alignPlot$matched)){
  #       if(length(which(alignPlot$mzStd==alignPlot$matched[i])>1)){
  #         lines(rep(alignPlot$matched[i], 2), c(0, -alignPlot$intStd[max(which(alignPlot$mzStd==alignPlot$matched[i]))]), col = "red", lwd=0.8)
  #       } else {
  #         lines(rep(alignPlot$matched[i], 2), c(0, -alignPlot$intStd[which(alignPlot$mzStd==alignPlot$matched[i])]), col = "red", lwd=0.8)
  #       }
  #       
  #   }
  # }
  
  
  #Setting up the two parts of left y-axes
  axis(4, at = yTicks, labels = formatC(abs(yTicksRight), format="e", digits=1), pos=xlim[2], ylab="intensity")
  axis(2, at = yTicks, labels = abs(yTicks), pos = xlim[1], 
       ylab = "intensity %")
  axis(1, pos = -125, at=seq(0,maxMZ,by=25))
  
  #Setting up the dividing line between the spectra
  lines(xlim, c(0, 0), lwd=0.8)
  #Setting up the rectangle around the spectra
  rect(xlim[1], -125, xlim[2], 125)
  
  
  #Text on the axes and in the box
  mtext("m/z", side = 1, line = 2)
  mtext("intensity (%)", side = 2, line = 2)
  
  #Defining where the text goes
  plot.window(xlim = c(0, 20), ylim = c(-10, 10))
  text(10, 9, "Sample")
  text(10, -9, "Std")
  
  #Left inner margin
  mtext(paste("Sim score:",round(alignPlot$simScore,2)), side=3, line=2, cex=1, adj=0)
  mtext(paste("Sample peaks matched:",length(alignPlot$matched),"/",length(alignPlot$mzSamp)), side=3, line=1, cex=1, adj=0)
  mtext(paste("Std peaks matched:",length(alignPlot$matched),"/",length(alignPlot$mzStd)), side=3, line=0, cex=1, adj=0)
  
  #Debug
  # message(paste0("CE: ", CE))
  # message(paste0("dPPM_1: ", alignPlot$sampMZPrec))
  # message(paste0("dPPM_2: ", alignPlot$stdMZPrec))
  # message(paste0("dPPM_3: ", round((as.numeric(alignPlot$sampMZPrec)-as.numeric(alignPlot$stdMZPrec))/(as.numeric(alignPlot$stdMZPrec)/10^6), 1)))
  # message("-------------------------------")
  
  #Title generation, depending on input values
  if(!is.null(alignPlot$sampMZPrec) && !is.null(alignPlot$ms2StdMetaData) && nrow(alignPlot$ms2StdMetaData) != 0 && ms2ID != ""){
    mtext(paste("Sample mz & RT", round(alignPlot$sampMZPrec,5), "@", round(alignPlot$rtSamp,1),"s"), side=3, line=2, cex=1, adj=0, outer=TRUE)
    mtext(paste("MS2ID:", ms2ID), side=3, line=2, cex=1, adj=1)
    if(!is.null(CE)){
      mtext(paste("CE:", CE), side=3, line=1, cex=1, adj=1)
    }
    mtext(paste("dPPM", round((as.numeric(alignPlot$sampMZPrec)-as.numeric(alignPlot$stdMZPrec))/(as.numeric(alignPlot$stdMZPrec)/10^6), 1)), side=3, line=0, cex=1, adj=1)
    mtext(paste("dRT:", round(as.numeric(alignPlot$rtSamp)-as.numeric(alignPlot$rtStd),0),"s"), side=3, line=3, cex=1, adj=1)
  } else {
    mtext(paste("CosineSim plot"), side=3, line=2, cex=1, adj=0, outer=TRUE)
    mtext(paste("MS2ID:", ms2ID), side=3, line=2, cex=1, adj=1)
    mtext(paste("dPPM", round((as.numeric(alignPlot$sampMZPrec)-as.numeric(alignPlot$stdMZPrec))/(as.numeric(alignPlot$stdMZPrec)/10^6), 1)), side=3, line=0, cex=1, adj=1)
    mtext(paste("dRT:", round(as.numeric(alignPlot$rtSamp)-as.numeric(alignPlot$rtStd),0),"s"), side=3, line=3, cex=1, adj=1)
  }
  #Right inner margin
  if(molName!=""){
    mtext(paste("MolName:", molName), side=3, line=2, cex=1, adj=1, outer=TRUE)
  }
  
  if(!is.na(adduct)){
    if(adduct!=""){
      mtext(paste("Adduct:",adduct),side=3, line=1, cex=1, adj=1, outer=TRUE)
    }
  }
  
  # print(wiw)
  # print(wiwList)
  if(!is.null(wiw) && !is.null(wiwList)){
    mzPrecDiff<-(dPPM/10^6)*alignPlot$sampMZPrec
    plotFileName<-wiwList[wiw,3] #which(wiwList$MS2MassMatch > alignPlot$sampMZPrec-mzPrecDiff & wiwList$MS2MassMatch < alignPlot$sampMZPrec+mzPrecDiff & wiwList$ReportFile==wiw),3]
    # if(length(plotFileName)<1){
    #   mzPrecDiff<-((dPPM*2)/10^6)*alignPlot$sampMZPrec
    #   plotFileName<-wiwList[which(wiwList$MS2MassMatch > alignPlot$sampMZPrec-mzPrecDiff & wiwList$MS2MassMatch < alignPlot$sampMZPrec+mzPrecDiff & wiwList$ReportFile==wiw),3]
    #   if(length(plotFileName)<1){
    #     return()
    #   }
    # }
    plotFileName<-strsplit(plotFileName, "\\.")[[1]][1]
    mtext(paste("Sample:", plotFileName), side=3, line=1, cex=1, adj=0, outer=TRUE)
  }
  
  #Exporting plot for further modification in other function or printing
  exportPlot<-recordPlot()
  return(exportPlot)
}