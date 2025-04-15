#### plotClust
# Makes pdf files with merged intensities from ramClust objects from 12 different settings to be compared to identify the best settings for ramClust grouping

plotClust=function(ram,clustnr,xcmsData,samps,dtime=15,dmz=.05) {
  if(missing(samps)) {
    nSamp=nrow(ram$SpecAbund)
    samps=1:nSamp
  } else nSamp=length(samps)
  whichFeats=which(ram$featclus==clustnr)
  peakMeta=cbind(ram$fmz,ram$frt)
  pkMetaGrp=peakMeta[whichFeats,]
  rtr=ram$clrt[clustnr]+c(-dtime,dtime)
  rtr[rtr<0]=0
  mzr=cbind(ram$fmz[whichFeats]-dmz,ram$fmz[whichFeats]+dmz)
  chr <- chromatogram(xcmsData, mz = mzr, rt = rtr)
  plot(0:1,0:1,type='n',axes=F,xlab='Retention time (s)', ylab='Intensity (AU)',main=paste0('RAM cluster ',clustnr,'; RT ',signif(ram$clrt[clustnr],5),'s'))
  box(bty='l')
  for (pk in 1:length(whichFeats)) {
    rts=ints=list()
    for (samp in 1:nSamp) {
      # rts[[samp]]=chr[pk,samps[samp]]@rtime
      # ints[[samp]]=chr[pk,samps[samp]]@intensity
      rts[[samp]]=chr[pk,samp]@rtime
      ints[[samp]]=chr[pk,samp]@intensity
    }
    nrts=min(sapply(rts,length))
    rts=sapply(rts,function(x) x[1:nrts])
    rts=rowMeans(rts)
    ints=sapply(ints,function(x) x[1:nrts])
    ints=rowMeans(ints,na.rm=T)
    par(new=T)
    plot(rts,ints,type='l',col=pk+1,ylim=c(0,max(ints,na.rm=T)),axes=F,xlab='',ylab='')
  }
  axis(1)
  legend('topright',legend = paste0('F',whichFeats,'@mz',signif(pkMetaGrp[,1],5)), lty=1,col=(1:length(whichFeats))+1,bty='n')
}