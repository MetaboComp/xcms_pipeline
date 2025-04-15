###########################################################################
# Post MSDial processing - drift correction, clustering of features etc.  #
###########################################################################

##### 250212 - Anton Ribbenstedt - v1.0 ####

# Loading libraries
library(openxlsx)
library(ggfortify)
library(StatTools)
library(WaveICA2.0)
library(batchCorr)
library(RAMClustR)
library(MSnbase)
library(mzR)
library(xcms)


#### Settings for script ####
dir <- "C:/Work/250211_TestProject/" #Path to directory of files

chrom <- "RP" #RP or HILIC
pol <- "NEG" #POS or NEG

instrument <- "MRT" #MRT, Synapt

#Filter removing all features which don't have MS2 data connected to them
onlyMS2FilesFilter <- F

#Filter removing features where 80th percentile is lower than specified area
#OBS Instrument dependent#
lowIntFilter <- 1000

#What this script does:
# -Removes all features with > 80% NAs (0 intensity or lower)
# -Removes all blanks from the final PT
# -






#### Reading original format MSDIAL export ####
dateToday <- substr(format(as.Date(Sys.Date()), "%Y%m%d"), 3, 1000)
chromPol <- paste0(chrom,"_",pol)
setwd(dir)
source.list<-list.files("./R/", full.names = T) #241028 Hackathon: Put in CMSITools - librarification
sapply(source.list, source)

#Loading all standalone functions


xlsxFile <- list.files(dir, pattern=".xlsx")

MSDIAL_Import <- read.xlsx(xlsxFile, sheet=1)


#### Removing all features without MS2 data ####
if(onlyMS2FilesFilter){
  MSDIAL_Import <- MSDIAL_Import[which(MSDIAL_Import[c(5:nrow(MSDIAL_Import)),8] ==T), ]
}

MSDIAL_DF <- MSDIAL_Import[c(5:nrow(MSDIAL_Import)), c(1:3)]
MSDIAL_DF <- cbind(MSDIAL_DF,
                   MSDIAL_Import[c(5:nrow(MSDIAL_Import)), c(36:(ncol(MSDIAL_Import)-2))])

colnames(MSDIAL_DF) <- c("ID", "rt", "mz",
                         MSDIAL_Import[4, c(36:(ncol(MSDIAL_Import)-2))])

if(!("Reworked" %in% getSheetNames(xlsxFile))){
  wb <- openxlsx::loadWorkbook(xlsxFile)
  addWorksheet(wb, sheetName="Reworked")
  writeData(wb, "Reworked", MSDIAL_DF)
  saveWorkbook(wb, xlsxFile, overwrite=T)
}


MSDIAL_PT <- MSDIAL_DF[,-c(1:3)]

for(i in 1:ncol(MSDIAL_PT)){
  MSDIAL_PT[,i] <- as.double(MSDIAL_PT[,i])
}

#### Arranging based on injOrder ####
injOrder <- c()
for(i in 1:ncol(MSDIAL_PT)){
  injOrder <- c(injOrder,
                as.integer(strsplit(colnames(MSDIAL_PT)[i], "_")[[1]][6]))
}
MSDIAL_PT <- MSDIAL_PT[,order(injOrder)]
injOrder <- injOrder[order(injOrder)]

#### Removing blanks ####
injOrder <- injOrder[-c(which(grepl("blank", colnames(MSDIAL_PT), ignore.case = T)),
                        which(grepl("MS2", colnames(MSDIAL_PT), ignore.case = T)),
                        which(grepl("meoh", colnames(MSDIAL_PT), ignore.case = T)),
                        which(grepl("mse", colnames(MSDIAL_PT), ignore.case = T)))]
# injOrder <- injOrder[-which(grepl("MS2", colnames(MSDIAL_PT), ignore.case = T))]
MSDIAL_PT <- MSDIAL_PT[,-c(which(grepl("blank", colnames(MSDIAL_PT), ignore.case = T)),
                           which(grepl("MS2", colnames(MSDIAL_PT), ignore.case = T)),
                           which(grepl("meoh", colnames(MSDIAL_PT), ignore.case = T)),
                           which(grepl("mse", colnames(MSDIAL_PT), ignore.case = T)))]

#### Setting rownames to alignment ID from MSDIAL & setting up meta-dfs####
rownames(MSDIAL_PT) <- MSDIAL_DF$ID
meta_feats <- cbind(MSDIAL_Import[5:nrow(MSDIAL_Import),1:3])
for(i in 1:ncol(meta_feats)){
  meta_feats[,i] <- as.double(meta_feats[,i])
}
colnames(meta_feats) <- c("ID", "rt", "mz")
meta_samps <- as.data.frame(cbind(injOrder, colnames(MSDIAL_PT), NA))
meta_samps[which(grepl("sQC", meta_samps[,2])),3] <- "sQC"
meta_samps[which(grepl("ltQC", meta_samps[,2])),3] <- "ltQC"
meta_samps[which(is.na(meta_samps[,3])),3] <- "sample"
meta_samps[,1] <- as.integer(meta_samps[,1])

# if(metaFilter == T){
#   meta_samps <- cbind(meta_samps, rep(NA, nrow(meta_samps)))
#   
#   for(i in 1:nrow(metaGroupInfo)){
#     meta_samps[grepl(metaGroupInfo[i,2], meta_samps[,2]),
#                4] <- metaGroupInfo[i,4]
#   }
#   
#   meta_samps[is.na(meta_samps[,4]),4] <- "QC"
# }


colnames(meta_samps) <- c("injOrder", "sampName", "group")


MSDIAL_PT <- t(MSDIAL_PT)

#### Filtering with too many NAs ####
MSDIAL_PT[MSDIAL_PT==0] <- NA

#Removing features with > 80% NA
whichFeatToRemove <- c()
for(i in 1:ncol(MSDIAL_PT)){
  if(sum(is.na(MSDIAL_PT[,i]))/nrow(MSDIAL_PT) > 0.8){
    whichFeatToRemove <- c(whichFeatToRemove,
                           i)
  }
}

stat_NAs_remove <- length(whichFeatToRemove)/ncol(MSDIAL_PT) * 100 #1.1% removed
MSDIAL_PT <- MSDIAL_PT[,-whichFeatToRemove]

#### Filtering with few samples with high peak and rest with very low peaks ####
MSDIAL_PT[is.na(MSDIAL_PT)] <- 0

whichFeatToRemove <- c()
#For each feature
for(i in 1:ncol(MSDIAL_PT)){
  #For each group
  if(quantile(MSDIAL_PT[, i],0.8) < lowIntFilter){ #Double check blank filter from linux-pipeline to evaluate which peaks to keep / remove
    whichFeatToRemove <- c(whichFeatToRemove,
                           i)
  }
}

stat_lowArea_remove <- length(whichFeatToRemove)/ncol(MSDIAL_PT)
MSDIAL_PT <- MSDIAL_PT[,-whichFeatToRemove]


MSDIAL_PT[MSDIAL_PT==0] <- NA


#### Imputation of features ####
#Missingness
sum(is.na(MSDIAL_PT))/(dim(MSDIAL_PT)[1] * dim(MSDIAL_PT)[2]) #Missingness: 0.5%
MSDIAL_PT <- mvImpWrap(MAT = MSDIAL_PT , method = "RF", guess = "minVar", forceZero = T)

saveRDS(MSDIAL_PT, paste0(dateToday, "_MSDIAL_PT_PostmvImp_Pos.rds"))

#### Making a PCA of ####
#PCA of pre-drift corr
pca_res <- prcomp(MSDIAL_PT, scale. = TRUE)
png(filename=paste0(dateToday, "_", chromPol,"_PreDriftCorr-Dots.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", size=3)  + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow")) #autoplot(pca_res, data=meta, colour="Type", size=3, label=T, label.size=3)
dev.off()
png(filename=paste0(dateToday, "_", chromPol,"_PreDriftCorr-Labels.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", shape=F, label=T, label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()
rownames(pca_res$x) <- meta_samps[,1]
png(filename=paste0(dateToday, "_", chromPol,"_PreDriftCorr-Numbers.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", shape=F, label=T, label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()

### WaveICA 2.0 ###
MSDIAL_PT_PostWaveICA <- WaveICA_2.0(data=MSDIAL_PT, Injection_Order=meta_samps$injOrder, alpha=0, Cutoff=0.1, K=10)
MSDIAL_PT_PostWaveICA <- as.data.frame(as.matrix(MSDIAL_PT_PostWaveICA[[1]]))

#PCA of post-drift corr
pca_res <- prcomp(MSDIAL_PT_PostWaveICA, scale. = TRUE)
png(filename=paste0(dateToday, "_", chromPol,"_PostWaveICA2-Dots.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", size=3)  + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()
png(filename=paste0(dateToday, "_", chromPol,"_PostWaveICA2-Labels.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", shape=F, label=T, label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()
rownames(pca_res$x) <- meta_samps[,1]
png(filename=paste0(dateToday, "_", chromPol,"_PostWaveICA2-Numbers.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", shape=F, label=T, label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()

### batchCorr ####
dir.create(paste0(getwd(), "/batchCorr/"))
drift_correction <- correctDrift(peakTable = MSDIAL_PT,
                                 injections = meta_samps$injOrder,
                                 sampleGroups = meta_samps$group,
                                 QCID = "sQC", RefID = "ltQC",
                                 G = seq(5,35, by = 3),
                                 reportPath = paste0(getwd(), "/batchCorr/"))

#PCA of post-drift corr
pca_res <- prcomp(drift_correction$TestFeatsFinal, scale. = TRUE)
png(filename=paste0(dateToday, "_", chromPol,"_PostBatchCorr-Dots.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", size=3)  + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()
png(filename=paste0(dateToday, "_", chromPol,"_PostBatchCorr-Labels.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", shape=F, label=T, label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()
rownames(pca_res$x) <- meta_samps[,1]
png(filename=paste0(dateToday, "_", chromPol,"_PostBatchCorr-Numbers.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", shape=F, label=T, label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()

#### Attempt to make PCA with final loadings of original data ####
pca_res <- prcomp(drift_correction$TestFeats[, which(colnames(drift_correction$TestFeats) %in% colnames(drift_correction$TestFeatsFinal))], scale.=TRUE) #rownames(pca_res$rotation)
png(filename=paste0(dateToday, "_", chromPol,"_PreDriftCorr-OnlyBatchCorrFeats-Dots.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", size=3)  + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()

### RamClustR ####
colnames(MSDIAL_PT_PostWaveICA) <- paste0(meta_feats[which(meta_feats$ID %in% colnames(MSDIAL_PT_PostWaveICA)), 3],
                                          "_",
                                          round((meta_feats[which(meta_feats$ID %in% colnames(MSDIAL_PT_PostWaveICA)), 2]*60),0))




########################################################################
#### Two versions of final PT generated, one per normalization mode ####
########################################################################
MSDIAL_PT_Final_Wave <- scaleAwayFromZero(MSDIAL_PT, MSDIAL_PT_PostWaveICA)
drift_correction$TestFeatsCorr <- scaleAwayFromZero(PTBefore = drift_correction$TestFeats,
                                                    PTAfter = drift_correction$TestFeatsCorr)
MSDIAL_PT_Final_batchCorr <- drift_correction$TestFeatsCorr[ ,which(colnames(drift_correction$TestFeatsCorr) %in% colnames(drift_correction$TestFeatsFinal))]

colnames(MSDIAL_PT_Final_batchCorr) <- paste0(meta_feats[which(meta_feats$ID %in% colnames(drift_correction$TestFeatsFinal)), 3],
                                              "_",
                                              round((meta_feats[which(meta_feats$ID %in% colnames(drift_correction$TestFeatsFinal)), 2]*60),0))

pca_res <- prcomp(MSDIAL_PT_PostWaveICA[, which(colnames(MSDIAL_PT_PostWaveICA) %in% colnames(MSDIAL_PT_Final_batchCorr))], scale.=TRUE) #rownames(pca_res$rotation)
png(filename=paste0(dateToday, "_", chromPol,"_PostWaveICA2-Dots.png"), width=2559, height=1376, res=120)
autoplot(pca_res, data=meta_samps, colour="group", size=3)  + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))
dev.off()

write.csv(as.data.frame(MSDIAL_PT_Final_Wave), file='peakTable_wave.csv')
write.csv(as.data.frame(MSDIAL_PT_Final_batchCorr), file='peakTable_batchCorr.csv')

stop("You have to check the PCA and decide if sQCs are within the spread of the samples. If they are, batchCorr should be applied.")

#### User chooses between WaveICA2.0 and batchCorr for the final peak table (PT) ####
final_PT_name <- "peakTable_wave.csv" #Either 'peakTable_wave.csv' or 'peakTable_batchCorr.csv'




#Select 5 random files to extract cluster data from in order to evaluate best settings
randommzML <- sample(list.files(, full.names = T), 5)

#Settings
expDes=defineExperiment(force.skip = T)
sr=c(.3,.4,.5)
st=c(.5,1,1.5,2)
maxt=c(5)
par=expand.grid(st=st,sr=sr,maxt=maxt)
str(par)
nClust=nSing=sizeMax=sizeMed=sizeMean=numeric(nrow(par))
nFeat=list()

files <- readMSData(randommzML, mode="onDisk")

for (i in 1:nrow(par)) {
  RRP=ramclustR(ms=final_PT_name,
                st = par$st[i],
                sr=par$sr[i],
                maxt = par$maxt[i],
                timepos = 2,
                sampNameCol = 1,
                featdelim = '_',
                ExpDes = expDes) #Important to note that a csv file is needed as input!
  nClust[i]=length(RRP$cmpd)
  nSing[i]=RRP$nsing
  sizeMax[i]=max(RRP$nfeat)
  sizeMed[i]=median(RRP$nfeat)
  sizeMean[i]=mean(RRP$nfeat)
  nFeat[[i]]=RRP$nfeat
  pdf(file=paste0('clusts_par',i,'.pdf'),width=15,height=8)
  par(mfrow=c(4,5),mar=c(4,4,2,0)+.5)
  clusts=round(c(2:6,seq(7,max(RRP$featclus),length.out = 15)))
  for (c in clusts) {
    plotClust(ram = RRP, clustnr = c, files=files)
  }
  dev.off()
}
cbind(par,nClust,nSing,sizeMax,sizeMean,sizeMed)




###########################################################################################
#### Evaluated output pdfs to determine which RamClustR settings to move fowrward with ####

paramRows <- cbind(par,nClust,nSing,sizeMax,sizeMean,sizeMed)
write.csv(paramRows, file = paste0(resultsPath4, "paramRows.csv"))

finalRamClust=ramclustR(ms=final_PT_name,
                        st = 0.35,
                        sr=0.5,
                        maxt = 5,
                        timepos = 2,
                        sampNameCol = 1,
                        featdelim = '_',
                        ExpDes = expDes) #Important to note that a csv file is needed as input!

finalRamClust$nsing #check how many features does not belong to any clusters
max(finalRamClust$featclus) #check how many clusters were created
Clusts <- finalRamClust$SpecAbund
colnames(Clusts) <- paste0("C", colnames(Clusts))
RamDF=cbind(Clusts,finalRamClust$MSdata[,finalRamClust$featclus==0])

#### Building chromPol ####
if(pol=="POS"){
  if(chrom=="RP"){
    chromPol <- "RP"
  } else {
    chromPol <- "HP"
  }
} else {
  if(chrom=="RP"){
    chromPol <- "RN"
  } else {
    chromPol <- "HN"
  }
}

##### 241028 Hackathon:
#This should be if/else chrompol
#HILIC/RP P/N
#RNRC = RP NEG RAMCLUST
colnames(RamDF) <- paste0(chromPol,'_',colnames(RamDF))
rownames(RamDF) <- rownames(read.csv(final_PT_name))
finalRamClust_temp <- finalRamClust
dateToday <- substr(format(as.Date(Sys.Date()), "%Y%m%d"), 3, 1000)
saveRDS(list(finalRamClust_temp, RamDF), file=paste0(dateToday, '_', chrom, '_', pol, '_Ram.rds'))

##### 241028 Hackathon:
#ChromPol paste checks
write.xlsx(as.data.frame(RamDF), paste0(dateToday, "_", chrom, "_", pol,"_FinalPT.xlsx"), rowNames=T)
write.csv2(as.data.frame(RamDF), paste0(dateToday, "_", chrom, "_", pol,"_FinalPT.csv"))
# save.image(paste0(resultsPath4, dateToday, "_AllData_RP_POS.RData"))

#### Building file which indexes which features ended up in which clusters
FeatOverview <- as.data.frame(cbind(finalRamClust_temp$featclus,
                                    round(finalRamClust_temp$fmz,5),
                                    round(finalRamClust_temp$frt,2),
                                    paste0(finalRamClust_temp$fmz, "_", finalRamClust_temp$frt),
                                    colMedians(finalRamClust_temp$MSdata),
                                    rep("",length(finalRamClust_temp$featclus))))
FeatOverview <- FeatOverview[order(FeatOverview[,1]),]
colnames(FeatOverview) <- c(paste0(chromPol,"C"),"mz","rt","mz_rt", "median_area","potential_type")

FeatOverview[,2] <- as.double(FeatOverview[,2])
FeatOverview[,3] <- as.double(FeatOverview[,3])
FeatOverview[,5] <- as.double(FeatOverview[,5])

if(pol=="POS"){
  adductList <- c(#"[M+H]+",
    "[M+Na]+",
    "[M+K]+",
    "[M+NH4]+",
    "[M+CH3OH]+",
    "[2M+H]+",
    "[M-H2O+H]+",
    "[M+2H]2+")
} else {
  adductList <- c(#"[M-H]-",
    "[M+FA-H]-",
    "[M+Na-2H]-",
    "[M+K-2H]-",
    "[M+Cl]-",
    "[2M-H]-",
    "[M-H2O-H]-",
    "[M-2H]2-")
}

#### Trying to figure out potential types for features: isotope, adduct, etc ####
uniqueClusts <- unique(FeatOverview[,1])
uniqueClusts <- as.integer(uniqueClusts[-which(uniqueClusts == "0")])
uniqueClusts <- uniqueClusts[order(uniqueClusts)]

#Formatting nicely
FeatOverview <- FeatOverview[order(as.integer(FeatOverview[,1])),]

for(i in c(0,uniqueClusts)){
  FeatOverview[which(as.integer(FeatOverview[,1]) == i),] <- FeatOverview[which(as.integer(FeatOverview[,1]) == i)[order(FeatOverview[which(as.integer(FeatOverview[,1]) == i), 2])], ]
}


for(i in uniqueClusts){
  tempClustDF <- FeatOverview[which(FeatOverview[,1] == as.character(i)),]
  max_area <- max(tempClustDF[,5])
  tempClustDF[which(tempClustDF[,5] == max_area),6] <- "Molecule"
  molMz <- tempClustDF[which(tempClustDF[,5] == max_area),2]
  
  #### Calculating potential adducts, -H/+H included in if criteria later  ####
  tempAdductMasses <- data.frame(matrix(ncol=2, nrow=length(adductList)))
  tempAdductMasses[,1] <- adductList
  if(pol == "POS"){
    tempAdductMasses[,2] <- c(#molMz+1.007276,
      molMz+22.989218,
      molMz+38.963158,
      molMz+18.033823,
      molMz+33.033489,
      (molMz*2)+1.007276,
      molMz+(17.00274-0.0005485833),
      (molMz/2)+1.007276)
  } else {
    tempAdductMasses[,2] <- c(#molMz-1.007276,
      molMz+44.998201,
      molMz+20.974666,
      molMz+36.948606,
      molMz+34.969402,
      (molMz*2)-1.007276,
      molMz-19.01839,
      (molMz/2)-1.007276)
  }
  
  #### Finding isotopes of molecule, in-source fragments and adducts
  for(j in 1:nrow(tempClustDF)){
    #### Checking if potential adduct ####
    whichAdduct <- which((tempAdductMasses[,2]-1.01) < tempClustDF[j,2]  & (tempAdductMasses[,2]+1.01) > tempClustDF[j,2])
    
    if(length(whichAdduct) > 0){
      whichAdduct <- whichAdduct[which.min(abs(tempAdductMasses[whichAdduct,2] - tempClustDF[j,2]))]
      tempClustDF[which(tempClustDF[,2] == molMz),6] <- ifelse(pol=="POS", "[M+H]+", "[M-H]-") 
      tempClustDF[j,6] <- tempAdductMasses[whichAdduct,1]
    }
    
    #Checking if isotope
    if(any(tempClustDF[j,2] > (tempClustDF[,2] + 1.001) & tempClustDF[j,2] < (tempClustDF[,2] + 1.008))){
      tempClustDF[j,6] <- "Isotope"
    }
  }
  
  FeatOverview[which(FeatOverview[,1] == as.character(i)),] <- tempClustDF
}

#Saving cluster breakdown data
write.xlsx(as.data.frame(FeatOverview), paste0(dateToday, "_", chrom, "_", pol, "_ClusterBreakdown.xlsx"), rowNames = F)
write.csv2(FeatOverview, paste0(dateToday,  "_", chrom, "_", pol, "_ClusterBreakdown.csv"), row.names = F)

#Saving last image
save.image(paste0(format(as.Date(Sys.Date()), "%Y%m%d"), "_", chrom, "_", pol,"_AllData.RData"))
