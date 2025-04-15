##################################################################
# Pure xcms pipeline from .mzML files to data fit for statistics #
##################################################################

##### 241028 Hackathon - Anton Ribbenstedt + Carl Brunius - Planning of and implementation of v1.6 ####
#- Compartmentalize all "pipes" into separate functions in one wrapper
#- Will have to slaughter IPORT to get the optimization and apply to LaMa finding optimizer
#- Tracking features of aligned and PPd peaks to features

#- Major changes:
#- rtAling within PP of each separate batch with LaMa
#- perform rtAlign on all combined batches
#- Will it drop initial rtAlign if doing it again? Ask Johannes
#- Add specific LaMas to PP object in step 2 already with an early "hard fill" function

##### 241105 Fixes:
# [v] Moving all parameter settings to beginning of script which will be saved in paste0(outPath,"paths_params_functions.RData")
# [v] New default settings for instruments etc.
# [v] rtAlign plots for each separate batch + only sQC all batches


#### Cleaning up prior to start ####
rm(list = ls())
gc()

options(digits = 10)

#### Libraries
library(xcms) # main algorithm for peak picking
library(CMSITools)
library(tidyverse) # for function %>%
library(StatTools) # For the imputation function, 'mvImp', script written by Carl Brunius, https://gitlab.com/CarlBrunius/StatTools
library(BiocParallel) # parallell computing
library(stringi)
library(openxlsx)
library(RColorBrewer)
library(scales)
library(dplyr)
library(batchCorr)
library(ranger)
library(doParallel)
library(MSnbase)
library(IPO)
library(IPO2)
library(RAMClustR)
library(fs)
library(data.table)
library(ggfortify)
library(miscTools)

#### Inputs for the entire "pipeline"
#### Should be an excel config file?7
#Paths and dataset meta data
setwd("/home/antonri/Documents/R/")
inPath <- "/home/antonri/Documents/Pipeline/Input/"
outPath <- "/home/antonri/Documents/Pipeline/Output/"
numThreads <- 14
chrom <- "" #RP / HILIC
pol <- "" #POS / NEG
instrument <- "" #MRT
injToRemove <- c("cond", "SQCcond", "WASH", "iterative", "MeOH", "SST") #Parts of filenames that are to be removed in LCMS_data
set.seed(100)

#Static cwp parameters
minFrac <- 0.2

#IPO RT parameters
minfrac_IPO <- 0.8 #A set value, minimum amount of samples peaks have to be found in to be considered a feature, default good 241028
gap_init_IPO <- c(0.2, 0.7) #If good chrom -> 1:1 alignment, if bad chrom -> two step process in deviating from 1:1 with different penalties, varies with response
gap_extend_IPO <- c(2,3)
response_IPO <- c(1,10) #% of all potential anchor-points; good chrom -> fewer %, bad chrom -> higher %
prof_step_IPO <- c(0.8,1.1) #Finding anchor points mz in bucket, bucket-width in mz
mzvid_IPO <- c(0.015, 0.035)

#### Automatically setting up cwp-starting parameters based on instrument and predetermined knowledge (faster processing if starting closer to the mark)
#### Predetermined knowledge:
# - Output from running IPO2 on datasets
# - Manually inspecting noise in void-volume of random files

optimVars <- c("min_peakwidth", "max_peakwidth", "mzdiff", "ppm" )
if(instrument == "MRT"){
  if(chrom == "RP"){
    if(pol == "POS"){
      # Starting values
      cwParam <- CentWaveParam(
        ppm = 6,
        peakwidth = c(1.7, 22),
        snthresh = 10,
        prefilter = c(5, 45000),
        mzdiff = 0.0017,
        noise = 20000
      )
      
      # Ranges for IPO2 optimization
      # Order of parameters: "min_peakwidth", "max_peakwidth", "mzdiff", "ppm"
      upper <- c(5, 45, Inf, 25) 
      lower <- c(0.5, 10, -Inf, 1)
      
      # Absolute threshold for instrument
      intThreshold <- 20000
      
    } else if (pol == "NEG"){
      #Predetermined knowledge
      #noise: 7500, prefilter: c(3, 20000), snthresh: 10
      #min_peakwidth: 1.650148446, max_peakwidth: 21.210341074, mzdiff: 0.001600235, ppm: 5.745059710 
      
      # Starting values
      cwParam <- CentWaveParam(
        ppm = 6,
        peakwidth = c(1.7, 22),
        snthresh = 10,
        prefilter = c(5, 45000),
        mzdiff = 0.0017,
        noise = 20000
      )
      
      # Ranges for IPO2 optimization
      # Order of parameters: "min_peakwidth", "max_peakwidth", "mzdiff", "ppm"
      upper <- c(5, 45, Inf, 25)
      lower <- c(0.5, 10, -Inf, 1)
      
      # Absolute threshold for instrument
      intThreshold <- 20000
    }
  } else if (chrom == "HILIC"){
    # Starting values
    cwParam <- CentWaveParam(
      ppm = 10,
      peakwidth = c(2.5, 30),
      snthresh = 10,
      prefilter = c(5, 15000),
      mzdiff = -0.001,
      noise = 10000
    )
    
    # Ranges for IPO2 optimization
    # Order of parameters: "min_peakwidth", "max_peakwidth", "mzdiff", "ppm"
    upper <- c(5, 45, Inf, 18)
    lower <- c(0, 5, -Inf, 5)
    
    # Absolute threshold for instrument
    intThreshold <- 10000
    
    
  }
} else if (instrument == "SynaptXS"){
  if(chrom == "RP"){
    if(pol == "POS"){
      
    } else if (pol == "NEG"){
      
    }
  } else if (chrom == "HILIC"){
    if(pol == "POS"){
      
    } else if (pol == "NEG"){
      
    }
  }
}


# Checking OS and setting up BiocParallell based on that
if(Sys.info()[1] == "Windows"){
  register(SnowParam(workers = numThreads, type = "SOCK"))
} else {
  register(MulticoreParam(workers = numThreads))
}

# Setting up paths for all parts of the pipeline and saving as .rds which is loaded for every part of dataset
inPath <- inPath
if (str_sub(inPath, -1, -1) == "/") { # so user can type path with or without / at the end. Needed without below when making the names_list.
  inPath <- str_sub(inPath, 1, -2)
}

if(!dir.exists(outPath)){
  dir.create(outPath)
}

outPath <- paste0(outPath, "/results/")
if(!dir.exists(outPath)){
  dir.create(outPath)
}

resultsPath1 <- paste0(outPath, "Pipe1_IPOPPOpt/") #241028 Hackathon: Functional folder names to describe content - look over names
if(!dir.exists(resultsPath1)){
  dir.create(resultsPath1)
}
resultsPath2 <- paste0(outPath, "Pipe2_PPBatches_IPOrt/")
if(!dir.exists(resultsPath2)){
  dir.create(resultsPath2)
}
resultsPath3 <- paste0(outPath, "Pipe3_rtAlign_BatchCorr/")
if(!dir.exists(resultsPath3)){
  dir.create(resultsPath3)
}
outPath_corrDrift <- paste0(outPath, "Pipe3_rtAlign_BatchCorr/corrDrift/")
if(!dir.exists(outPath_corrDrift)){
  dir.create(outPath_corrDrift)
}
resultsPath4 <- paste0(outPath, "Pipe4_batchAlign_ramClust/")
if(!dir.exists(resultsPath4)){
  dir.create(resultsPath4)
}

#Loading all functions
source.list<-list.files("./R/", full.names = T) #241028 Hackathon: Put in CMSITools - librarification
sapply(source.list, source)

#################################
#### Part 0 :: Double checking files prior to starting pipeline
#################################

#Removing "_CLMC" in file name if any
allFiles <- list.files(inPath,
                       recursive = T,
                       full.names = T)

if(any(grepl("_CLMC", allFiles))){
  for(i in 1:length(allFiles)){
    file.rename(allFiles[i],
                gsub("_CLMC", "", allFiles[i]))
  }
}

#Removing any spaces
allFiles <- list.files(inPath,
                       recursive = T,
                       full.names = T)

if(any(grepl(" ", allFiles))){
  for(i in 1:length(allFiles)){
    file.rename(allFiles[i],
                gsub(" ", "-", allFiles[i]))
  }
}

#If batchnames in filenames are not in 00 format, change it
allFiles <- list.files(inPath,
                       recursive = T,
                       full.names = T)

for(i in 1:length(allFiles)){
  sampBatch <- strsplit(basename(allFiles[i]), "_")[[1]][2]
  batchNumb <- strsplit(strsplit(sampBatch, "B")[[1]][2], "W")[[1]][1]
  if(length(batchNumb) < 2){
    newName <- gsub("B[0-9]W",
                    paste0("B0", as.integer(batchNumb), "W"),
                    allFiles[i])
    file.rename(allFiles[i],
                newName)
  }
}

#Making a data.table of summary statistics of files in each batch folder
allFiles <- list.files(inPath,
                       recursive = T)
batchFolders <- unique(path_dir(allFiles))

fileSummary <- data.table(matrix(ncol=8))
setnames(fileSummary,
         paste0("V",c(1:8)),
         c("Batch", "sQCs", "ltQCs", "Blanks", "CondPlasma",  "Samples", "nTotal", "outliersFileSize"))
fileSummary <- fileSummary[, lapply(.SD, as.integer), by=Batch]
fileSummary <- fileSummary[, Batch:=as.character(Batch)]

for(i in 1:length(batchFolders)){
  batchFiles <- allFiles[grepl(batchFolders[i], allFiles)]
  nSQCs <- sum(grepl("sQC",
                     batchFiles,
                     ignore.case = T))
  nltQCs <- sum(grepl("ltQC",
                      batchFiles,
                      ignore.case = T))
  nBlanks <- sum(grepl("blank",
                       batchFiles,
                       ignore.case = T))
  nCond <- sum(grepl("cond",
                     batchFiles,
                     ignore.case = T))
  nSamples <- length(batchFiles) - nSQCs - nltQCs - nBlanks - nCond
  
  #Calculating potential outliers in too low file size
  x <- file.size(paste0(inPath, "/", batchFiles))
  qnt <- quantile(x,
                  probs=c(.1),
                  na.rm=T)
  H <- 1.5 * IQR(x,
                 na.rm=T)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  nOutliers <- sum(is.na(y))
  
  
  fileSummary <- rbind(fileSummary,
                       cbind(Batch=batchFolders[i],
                             sQCs=nSQCs,
                             ltQCs=nltQCs,
                             Blanks=nBlanks,
                             CondPlasma=nCond,
                             Samples=nSamples,
                             nTotal=length(batchFiles),
                             outliersFileSize=nOutliers))
}
fileSummary <- fileSummary[-1,]

fwrite(fileSummary,
       file=paste0(format(Sys.Date(),format="%y%m%d"), "_FileSummary.csv"))  ####241028 Hackathon: Put all in support functions


##### Reading all files and making a meta data df to use throughout the pipeline
LCMS_data <- CMSITools::getFiles(inPath, pattern=".mzML") %>% CMSITools::cleanFiles(injToRemove)    ##### 241028 Hackathon: Update for MRTs

# LCMS_data is needed for meta input for getBatch()                                                                         ##### 241028 Hackathon: Should be available since before redundant
if (chrom == 'RP') {
  if (pol == 'POS') {
    LCMS_data <- LCMS_data %>% CMSITools::getRP() %>% CMSITools::getPOS()
  }
  if (pol == 'NEG') {
    LCMS_data <- LCMS_data %>% CMSITools::getRP() %>% CMSITools::getNEG()
  }
} else if (chrom == 'HILIC') {
  if (pol == 'POS') {
    LCMS_data <- LCMS_data %>% CMSITools::getHILIC() %>% CMSITools::getPOS()
  }
  if (pol == 'NEG') {
    LCMS_data <- LCMS_data %>% CMSITools::getHILIC() %>% CMSITools::getNEG()
  }
}
LCMS_data <- add_group_column(LCMS_data)


setwd(outPath)
save.image(paste0(outPath,"paths_params_functions.RData")) #241028 Hackathon: Better name?

stop("Manual interpretation required")



#################################
#### Part 1 :: IPO2 optimization of peak picking parameters :: Per batch ####
#################################

# Extracting only sQCs for PP optimization
if (chrom == "RP") { 
  sQCs <- LCMS_data %>% getRP()
} else if (chrom == "HILIC") { 
  sQCs <- LCMS_data %>% getHILIC() 
}
if (pol == "POS") { 
  sQCs <- sQCs %>% getPOS() %>% getQC("sQC")
} else if (pol == "NEG") { 
  sQCs <- sQCs %>% getNEG() %>% getQC("sQC") 
}  

# Optimizing peak picking in all batches
batches <- unique(sQCs$batch)
batchNum <- unique(stri_extract_first_regex(sQCs$batch, "[0-9]+"))

for(i in 1:length(batches)){
  print("\n\n\n")
  print(batchNum[i])
  print("\n\n\n")
  sink(paste0(resultsPath1, "IPO_log_", batchNum[i],".txt"))
  startTime <- Sys.time()
  sQCs_batch_specific <- sQCs[sQCs$batch==batches[i],]
  sQCs_batch_names <- sQCs_batch_specific %>% extractNames() %>% sample(min(5, length(.)))
  
  # Running optimization on 5 random sQCs from batch
  optParams <- optimXCMS(sQCs_batch_names,
                         cwParam,
                         optimVars,
                         upper=upper,
                         lower=lower,
                         "NLOPT_LN_NELDERMEAD",
                         verbose = T)
  
  names(optParams$solution) <- optimVars
  endTime <- Sys.time() - startTime
  print(paste0("Time to compute: ", endTime))
  sink() #Ending capture of all console output during optimization
  
  saveRDS(optParams, paste0(resultsPath1, "batch", batchNum[i], "_opt_params.rds"))
  
  #Saving optimized settings from optParams into a txt document
  write.csv(as.data.frame(optParams$solution), paste0(resultsPath1, "/batch", batchNum[i], "_opt_params.csv"), row.names=T)
  
  #Peak picking the sQCs with optimized parameters
  pd <- data.frame(sample_name = sQCs_batch_names,
                   sample_group = c(rep("sQC", length(sQCs_batch_names))),
                   stringsAsFactors = FALSE)
  
  cwp <- CentWaveParam(peakwidth = c(optParams$solution[1], optParams$solution[2]), noise = cwParam@noise,
                       prefilter = cwParam@prefilter,
                       ppm = optParams$solution[4],
                       mzdiff = optParams$solution[3])
  
  
  #### Consider making as xcmsSet in order to use old IPO ####
  sQCs_MSExp <- readMSData(files = sQCs_batch_names,
                           # pdata = pd,
                           mode="onDisk")
  # sQCs_MSExp <- filterSpectra(sQCs_MSExp, filterEmptySpectra)
  sQCs_MSExp <- findChromPeaks(sQCs_MSExp, param=cwp)
  
  saveRDS(sQCs_MSExp, file=paste0(resultsPath1, "batch", batchNum[i], "_sQCs.rds"))
}

# Saving the non-optimized parameters for future use                                            #####241028 Hackathon: Remove? Already present in .RData file
outFile <- "/static_params.txt"
if (!file.exists(paste0(resultsPath1, outFile))) {
  sink(paste0(resultsPath1, outFile))
  cat('noise,', cwParam@noise, "\n")
  cat('prefilter,', cwParam@prefilter[1], ',', cwParam@prefilter[2], "\n")
  cat('snthresh,', cwParam@snthresh, "\n")
  cat('polarity,', pol, "\n")
  cat('chromatography,', chrom, "\n")
  sink()
}

# Making a .csv file with collected optimized values from all batches                            #####241028 Hackathon: Support function
optParamsDF <- data.frame(matrix(ncol=5))
colnames(optParamsDF) <- c("batch",
                           "min_peakwidth",
                           "max_peakwidth",
                           "ppm",
                           "mzdiff")
for(batch in batches){
  i <- which(batches %in% batch)
  
  batchNum <- as.character(stri_extract_first_regex(batch, "[0-9]+"))
  # if(length(batchNum) < 2)){
  #   batchNum <- 
  # }
  
  batch_param_file <- read.csv(paste0(resultsPath1, "/batch", batchNum, "_opt_params.csv"))
  
  for (j in 1:nrow(batch_param_file)){
    param <- batch_param_file[j,1]
    value1 <- batch_param_file[j,2]
    if (param == "min_peakwidth") {
      min_peakwidth <- as.numeric(value1)
    } else if (param == 'max_peakwidth') {
      max_peakwidth <- as.numeric(value1)
    } else if (param == 'mzdiff') {
      mzdiff <- as.numeric(value1)
    } else if (param == 'ppm') {
      ppm <- as.numeric(value1)
    }
  }
  
  optParamsDF[i,] <- c(batch,
                       min_peakwidth,
                       max_peakwidth,
                       ppm,
                       mzdiff)
}

write.csv(optParamsDF, paste0(resultsPath1,"/opt_params.csv"))
save.image(paste0(resultsPath1, format(as.Date(Sys.Date()), "%Y%m%d"), "_Part1.RData"))


#################################
#### Part 2 :: Peak picking with optimized parameters :: Per batch ####
#################################

#Cleaning up prior to starting
rm(list = ls())
gc()

load("paths_params_functions.RData")
setwd(outPath)

# Static Params file from IPO2 optimization read
# static_param_file <- paste0(resultsPath1, "/", "static_params.txt")                    ##### 241028 Hackathon: Already in paths_... .RData
# lines <- readLines(static_param_file) 
# 
# for (line in lines) {
#   line <- trimws(line)
#   param <- word(line, 1, sep = "\\s*,\\s*")
#   value1 <- word(line, 2, sep = "\\s*,\\s*")
#   value2 <- word(line, 3, sep = "\\s*,\\s*")
#   if (param == "noise") {
#     noise <- as.numeric(value1)
#   } else if (param == 'prefilter') {
#     prefilt <- as.numeric(c(value1, value2))
#   } else if (param == 'snthresh') {
#     snthr <- as.numeric(value1)
#   } else if (param == 'polarity') {
#     pol <- value1
#   } else if (param == 'chromatography') {
#     chrom <- value1
#   }
# }

# Reading file with all optimized parameters from IPO2 here
opt_param_file <- paste0(resultsPath1, "/", "opt_params.csv")
opt_param_df <- read.csv(opt_param_file, sep=',')

# Peak picking for each unique batch specified in the opt_params.csv outFile from IPO2 optimization
batches <- opt_param_df$batch

##### 241028 Hackathon: Move LCMS_data out here and specify new name within loop

for (batch in batches) {
  
  
  # Collecting optimized parameters
  peakwidth <- as.numeric(c(opt_param_df$min_peakwidth[opt_param_df$batch == batch], opt_param_df$max_peakwidth[opt_param_df$batch == batch]))
  ppm <- as.numeric(opt_param_df$ppm[opt_param_df$batch == batch])
  mzdiff <- opt_param_df$mzdiff[opt_param_df$batch == batch]
  
  # Enforcing double digit single digit system for sorting-of-data purposes
  LCMS_data_batch <- LCMS_data[grep(batch, LCMS_data$batch),] # returns only the rows of data matching the batch no
  LCMS_data_names <- extractNames(LCMS_data_batch)
  
  # Preparing data structure for reading all files related to 
  pd <- data.frame(
    sample_name = sub(basename(LCMS_data_names), pattern = ".mzML",
                      replacement = "", fixed = TRUE),
    sample_group = LCMS_data_batch$group, stringsAsFactors = FALSE)
  
  raw_data <- readMSData(files = LCMS_data_names,
                         pdata = new("NAnnotatedDataFrame", pd),
                         mode = "onDisk")
  
  raw_data <- filterEmptySpectra(raw_data)
  
  # runs the peak picking
  cwp <- CentWaveParam(peakwidth=peakwidth,
                       ppm=ppm,
                       mzdiff=mzdiff,
                       prefilter=cwParam@prefilter,
                       noise=cwParam@noise,
                       snthresh=cwParam@snthresh)                    ##### 241028 Hackathon: Make sure right variable names after removing static_params.txt [v]
  
  xdata <- findChromPeaks(raw_data,
                          param=cwp)
  
  batchNum <- as.character(stri_extract_first_regex(opt_param_df$batch, "[0-9]+"))[which(batches %in% batch)]
  saveRDS(xdata, file= paste0(resultsPath2, "xcms_", pol, "_PreBW_", batchNum, ".rds"))
}

save.image(file=paste0(resultsPath2, format(as.Date(Sys.Date()), "%Y%m%d"), "_Part2.RData"))






#################################
#### Part 3 :: IPO retAlign parameter optimization :: Unified dataset ####
#################################

##### 241028 Hackathon: If basing on LaMas, this will be completely changed
##### 241028 Hackathon: If we want to optimize parameters for LaMa alignment, it should be here!

# Cleaning up prior to start
rm(list = ls())
gc()

load("paths_params_functions.RData")
setwd(outPath)

allRetCor <- list.files(resultsPath1, full.names=T, pattern = ".rds")
allRetCor <- allRetCor[grepl("sQC", allRetCor)]


##### 241028 Hackathon: From here to 503 modify to read every batch
xset <- readRDS(allRetCor[1])

##### 241028 Hackathon: This will be redundant when running all batches
# if(!is.null(batchDist) && batchDist < 50 && batchDist != 999){                                                                    
#   whichBatchesToRun <- seq(batchFrom, length(allRetCor), by=batchDist)
#   #Adding last batch if not already included
#   if(whichBatchesToRun[length(whichBatchesToRun)] != length(allRetCor)){
#     whichBatchesToRun <- c(whichBatchesToRun,
#                            length(allRetCor))
#   }
# } else if(length(allRetCor) > 1) {
#   # print(allRetCor)
#   whichBatchesToRun <- seq(2, length(allRetCor), by=1)
# }
# 
# if((length(whichBatchesToRun)-1) == 0){
#   
# } else if((whichBatchesToRun[length(whichBatchesToRun)] - whichBatchesToRun[length(whichBatchesToRun)-1]) < floor(batchDist/2)){
#   whichBatchesToRun <- whichBatchesToRun[-(length(whichBatchesToRun)-1)]
# }

# print(allRetCor)
# print("WhichBatchesToRun")

##### 241028 Hackathon: From 2 and upward since we're running all batches
#Combining all batches                                                            
if(length(allRetCor) > 1){
  for(i in 2:length(allRetCor)){ #whichBatchesToRun
    xset <- c(xset, readRDS(allRetCor[i]))
  }
}

# Setting parameters for rtAlign                                                  
##### 241028 Hackathon: Will be redundant when using LaMas
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$minfrac <- minfrac_IPO #a set value, minimum amount of samples peaks have to be found in to be considered a feature
retcorGroupParameters$gapInit <- gap_init_IPO
retcorGroupParameters$gapExtend <- gap_extend_IPO
retcorGroupParameters$response <- response_IPO
retcorGroupParameters$profStep <- prof_step_IPO
retcorGroupParameters$mzwid <- mzvid_IPO

xset <- as(xset, 'xcmsSet')

numThreads = ceiling(numThreads/4)

sink(file=paste0(resultsPath2, "IPO_RT_log.txt"))
time.RetGroup <- system.time({ # measuring time
  resultRetcorGroup <-
    optimizeRetGroup(xset = xset, #object after peak picking
                     params = retcorGroupParameters, #parameters
                     # BPPARAM = bpparam(),
                     nSlaves = numThreads,
                     subdir = NULL,
                     plot = TRUE)
})
print(paste0("Time taken: ", time.RetGroup))
sink()

saveRDS(resultRetcorGroup, paste0(resultsPath2, "/retcorBestSettings_global.rds"))

outFile <- paste0(resultsPath2, "retcorBestSettings_global.txt")
if (!file.exists(outFile)) {
  sink(outFile)
  cat('gapInit,', resultRetcorGroup$best_settings$gapInit, "\n")
  cat('gapExtend,', resultRetcorGroup$best_settings$gapExtend, "\n")
  cat('profStep,', resultRetcorGroup$best_settings$profStep, "\n")
  cat('response,', resultRetcorGroup$best_settings$response, "\n")
  cat('bw,', resultRetcorGroup$best_settings$bw, "\n")
  cat('mzwid,', resultRetcorGroup$best_settings$mzwid, "\n")
  cat('minfrac,', resultRetcorGroup$best_settings$minfrac, "\n")
  sink()
}

save.image(paste0(resultsPath2, format(as.Date(Sys.Date()), "%Y%m%d"), "_Part3.RData"))






#################################
#### Part 4 :: Applying retAlign optimized parameters :: Unified dataset ####
#################################

#### Cleaning up prior to start ####
rm(list = ls())
gc()

load("paths_params_functions.RData")
setwd(outPath)

##### 241028 Hackathon: This only needed when IPO, will need other input when doing LaMas
# Read rtAlign optimized parameters from .txt-files                                   
retcor_param_file <- paste0(resultsPath2, "retcorBestSettings_global.txt")

lines <- readLines(retcor_param_file) 
for (line in lines) {
  line <- trimws(line)
  param <- word(line, 1, sep = "\\s*,\\s*")
  value <- word(line, 2, sep = "\\s*,\\s*")
  if (param == "response") {
    response <- as.numeric(value)
  } else if (param == 'gapInit') {
    gapInit <- as.numeric(value)
  } else if (param == 'gapExtend') {
    gapExtend <- as.numeric(value)
  } else if (param == 'minfrac') {
    minFraction <- as.numeric(value)
  } else if (param == 'bw') {
    bw <- as.numeric(value)
  } else if (param == 'profStep') {
    binSizeObiwarp <- as.numeric(value)
  } else if (param == 'mzwid') {
    binSize <- as.numeric(value)
  }
}

# Create list of input files here:
inFiles <- list.files(resultsPath2, "*PreBW*") # input rds files from script 3 output

# create merged rds-object:
xdata <- readRDS(paste0(resultsPath2, inFiles[1]))
if(length(inFiles) > 1){
  for (file in inFiles[2:length(inFiles)]) {
    rds_obj <- readRDS(paste0(resultsPath2, file))
    xdata <- c(xdata, rds_obj)
  }
}

saveRDS(xdata, paste0(resultsPath3, "xdata-mergedPreAdjustRtime.rds"))

#Alignment ---------------------------
# register(SerialParam()) #Used to be needed due to memory errors                              ##### 241028 Hackathon: LaMa warping hard set parameters in beginning of script
pgp <- ObiwarpParam(binSize=binSizeObiwarp,
                    response=response,
                    gapInit=gapInit,
                    gapExtend=gapExtend,
                    rtimeDifferenceThreshold = 50)
xdata <- adjustRtime(xdata, param = pgp)
saveRDS(xdata, paste0(resultsPath3, "xdata-PostAdjustRtime.rds"))



pdp_pregroup <- PeakDensityParam(sampleGroups = xdata@phenoData@data$sample_group,
                                 minFraction = minFrac,
                                 bw = bw, 
                                 binSize = binSize)
xdata <- groupChromPeaks(xdata, param = pdp_pregroup)
# print(dim(getTable(xdata)))

#Printing summary statistics regarding the grouped peak table
sink(paste0(resultsPath3, "postGroupChromPeaks.txt"))
png(paste0(resultsPath3, "getTable.png"))
getTable(xdata)
sink()
dev.off()

#Save the corrected grouped object
saveRDS(xdata, paste0(resultsPath3, "XCMS_Grouped.rds"))

#Define Colors for rtAlign plot-------------------                                              ##### 241028 Hackathon: Support function for plotting the entire experiment
batches <- lapply(xdata@phenoData@data$sample_name, get_batch_nos)
batches <- unlist(batches) %>% unique
max_batches <- length(batches)

group_colors <- paste0(hue_pal()(max_batches))
col_vector <- c(1:max_batches) # Added by FJ 
names(group_colors) <- paste0("B", col_vector) #Assign a different color to each batch

sampleColVec <- 1:length(xdata@phenoData@data$sample_name)
for(i in 1:length(batches)){
  sampleColVec[grepl(paste0("B", batches[i], "W"), xdata@phenoData@data$sample_name)] <- group_colors[i]
}

whichFileNames_sQC <- xdata@phenoData@data$sample_name[which(grepl("sQC", xdata@phenoData@data$sample_group))]
sampleColVec_sQC <- sampleColVec[which(grepl("sQC", xdata@phenoData@data$sample_group))]

#Create Alignment Plot For sQCs in all batches--------------
jpeg(paste0(resultsPath3, "alignmentPlot_allBatches.jpg"))
plotAdjustedRtime(MSnbase::filterFile(xdata,
                                      file=which(xdata@phenoData@data[,1] %in% whichFileNames_sQC)),
                  col=sampleColVec_sQC,
                  lwd = 2 #line width of curve
)
dev.off()

#Create Alignment Plot For Each Separate Batch--------------
for(i in 1:max_batches){
  if(i < 10){
    currBatch <- as.character(paste0("0",i))
  } else {
    currBatch <- as.character(i)
  }
  
  whichFileNames_batch <- which(grepl(paste0("B",currBatch,"W"), xdata@phenoData@data$sample_name)) #xdata@phenoData@data$sample_name[which(grepl(currBatch, xdata@phenoData@data$sample_name))]
  sampleColVec_batch <- sampleColVec[which(grepl(paste0("B",currBatch,"W"), xdata@phenoData@data$sample_name))]
  
  jpeg(paste0(resultsPath3, "alignmentPlot_B", currBatch,".jpg"))
  plotAdjustedRtime(filterFile(xdata,
                               whichFileNames_batch),
                    col=sampleColVec_batch,
                    lwd = 2 #line width of curve
  )
  dev.off()
}





#Extract the feature table-------------------
FT_NOFill <- getTable(xdata)
saveRDS(FT_NOFill, file = paste0(resultsPath3, "FT_NOFill.rds"))

#Make tracking document for all features-------------------                                             ##### 241028 Hackathon:

trackFeats <- data.table(data.frame(Feature=colnames(FT_NOFill),
                                    Grouping=rep(T, ncol(FT_NOFill)),
                                    PostBlankFilter=rep(T, ncol(FT_NOFill)),
                                    PostBatchCorr=rep(T, ncol(FT_NOFill))))

saveRDS(trackFeats, paste0(resultsPath3,"trackFeats.rds"))

#Hard filling-------------------
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())

Missing_beforeFill<- apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
                           FUN = function(z) sum(is.na(z)))
Missing_afterFill<- apply(featureValues(xdata), MARGIN = 2,
                          FUN = function(z) sum(is.na(z)))
Miss <-  cbind.data.frame(Missing_beforeFill, Missing_afterFill)

sink(paste0(resultsPath3, "Missingness_hardfill.txt"))
FT_HardFill <- getTable(xdata)
sink()

saveRDS(FT_HardFill, file = paste0(resultsPath3, "FT_HardFill.rds"))                                 ##### 241028 Hackathon: Check what's actually necessary
write.csv(Miss, file = paste0(resultsPath3, "Missingness.csv"))
saveRDS(xdata, file = paste0(resultsPath3, "xdata_filled.rds"))

#Blank filter preparations here-------------------                                                   ##### 241028 Hackathon: Flag high blank features here but do not remove
#Checking if any blanks in samples

##### 241028 Hackathon:
#- 1. Within each batch: Check if 80th percentile of each feature is > 3x blank value   --- Parameters to be set at beginning, but these are default values
#- 2. Flag shitty features per batch
#- 3. Remove feature if > 50% flags amongst all batches    --- Parameters to be set at beginning, but these are default values

if(sum(grepl("Blank", rownames(FT_HardFill), ignore.case = T)) > 0){
  #Removing blank samples from DF
  blankDF <- as.data.frame(FT_HardFill[grepl("Blank", rownames(FT_HardFill), ignore.case = T), ])
  if(dim(blankDF)[1] > dim(blankDF)[2]){
    blankDF <- t(blankDF)
  }
  FT_HardFill <- FT_HardFill[-which(grepl("Blank", rownames(FT_HardFill), ignore.case = T)),]
  
  #Checking which batches have blanks (not all do since MRTs use 1 blank for two batches)
  if(length(batches) > 1){
    batchVec <- c()
    for(i in 1:length(rownames(blankDF))){
      batchVec <- c(batchVec,
                    as.integer(gsub("B", "", strsplit(strsplit(rownames(blankDF)[i], "_")[[1]][2],
                                                      "W")[[1]][1]
                    )))
    }
    batchVec_unique <- unique(batchVec)
  } else {
    batchVec <- as.integer(batches)
    batchVec_unique <- batches
  }
  
  
  hasBlank <- rep(T, max_batches)
  for(i in 1:max_batches){
    if(!any(grepl(i, batchVec_unique))){
      hasBlank[i] <- F
    }
  }
} else {
  hasBlank <- rep(F, max_batches)
}

#Making a DF which contains true or false for each feature in each batch
#To be used to determine if we keep a feature or not
featPoorInBlankDF <- featLowIntDF <- as.data.frame(matrix(data=F, nrow=max_batches, ncol=ncol(blankDF)))
rownames(featPoorInBlankDF) <- rownames(featLowIntDF) <- c(1:max_batches)
colnames(featPoorInBlankDF) <- colnames(featLowIntDF) <- colnames(blankDF)

#Actually comparing blanks with real samples
for(i in 1:max_batches){
  if(i < 10){
    currBatch <- as.character(paste0("0",i))
  } else {
    currBatch <- as.character(i)
  }
  
  #Subsetting blanks and checking if blanks from current batch or which was closest previous batch
  if(any(hasBlank)){
    if(hasBlank[i]){
      batchBlankDF <- as.data.frame(blankDF[grepl(as.integer(currBatch), batchVec),])
      if(dim(batchBlankDF)[2] < dim(batchBlankDF)[1]){
        batchBlankDF <- t(batchBlankDF)
      }
    } else {
      
      whichBeforeBatch <- i - 1
      while(whichBeforeBatch != 0){
        if(hasBlank[whichBeforeBatch]){
          if(whichBeforeBatch < 10){
            currBatch <- as.character(paste0("0", whichBeforeBatch))
          } else {
            currBatch <- as.character(whichBeforeBatch)
          }
          
          batchBlankDF <- blankDF[grepl(currBatch, batchVec), ]
          break
        } else {
          whichBeforeBatch <- whichBeforeBatch - 1
        }
      }
      if(whichBeforeBatch == 0){
        stop("No previous batch to use!!! Double check what's going on")
      }
    }
    
    #Taking medians of each column turning it into a vector (if only one blank will be the same value)
    batchBlankDF <- colMedians(batchBlankDF)
    
    batchBlankDF[is.na(batchBlankDF)] <- 0
  }
  
  #Subsetting current batch
  batchDF <- FT_HardFill[which(grepl(paste0("B",currBatch,"W"), rownames(FT_HardFill))),]
  batchDF[is.na(batchDF)] <- 0
  
  
  
  #Looping through all features and comparing against current batch
  for(l in 1:ncol(batchDF)){
    
    #Blank flagging
    if(any(hasBlank)){
      if(quantile(batchDF[,l], probs=0.8) < (batchBlankDF[l]*3)){
        featPoorInBlankDF[i,l] <- T
      }
    }
    
    #Intensity filter flagging
    if(length(intThreshold) > 0){
      if(quantile(batchDF[,l], probs=0.8) < intThreshold){
        featLowIntDF[i,l] <- T
      }
    }
  }
}

#Remove features where 80th percentile < 3x blank intensity in > 50% of batches
if(any(hasBlank)){
  toRemoveBlank <- which((colSums(featPoorInBlankDF)/nrow(featPoorInBlankDF)) > 0.5)
  LCMS_data <- LCMS_data[-which(grepl("blank", LCMS_data$sample, ignore.case = T)), ]
} else {
  toRemoveBlank <- NULL
}

#Remove features where 80th percentile < absolute threshold in > 50% of batches
if(length(intThreshold) > 0){
  toRemoveInt <- which((colSums(featPoorInBlankDF)/nrow(featPoorInBlankDF)) > 0.5)
} else {
  toRemoveInt <- NULL
}

#If any of the filters applied
toRemoveTotal <- unique(c(toRemoveBlank, toRemoveInt))

#If any non-null filters
if(length(toRemoveTotal) > 0){
  FT_HardFill <- FT_HardFill[,-toRemoveTotal]
  saveRDS(FT_HardFill, file = paste0(resultsPath3, "FT_HardFill_Filtered.rds"))
  filtered <- T
  
  #Documenting number of features removed in each filter
  sink(paste0(resultsPath3, "blankAndIntFilterRemoval.txt"))
  print(paste0("Blank filter % removed: ", length(toRemoveBlank)/ncol(FT_HardFill)))
  print(paste0("Int filter % removed: ", length(toRemoveInt)/ncol(FT_HardFill)))
  sink()
  
  trackFeats$PostBlankFilter[toRemoveTotal] <- !(as.logical(toRemoveTotal))
  
  saveRDS(trackFeats, paste0(resultsPath3,"trackFeats.rds"))
} else {
  filtered <- F
  sink(paste0(resultsPath3, "blankAndIntFilterRemoval.txt"))
  print("No filters applied")
  sink()
}

saveRDS(list(filtered,
             toRemoveBlank,
             toRemoveInt,
             toRemoveTotal), file=paste0(resultsPath3, "BlankFilter.rds"))


#Imputation-------------------
sink(paste0(resultsPath3, "mvImpWrap_log.txt"))
if (sum(Missing_afterFill) != 0) {
  FT_fill_imputation <- mvImpWrap(MAT = FT_HardFill , method = "RF", guess = "minVar", forceZero = TRUE)
} else {
  FT_fill_imputation <- FT_HardFill
}
sink()

#Recording % of missing after imputation
sink(paste0(resultsPath3, "percentageMissingAfterImputation.txt"))
print(paste0("Less than zero: ", sum(FT_fill_imputation < 0) / length(FT_fill_imputation)))
print(paste0("Equals zero: ", sum(FT_fill_imputation == 0) / length(FT_fill_imputation)))
print(paste0("Less than or equal zero: ", sum(FT_fill_imputation <= 0) / length(FT_fill_imputation)))
sink()

FT_fill_imputation[FT_fill_imputation < 0] <- 0  

write.csv(FT_fill_imputation, file = paste0(resultsPath3, "FT_fill_imputation.csv"))
saveRDS(FT_fill_imputation, file = paste0(resultsPath3, "FT_fill_imputation.rds"))
saveRDS(LCMS_data, file = paste0(resultsPath3, "LCMS_data_no_blanks.rds"))

save.image(paste0(resultsPath3, format(as.Date(Sys.Date()), "%Y%m%d"), "_Part4.RData"))

#################################
#### Part 5 :: Align batches and batchCorr :: Unified dataset and per batch ####
#################################

#### Cleaning up in-between ####
rm(list = ls())
gc()

load("paths_params_functions.RData")
setwd(outPath)
LCMS_data <- readRDS(paste0(resultsPath3, "LCMS_data_no_blanks.rds"))

# Loading all necessary datasets
#Checking if any filters applied and in that case reading filtered hard-fill dataset
filters <- paste0(resultsPath3, "BlankFilter.rds")
if(filters[[1]] == T){
  FT_HardFill <- readRDS(paste0(resultsPath3, "FT_HardFill_Filtered.rds"))
} else {
  FT_HardFill <- readRDS(paste0(resultsPath3, "FT_HardFill.rds"))
}

FT_NOFill <- readRDS(paste0(resultsPath3, "FT_NOFill.rds"))

FT_fill_imputation <- readRDS(paste0(resultsPath3, "FT_fill_imputation.rds"))

filenames <- rownames(FT_fill_imputation)
batch_codes <- get_batch_codes(filenames)
groups <- get_groups(filenames)

dir.create(paste0(resultsPath3,"/AlignBatches"))

#### Logic for one vs several batches and if no align candidates
if(length(unique(batch_codes)) > 1){
  peakIn <- peakInfo(PT = FT_HardFill, sep = '@', start = 1)
  alignBat <- alignBatches(peakInfo = peakIn, PeakTabNoFill = FT_HardFill, PeakTabFilled = FT_fill_imputation, batches = batch_codes, sampleGroups = groups,
                           selectGroup = "sQC", reportPath = paste0(resultsPath3,"/AlignBatches")) # was QC from gabe
  if(!is.null(alignBat)){
    sink(paste0(resultsPath3,"PTAlignment.txt"))
    print("Data was aligned")
    sink()
    PT_alignBat=alignBat$PTalign
  } else {
    sink(paste0(resultsPath3,"PTAlignment.txt"))
    print("Data was NOT aligned, using imputation PT")
    sink()
    PT_alignBat=FT_fill_imputation
  }
} else {
  sink(paste0(resultsPath3,"PTAlignment.txt"))
  print("Only one batch, using imputation PT")
  sink()
  PT_alignBat <- FT_fill_imputation
}

# --- correct drifts ---
batch_codes_unique <- unique(batch_codes)
batch_codes_unique <- batch_codes_unique[order(as.integer(unique(get_batch_nos(filenames))))]
correctedBatches <- list()

##### 241028 Hackathon:
#- Modifications of batchCorr
#- If sample before first sQC or after last sQC:
#- 1. Drop samples
#- 2. Apply closest correction factors
#- 3. Add artificial sQC by copying closest sQC


#Checking what the highest sQC number is through all batches                                        ##### 241028 Hackathon: This should be specified per batch before starting pipe
# sQCNames <- rownames(FT_fill_imputation)[grepl("sQC", rownames(FT_fill_imputation))]              ##### 241106 Anton: Realized that this will be completely redundant once the sQC script kicks in
# highestQCInt <- 1                                                                                 #####               Potential problem when two or more sQCs are missing....
# highestQC
# for(i in 1:length(sQCNames)){
#   tempQC <- strsplit(sQCNames[i],
#                      "_")[[1]][5]
#   if(as.integer(gsub("sQC", "", tempQC)) > highestQCInt){
#     highestQC <- tempQC
#   }
# }

#Actual batch correction
for (batch_code in batch_codes_unique) {
  batch_no <-  get_batch_nos(batch_code)
  print("--------------------------------------------------")
  print(paste0("Drift correction for batch ", batch_no, ":"))
  
  selected_batch <- retrieveBatch(peakTable = PT_alignBat, meta = LCMS_data, select = batch_code)
  batchDir <- paste0(resultsPath3, batch_code)
  if(!dir.exists(batchDir)){
    dir.create(batchDir)
  }
  
  #### 241106 Anton: Redundant as blank filter has removed all the blanks if any are present
  #Remove potential blanks
  # if(any(grepl("Blank", selected_batch$meta$sample, ignore.case = T))){
  #   selected_batch$peakTable <- selected_batch$peakTable[-which(grepl("Blank", rownames(selected_batch$peakTable), ignore.case = T)), ]
  #   selected_batch$meta <- selected_batch$meta[-which(grepl("Blank", selected_batch$meta$sample, ignore.case = T)), ]
  # }
  
  #Saving non-normalized batch PCAs
  ##### 241028 Hackathon:
  #- How to look at post-correction PCAs?
  #- 1. Project into old PCA (will give a better view of how the situation improved)
  #- 2. Make a new PCA (might look like there still are drift effects)
  pca_res <- prcomp(selected_batch$peakTable, scale. = TRUE)
  png(filename=paste0(batchDir, "/", format(as.Date(Sys.Date()), "%Y%m%d"), "_", chrom, "_", pol,"_PreDriftCorr-Dots.png"), width=2559, height=1376, res=120)
  print(autoplot(pca_res, data=selected_batch$meta, colour="group", size=3)  + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))) #autoplot(pca_res, data=meta, colour="Type", size=3, label=T, label.size=3)
  dev.off()

  rownames(pca_res$x) <- rownames(selected_batch$peakTable)
  png(filename=paste0(batchDir, "/", format(as.Date(Sys.Date()), "%Y%m%d"), "_", chrom, "_", pol,"_PreDriftCorr-Labels.png"), width=2559, height=1376, res=120)
  print(autoplot(pca_res, data=selected_batch$meta, colour="group", shape=F, label=T, label.label="filename", label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow")))
  dev.off()

  rownames(pca_res$x) <- selected_batch$meta$sample
  png(filename=paste0(batchDir, "/", format(as.Date(Sys.Date()), "%Y%m%d"), "_", chrom, "_", pol,"_PreDriftCorr-Numbers.png"), width=2559, height=1376, res=120)
  print(autoplot(pca_res, data=selected_batch$meta, colour="group", shape=F, label=T, label.label="injection", label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow")))
  dev.off()

  # Let user decide on G-parameters
  sink(file=paste0(batchDir,"/",batch_code,"_output.txt"))
  drift_correction <- correctDrift(peakTable = selected_batch$peakTable, injections = selected_batch$meta$injection, sampleGroups =
                                     selected_batch$meta$group, QCID = "sQC", RefID = "ltQC", G = seq(1,35, by = 3), reportPath=batchDir)
  sink()

  correctedBatches[[as.integer(batch_no)]] <- drift_correction

  ##### 241028 Hackathon:
  #- Add post-correction plotting of batches
  #- Calle knows person who has done automatic correction choosing algorithm which will apply the most suitable
  #- Implement that instead of batchCorr if library is available

  ##### 241106 Anton: How to project new numbers into old PCA?
  pca_res <- prcomp(drift_correction$TestFeatsFinal, scale. = TRUE)
  png(filename=paste0(batchDir, "/", format(as.Date(Sys.Date()), "%Y%m%d"), "_", chrom, "_", pol,"_PostBatchCorr-Dots.png"), width=2559, height=1376, res=120)
  print(autoplot(pca_res, data=selected_batch$meta, colour="group", size=3)  + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow"))) #autoplot(pca_res, data=meta, colour="Type", size=3, label=T, label.size=3)
  dev.off()

  png(filename=paste0(batchDir, "/", format(as.Date(Sys.Date()), "%Y%m%d"), "_", chrom, "_", pol,"_PostBatchCorr-Labels.png"), width=2559, height=1376, res=120)
  print(autoplot(pca_res, data=selected_batch$meta, colour="group", shape=F, label=T, label.label="filename", label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow")))
  dev.off()

  # rownames(pca_res$x) <- LCMS_data$sample
  png(filename=paste0(batchDir, "/", format(as.Date(Sys.Date()), "%Y%m%d"), "_", chrom, "_", pol,"_PostBatchCorr-Numbers.png"), width=2559, height=1376, res=120)
  print(autoplot(pca_res, data=selected_batch$meta, colour="group", shape=F, label=T, label.label="injection", label.size=3) + scale_color_manual(values = c("black", "red", "green", "orange", "blue", "yellow")))
  dev.off()
}

# Save meta + correctedBatches and transfer PDFs to their respective
saveRDS(LCMS_data, file=paste0(resultsPath3, "meta.rds"))
saveRDS(correctedBatches, paste0(resultsPath3, "corrected.rds"))
save.image(paste0(resultsPath3, format(as.Date(Sys.Date()), "%Y%m%d"), "_Part5.RData"))

#################################
#### Part 6 :: Align batches and batchCorr :: Unified dataset and per batch ####
#################################

#### Cleaning up in-between
rm(list = ls())
gc()

load("paths_params_functions.RData")
setwd(outPath)

# Loading batchCorr objects list
correctedBatches <- readRDS(paste0(resultsPath3, "corrected.rds"))
xcmsData <- readRDS(paste0(resultsPath3, "XCMS_Grouped.rds"))

# Collecting meta data from rownames in FT_fill_imputation
filenames <- rownames(readRDS(paste0(resultsPath3, "FT_fill_imputation.rds")))
batch_codes <- get_batch_codes(filenames)
groups <- get_groups(filenames)
batch_codes_unique <- unique(batch_codes)
batch_codes_unique <- batch_codes_unique[order(as.integer(unique(get_batch_nos(filenames))))]

# Taking care of negative values from batchCorr
flags<-list()
for(i in 1:length(batch_codes_unique)){
  print(paste0("Batch", i))
  print(paste0("Any non-corr < 0: ", any(correctedBatches[[i]]$TestFeats < 0)))
  print(paste0("Any corr < 0: ", any(correctedBatches[[i]]$TestFeatsCorr < 0)))
  print(paste0("Any final == 0: ", any(correctedBatches[[i]]$TestFeatsFinal == 0)))
  print(paste0("Any final < 0: ", any(correctedBatches[[i]]$TestFeatsFinal < 0)))
  
  
  if(any(correctedBatches[[i]]$TestFeatsFinal <= 0)){
    correctedBatches[[i]]$TestFeatsFinal[correctedBatches[[i]]$TestFeatsFinal <= 0] <- NA
    tempDF <- mvImpWrap(MAT = correctedBatches[[i]]$TestFeatsFinal, method = "RF", guess = "minVar", forceZero = TRUE)
    
    if(!is.null(tempDF)){
      correctedBatches[[i]]$TestFeatsFinal <- tempDF
    } else {
      correctedBatches[[i]]$TestFeatsFinal[is.na(correctedBatches[[i]]$TestFeatsFinal)] <- 0
    }
  }
  
  #Flagging features where every value in batch is <= 0-values in whole features
  ##### 241028 Hackathon:
  #- This is kind of sketchy, should be done differently most likely
  #- Fix in future
  zeroProp <- apply(correctedBatches[[i]]$TestFeatsFinal, 2, function(x) sum(x == 0)/length(x))
  if(any(zeroProp == 1)){
    flags[[i]] <- zeroProp[which(zeroProp==1)]
  } else {
    flags[[i]] <- "None"
  }
  
  #Checking if any 0s and replacing with normally distributed low numbers
  if(any(correctedBatches[[i]]$TestFeatsCorr == 0)){
    correctedBatches[[i]]$TestFeatsCorr[correctedBatches[[i]]$TestFeatsCorr == 0] <- 10^rnorm(n = sum(correctedBatches[[i]]$TestFeatsCorr == 0), mean = 0, sd = 1)
  }
  
  ##### 241028 Hackathon:
  #- In dev-branch of batchCorr fix the incorrect feature removal which led to varying dimensions and other things
  # 
  if(dim(correctedBatches[[i]]$TestFeats)[1] == dim(correctedBatches[[i]]$TestFeatsCorr)[1] &&
     dim(correctedBatches[[i]]$TestFeats)[2] == dim(correctedBatches[[i]]$TestFeatsCorr)[2]){
    correctedBatches[[i]]$TestFeatsCorr <- scaleAwayFromZero(PTBefore = correctedBatches[[i]]$TestFeats,
                                                             PTAfter = correctedBatches[[i]]$TestFeatsCorr)
  } else {
    print(paste0("Batch", i, ": DIMS ARE NOT CORRECT!"))
  }
  # 
  # print(any(correctedBatches[[i]]$TestFeatsCorr == 0))
}

#Only aligning and normalizing if there are more than one batch
if(length(correctedBatches) > 1){
  # Merging all batches
  mergedData <- mergeBatches(correctedBatches)
  
  # Making a batch + sample group dataframe
  batches <- get_batch_codes(rownames(mergedData$peakTableCorr))
  types <- get_groups(rownames(mergedData$peakTableCorr))
  
  # Normalizing batch areas
  normData <- normBatch(peakTableCorr = mergedData$peakTableCorr,
                        batches = batches, 
                        sampleGroup = types,
                        refGroup = 'ltQC',
                        population = 'sample')
  PTnorm <- normData$peakTable
} else {
  #If only 1 batch, transfer directly
  PTnorm <- correctedBatches[[1]]$TestFeatsFinal
}

# Tracking which features were removed during merging
trackFeats <- readRDS(paste0(resultsPath3,"trackFeats.rds"))
trackFeats$PostBatchCorr[which(!(trackFeats$Feature %in% colnames(PTnorm)))] <- rep(F, length(which(!(trackFeats$Feature %in% colnames(PTnorm)))))
write.xlsx(trackFeats, paste0(resultsPath4,"trackFeats.xlsx"))
saveRDS(trackFeats, paste0(resultsPath3,"trackFeats.rds"))


# Selecting unique batches and replacing norm-values with NA to impute
for(batch_code in batch_codes_unique){
  i <- which(batch_codes_unique %in% batch_code)
  if(length(flags[[i]]) > 1){
    PTnorm[grepl(batch_code, rownames(PTnorm)),which(colnames(PTnorm) %in% names(flags[[i]]))] <- NA
  } else if (flags[[i]] != "None"){
    PTnorm[grepl(batch_code, rownames(PTnorm)),which(colnames(PTnorm) %in% names(flags[[i]]))] <- NA
  }
}

# Imputing the missing values for features which had entire batch missingness
if(any(is.na(PTnorm))){
  tempDF <- mvImpWrap(MAT = PTnorm, method = "RF", guess = "minVar", forceZero = TRUE)
  if(!is.null(tempDF)){
    PTnorm <- tempDF
  }
}


#Setting all remaining 0s or <0s to 1
if(any(PTnorm <= 1)){
  PTnorm[PTnorm < 1] <- 1
}

# Saving intermediary files
if(length(correctedBatches) > 1){
  saveRDS(normData, paste0(resultsPath4,"normData.rds"))
}
saveRDS(PTnorm, paste0(resultsPath4, "batchCorrFinal.rds"))


############################################
# Setting up ramClust optimization variables
############################################

peakTabFname <- "peakTable.csv"
write.csv(PTnorm,file=peakTabFname)

expDes=defineExperiment(force.skip = T)
sr=c(.2,.3,.4,.5)
st=c(0.1,0.5,1,1.5,2)
maxt=c(5)
par=expand.grid(st=st,sr=sr,maxt=maxt)
str(par)
nClust=nSing=sizeMax=sizeMed=sizeMean=numeric(nrow(par))
nFeat=list()
nSamp <- nrow(PTnorm)
samps=sample(1:nSamp,50) #choosing 20 random samples that will be used during plotting, could be increased for further assurance (at heavy computational cost).
samps <- samps[order(samps)]

# Subsetting xcmsData object to only include the 20 samples which were randomly selected
xcmsData_Files <- filterFile(
  xcmsData,
  samps,
  keepAdjustedRtime = hasAdjustedRtime(xcmsData),
  keepFeatures = TRUE
)

##### 241028 Hackathon:
#- Objective function to determine best parameters
#- Plots might be hard to objectively classify
#- Possible solution
#- Take EICs from clusters and check internal correlation within clusters


##### 241028 Hackathon:
#RRP = Ramclust Reversed Phase
#Rename all variables which are called "RRP" -> "RamClustChromPol"
for (i in 1:nrow(par)) {
  RRP=ramclustR(ms=peakTabFname,
                st = par$st[i],
                sr=par$sr[i],
                maxt = par$maxt[i],
                timepos = 2,
                sampNameCol = 1,
                featdelim = '@',
                ExpDes = expDes) #Important to note that a csv file is needed as input!
  nClust[i]=length(RRP$cmpd)
  nSing[i]=RRP$nsing
  sizeMax[i]=max(RRP$nfeat)
  sizeMed[i]=median(RRP$nfeat)
  sizeMean[i]=mean(RRP$nfeat)
  nFeat[[i]]=RRP$nfeat
  pdf(file=paste0(resultsPath4,'clusts_par',i,'.pdf'),width=15,height=8)
  par(mfrow=c(4,5),mar=c(4,4,2,0)+.5)
  clusts=round(c(2:6,seq(7,max(RRP$featclus),length.out = 15)))
  for (c in clusts) {
    ##### 241028 Hackathon:
    #- Change plotclust to be faster
    #- Take plotEIC function code and add in plotClust
    
    plotClust(ram = RRP,clustnr = c,xcmsData = xcmsData_Files ,samps = samps, dtime=30, dmz=0.1)
  }
  dev.off()
}
paramRows <- cbind(par,nClust,nSing,sizeMax,sizeMean,sizeMed)
write.csv(paramRows, file = paste0(resultsPath4, "paramRows.csv"))

sink(paste0(resultsPath4, "SampsForRamClustEval.txt"))
print(samps)
print(xcmsData_Files)
sink()

saveRDS(xcmsData_Files, paste0(resultsPath4, "xcmsData_Files_RamClust.rds"))

##### 241028 Hackathon:
#- Find representative feature for each cluster and print in a document
#- When clusters established look for:
#- Possible isotopes
#- Diff mz between non-isotopes
#- Suggested mol weight mz for diffs (water, NH4, etc.)
#- Fiehn, phenome center, make a big list and compare against

#### Having determined real final settings for ramClust ####
paramRows <- cbind(par,nClust,nSing,sizeMax,sizeMean,sizeMed)
write.csv(paramRows, file = paste0(resultsPath4, "paramRows.csv"))

finalRamClust=ramclustR(ms=peakTabFname,
                        st = 0.35,
                        sr=0.5,
                        maxt = 5,
                        timepos = 2,
                        sampNameCol = 1,
                        featdelim = '@',
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
rownames(RamDF) <- rownames(PTnorm)
finalRamClust_temp <- finalRamClust
dateToday <- substr(format(as.Date(Sys.Date()), "%Y%m%d"), 3, 1000)
saveRDS(list(finalRamClust_temp, RamDF), file=paste0(resultsPath4, dateToday, '_', chrom, '_', pol, '_Ram.rds'))

##### 241028 Hackathon:
#ChromPol paste checks
write.xlsx(as.data.frame(RamDF), paste0(resultsPath4, dateToday, "_", chrom, "_", pol,"_FinalPT.xlsx"), rowNames=T)
write.csv2(as.data.frame(RamDF), paste0(resultsPath4, dateToday, "_", chrom, "_", pol,"_FinalPT.csv"))
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
    if(any(tempClustDF[j,2] > (tempClustDF[,2] + 1.001) & tempClustDF[j,2] < (tempClustDF[,2] + 1.008))){
      tempClustDF[j,6] <- "Isotope"
    }
    
    #### Checking if potential adduct ####
    #### Anton 241129 - Should probably try all possible combinations ####
    whichAdduct <- which((tempAdductMasses[,2]-1.1) < tempClustDF[j,2]  & (tempAdductMasses[,2]+1.1) > tempClustDF[j,2])
    
    if(length(whichAdduct) > 0){
      whichAdduct <- whichAdduct[which.min(abs(tempAdductMasses[whichAdduct,2] - tempClustDF[j,2]))]
      tempClustDF[which(tempClustDF[,2] == molMz),6] <- "[M+H]+" 
      tempClustDF[j,6] <- tempAdductMasses[whichAdduct,1]
    }
  }
  
  FeatOverview[which(FeatOverview[,1] == as.character(i)),] <- tempClustDF
}


#Saving cluster breakdown data
write.xlsx(as.data.frame(FeatOverview), paste0(resultsPath4, dateToday, "_RP_POS_ClusterBreakdown.xlsx"), rowNames = F)
write.csv2(FeatOverview, paste0(resultsPath4, dateToday, "_RP_POS_ClusterBreakdown.csv"), row.names = F)

#Saving last image
save.image(paste0(resultsPath4, format(as.Date(Sys.Date()), "%Y%m%d"), "_Part6.RData"))

#Saving tracking information
write.xlsx(trackFeats,
           paste0(resultsPath4, dateToday, "_TrackFeats_", chromPol,".xlsx"))

##### 241028 Hackathon:
#Write tracking .xlsx for feature tracking
#Save in SQL format
#Dev