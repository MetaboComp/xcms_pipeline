#### normBatch
# Normalizing batches using rownames instead of supplied vectors

normBatch <- function (peakTableCorr, peakTableOrg, batches, sampleGroup, 
                       refGroup = "QC", population = "all", CVlimit = 0.3, FCLimit = 5, 
                       medianZero = "min") {
  if (missing(peakTableOrg)) 
    peakTableOrg <- peakTableCorr
  if (!(identical(dim(peakTableCorr), dim(peakTableOrg)) & 
        identical(colnames(peakTableCorr), colnames(peakTableOrg)))) 
    stop("Mismatch between peakTableCorr and peakTableOrg")
  if (missing(medianZero)) 
    medianZero <- "min"
  if (population == "all") 
    popSample <- rep(TRUE, nrow(peakTableCorr))
  else {
    if (!population %in% sampleGroup) 
      (stop("population identifier needs to be present in sampleGroups\nConsider setting population=\"all\""))
    else {
      popSample <- sampleGroup == population
    }
  }
  peakTableNormalized <- peakTableCorr
  uniqBatch <- unique(batches)
  nBatch <- length(uniqBatch) #nBatch <- 4
  nFeat <- ncol(peakTableCorr)
  CVMatrix = matrix(nrow = nBatch, ncol = nFeat)
  rownames(CVMatrix) <- uniqBatch
  colnames(CVMatrix) <- colnames(peakTableCorr)
  RefMeanIntensity <- CVMatrix
  RefNormMatrix <- matrix(FALSE, nrow = nBatch, ncol = nFeat)
  MeanIntensityRatios <- matrix(1, nrow = nBatch, ncol = nBatch)
  for (b in 1:nBatch) {
    batch <- uniqBatch[b]
    peakTableBatch <- peakTableCorr[grepl(batch, rownames(peakTableCorr)), ]
    peakTableBatch <- peakTableBatch[grepl(refGroup, rownames(peakTableBatch)), ]
    CVMatrix[b, ] = ifelse(batchCorr:::cv(peakTableBatch) <= CVlimit, 
                           TRUE, FALSE)
    RefMeanIntensity[b, ] = apply(peakTableBatch, 2, mean)
  }
  for (b in 1:(nBatch - 1)) {
    for (bb in (b + 1):nBatch) {
      MeanIntensityRatios[bb, b] = mean(RefMeanIntensity[bb, 
      ])/mean(RefMeanIntensity[b, ])
      MeanIntensityRatios[b, bb] = 1/MeanIntensityRatios[bb, 
                                                         b]
    }
  }
  for (feat in 1:nFeat) { #826
    featureIntensityRatios = matrix(1, nrow = nBatch, ncol = nBatch)
    for (b in 1:(nBatch - 1)) {
      for (bb in (b + 1):nBatch) {
        featureIntensityRatios[bb, b] = RefMeanIntensity[bb, 
                                                         feat]/RefMeanIntensity[b, feat]
        featureIntensityRatios[b, bb] = 1/featureIntensityRatios[bb, 
                                                                 b]
      }
    }
    candidates = abs(log(featureIntensityRatios/MeanIntensityRatios)) <= 
      log(FCLimit)
    
    # if(any(is.nan(abs(log(featureIntensityRatios/MeanIntensityRatios))))){
    #   print(paste0("Feat: ", feat))
    #   print(abs(log(featureIntensityRatios/MeanIntensityRatios)))
    #   print(featureIntensityRatios)
    #   print(MeanIntensityRatios)
    #   
    # }
    
    candidates[is.na(candidates)] = FALSE
    for (b in 1:nBatch) {
      if (CVMatrix[b, feat] == FALSE | is.na(CVMatrix[b, 
                                                      feat])) {
        candidates[, b] = candidates[b, ] = FALSE
      }
    }
    refBatch = min(which(colSums(candidates) == max(colSums(candidates))))
    refCorrFlags = candidates[, refBatch]
    refCorrIndex = which(refCorrFlags)
    refCorrIndex = refCorrIndex[refCorrIndex != refBatch]
    refIntensity = RefMeanIntensity[refBatch, feat]
    for (b in refCorrIndex) {
      correctionFactor = refIntensity/RefMeanIntensity[b, 
                                                       feat]
      peakTableNormalized[grepl(uniqBatch[b], rownames(peakTableCorr)), feat] = peakTableCorr[grepl(uniqBatch[b], rownames(peakTableCorr)), feat] * correctionFactor
    }
    RefNormMatrix[, feat] = refCorrFlags
    if (length(refCorrIndex) + 1 != nBatch) {
      refCorrected = peakTableNormalized[batches %in% 
                                           uniqBatch[c(refBatch, refCorrIndex)] & popSample, 
                                         feat]
      refCorrMedian = median(refCorrected)
      WhichPOPCorr = which(!refCorrFlags)
      for (n in WhichPOPCorr) {
        populationMedian <- median(peakTableOrg[grepl(uniqBatch[n], rownames(peakTableOrg)) & popSample, feat])
        if (populationMedian == 0) {
          if (medianZero == "min") {
            populationMedian <- peakTableOrg[grepl(uniqBatch[n], rownames(peakTableOrg)) & popSample, feat]
            populationMedian <- min(populationMedian[populationMedian != 
                                                       0])
          }
          else if (medianZero == "mean") {
            populationMedian <- mean(peakTableOrg[grepl(uniqBatch[n], rownames(peakTableOrg)) & popSample, feat])
          }
          else stop("Other options not included at present.")
        }
        correctionFactor = refCorrMedian/populationMedian
        peakTableNormalized[grepl(uniqBatch[n], rownames(peakTableNormalized)), 
                            feat] = peakTableOrg[grepl(uniqBatch[n], rownames(peakTableOrg)), 
                                                 feat] * correctionFactor
      }
    }
  }
  return(list(peakTable = peakTableNormalized, refCorrected = RefNormMatrix))
}