# get the entire peak shapes (intensity) with time (rt) for each detected peakgroup
setMethod('getMultiplexSpectrum', 'multiplexSpectrum', function(obj, files, dFeatures, ppm.ms2.mtp, int.filter.ms2 = 100, isFWHM=TRUE, span = 0.3, rerun = FALSE) {
  # browser()
  cat('Processing raw data ... \n')
  smpNames <- gsub('(?i)\\.mzxml$' , '', basename(files))
  smp.xRaw <- sapply(files, function(fn) {
    xcmsRaw(fn, includeMSn = TRUE)
  })
  names(smp.xRaw) <- smpNames
  
  cat('initializing dData ... \n')
  if ((!rerun) & file.exists('dData.RData')) {
    load('dData.RData')
  } else {
    dData <- sapply(smpNames, function(smp) {
      diaData(new('diaData'), smp.xRaw[[smp]])
    })
    
    dData <- sapply(dData, function(dD) {
      checkMissingScans(dD, dData)
    })
    # browser()
    save(dData, file = 'dData.RData')
  }
  
  cat('extracting multiplex spectra ... \n')
  # scanIndex: peak index of matched precursor ions with scan time
  scanIdx <- lapply(seq_along(dFeatures@peakgroup$corrected), function(idx.ft){
    lapply(seq_along(dData), function(idx.smp) {
      smp <- dData[[idx.smp]]
      which(dData[[idx.smp]]@ms1$scantime >= dFeatures@peakgroup$raw[[idx.ft]][idx.smp, 'rtmin'] &
              dData[[idx.smp]]@ms1$scantime <= dFeatures@peakgroup$raw[[idx.ft]][idx.smp, 'rtmax'])
    })
  })
  names(scanIdx) <- names(dData)
  
  featureMZ <- lapply(dFeatures@peakgroup$raw, function(pg) {
    mz <- pg[, 'mz']
    names(mz) <- names(dData)
    mz
  })
  
  # get precursor(ms1) peaks
  featurePeaks <- lapply(seq_along(dFeatures@peakgroup$raw), function(idx.pg) {
    pg.raw <- dFeatures@peakgroup$raw[[idx.pg]]
    # pg.corrected <- dFeatures@peakgroup$corrected[[idx.pg]]
    peaks.all <- lapply(seq(nrow(pg.raw)), function(idx.smp) {
      idx <- scanIdx[[idx.pg]][[idx.smp]]
      mz <- featureMZ[[idx.pg]][[idx.smp]]
      # mzRange <- ppmRange(mz, ppm.ms1)
      mzRange <- as.numeric(pg.raw[idx.smp, c('mzmin', 'mzmax')])
      scan.ms1 <- dData[[idx.smp]]@ms1$scan
      peaks.ms1 <- precursorMZ(scan.ms1, idx, mz, mzRange)
    })
    names(peaks.all) <- names(dData)
    peaks.all
  })
  # browser()
  
  # get product(ms2) spectrums (mtpSpectrums) under each feature of different scan time
  mtpSpectrums.scanIdx.ms2 <- lapply(seq_along(dFeatures@peakgroup$raw), function(idx.pg) {
    pg.raw <- dFeatures@peakgroup$raw[[idx.pg]]
    # pg.corrected <- dFeatures@peakgroup$corrected[[idx.pg]]
    scans.all <- lapply(seq(nrow(pg.raw)), function(idx.smp) {
      idx.ms1 <- scanIdx[[idx.pg]][[idx.smp]]
      mz <- featureMZ[[idx.pg]][idx.smp]
      mzRange <- dData[[idx.smp]]@ms1$mzrange
      info.SWATH <- getSWATHinfo(length(dData[[idx.smp]]@ms1$scan), length(dData[[idx.smp]]@ms2$scan), dData[[idx.smp]]@ms1$mzrange)
      idx.SWATH <- getIntervalIndex(mz, info.SWATH$window)
      idx.ms2 <- (idx.ms1 - 1) * (info.SWATH$nSWATH) + idx.SWATH
    })
    names(scans.all) <- names(dData)
    scans.all
  })
  # browser()
  
  obj@featureMZ <- featureMZ
  obj@featurePeaks <- featurePeaks
  obj@mtpSpectrums <- mtpSpectrums.scanIdx.ms2
  
  # dmz.ptm <- as.numeric(as.matrix(read.csv(fp.mod, stringsAsFactors = F))[,2])
  
  # smoothing ms1 peak with loess
  features.smooth <- lapply(seq_along(obj@featureMZ), function(idx.ft) {
    ft.smooth <- lapply(seq_along(obj@featureMZ[[idx.ft]]), function(idx.smp) {
      n.spec <- nrow(obj@featurePeaks[[idx.ft]][[idx.smp]])
      if (n.spec > 3) {
        smoothLoess(obj@featurePeaks[[idx.ft]][[idx.smp]], span = span)
      } else
        obj@featurePeaks[[idx.ft]][[idx.smp]]
    })
    names(ft.smooth) <- names(obj@featureMZ[[1]])
    ft.smooth
  })
  # browser()
  # match ms2 spectrums with lib (refering to mtp in highest int of smoothed ms1 peak) mz of spec
  match.spectrums <- lapply(seq_along(obj@featureMZ), function(idx.ft) {
    # idx.ft <- 628
    match.smps <- lapply(seq_along(obj@featureMZ[[idx.ft]]), function(idx.smp) {
      scan.ms2 <- dData[[idx.smp]]@ms2$scan
      spec.exp <- fragmentMZ(scan.ms2, mtpSpectrums.scanIdx.ms2[[idx.ft]][[idx.smp]])
      # spec.exp <- obj@mtpSpectrums[[idx.ft]][[idx.smp]]
      idx.lib.mtp <- which.max(features.smooth[[idx.ft]][[idx.smp]][,2])
      idx.lib.mtp.intfilter <- which(spec.exp[[idx.lib.mtp]][,2] >= int.filter.ms2)
      precursor.mz <- obj@featureMZ[[idx.ft]][[idx.smp]]
      lib.mz <- spec.exp[[idx.lib.mtp]][idx.lib.mtp.intfilter,1]
      idx.lib.mtp.mzfilter <- idx.lib.mtp.intfilter[lib.mz <= precursor.mz + 25e-6 * max(400, precursor.mz)]
      lib.mz.mtp.mzfilter <- spec.exp[[idx.lib.mtp]][idx.lib.mtp.mzfilter,1]
      lib.mz.diff <- diff(lib.mz.mtp.mzfilter)
      is.filter.isotope <- TRUE
      
      while(is.filter.isotope) {
        idx.remove <- which(lib.mz.diff >= 1.0003 & lib.mz.diff <= 1.0063) + 1
        if(length(idx.remove) > 0) {
          lib.mz.mtp.mzfilter <- lib.mz.mtp.mzfilter[-idx.remove]
          lib.mz.diff <- diff(lib.mz.mtp.mzfilter)
        } else {
          lib.mz.mtp.final <- lib.mz.mtp.mzfilter
          is.filter.isotope <- FALSE
        }
      }
      
      lib.mzRange <- t(sapply(lib.mz.mtp.final, function(mz) {
        ppmRange(mz, ppm.ms2.mtp)
      }))
      
      lapply(obj@featureMZ[[idx.ft]][[idx.smp]], function(idx.lib) {
        if (length(lib.mz.mtp.final) > 0) {
          spectra.pre <- searchSpectra(spec.exp, lib.mzRange, lib.mz.mtp.final)
          spectra <- lapply(spectra.pre, function(spec.peak) {
            do.call(cbind, spec.peak)
          })
        } else {
          spectra <- list()
        }
        spectra
      })
    })
    names(match.smps) <- names(obj@featureMZ[[1]])
    match.smps
  })
  
  match.spectrums.smooth <- lapply(seq_along(match.spectrums), function(idx.ft) {
    match.smooth <- lapply(seq_along(match.spectrums[[idx.ft]]), function(idx.smp) {
      if (!is.null(match.spectrums[[idx.ft]][[idx.smp]])) {
        lapply(match.spectrums[[idx.ft]][[idx.smp]], function(specs) {
          lapply(specs, function(spec) {
            if (nrow(spec) > 3) {
              smoothLoess(spec, span = span)
            } else {
              spec
            }
          })
        })
      }
    })
    names(match.smooth) <- names(obj@featureMZ[[1]])
    match.smooth
  })
  
  score.feature <- getFeatureSimilarity(match.spectrums.smooth, features.smooth, isFWHM)
  # score.spectrum <- getSpectrumSimilarity(match.spectrums.smooth, lib@spectrums, match.features.smooth, score.feature)
  # score.accuMS <- getAccurateMSscore(obj@featureMZ, obj@featureMZ, lib@precursorMZ, ppm)
  # browser()
  
  # obj@featureMZ <- featureMZ
  # obj@featurePeaks <- featurePeaks
  # obj@mtpSpectrums <- mtpSpectrums.scanIdx.ms2
  
  # return(obj)
  obj@featureMatch <- obj@featurePeaks
  obj@featureMatchSmooth <- features.smooth
  obj@spectrumMatch <- match.spectrums
  obj@featureScore <- score.feature
  # obj@spectrumScore <- score.spectrum
  # obj@accurateMSScore <- score.accuMS
  return(obj)
})


setMethod('plotSpecEIC', 'multiplexSpectrum', function(obj) {
  # browser()
  samplenames <- names(obj@featureMatch[[1]])
  dirSpecEIC <- file.path('SpecEICs', samplenames)
  names(dirSpecEIC) <- samplenames
  sapply(dirSpecEIC, function(fp) {
    if (!file.exists(fp)) {
      dir.create(fp, recursive = TRUE)
    }
  })
  
  sapply(seq_along(obj@featureMatch), function(idx.ft) {
    ft <- obj@featurePeaks[[idx.ft]]
    spec.ft <- obj@spectrumMatch[[idx.ft]]
    sapply(samplenames, function(nm.smp) {
      spec.smp <- spec.ft[[nm.smp]][[1]]
      if (!is.null(spec.smp)) {
        # nm.id <- names(spec.smp)
        # sapply(seq_along(nm.id), function(idx.id) {
        #   spec.id <- spec.smp[[idx.id]]
        fn.smp <- file.path(dirSpecEIC[nm.smp], paste0(paste('featureIdx', idx.ft, sep = '_'), '.png'))
        clr <- rainbow(length(spec.smp), alpha = 1)
        
        png(file = fn.smp, width = 960, height = 480)
        # par(mfrow = c(1,2))
        int.max <- max(max(do.call(rbind, spec.smp)[,2]), ft[[nm.smp]][,2])
        plot(ft[[nm.smp]][,2], ylim = c(0, int.max), type = 'b', pch = 19,cex = 1, main = "EIC", ylab = "intensity", xlab = "index")
        sapply(seq_along(spec.smp), function(idx.spec) {
          points(spec.smp[[idx.spec]][,2], type = 'b', cex = 0.8, col = clr[idx.spec])
          lines(spec.smp[[idx.spec]][,2], type = 'b', cex = 0.8, col = clr[idx.spec])
        })
        dev.off()
      }
    })
  })
})

# find the precursor peak shape according to scans
setGeneric('precursorMZ', function(scan.ms1, scan.idx, mz, mzRange) {
  mz <- unname(mz)
  peak <- c()
  for (idx in scan.idx) {
    scan.ms1.single <- scan.ms1[[idx]]
    idx.match <- which(scan.ms1.single[, 'mz'] >= mzRange[1] & scan.ms1.single[, 'mz'] <= mzRange[2])
    peak <- rbind(peak, switch(as.character(length(idx.match)),
                               "0" = c('mz' = mz, 'intensity' = 0),
                               '1' = scan.ms1.single[idx.match,],
                               scan.ms1.single[idx.match[which.max(scan.ms1.single[idx.match,'intensity'])],]
    )
    )
  }
  return(peak)
})

# get ms2 data for each scan with given scan indexes
setGeneric('fragmentMZ', function(scan.ms2, idx.ms2) {
  fragmentMZ <- lapply(idx.ms2, function(idx) {
    scan.ms2[[idx]]
  })
  names(fragmentMZ) <- as.character(idx.ms2)
  return(fragmentMZ)
})

# output a data range with given reference data and ppm
setGeneric('ppmRange', function(ref, ppm, spliter = 400) {
  dev <- ppm * 1e-6
  c(ref - max(dev*spliter/2, dev*ref/2), ref + max(dev*spliter/2, dev*ref/2))
})

# get info for SWATH setup
setGeneric('getSWATHinfo', function(lenMS1 = {}, lenMS2 = {}, mzRange = {}, fp = '../SWATHsetup.csv') {
  infoSWATH <- list()
  if (file.exists(fp)) {
    swath.window <- read.csv(fp)
    nSWATH <- nrow(swath.window)
    overlapped <- 1 + which(swath.window[-1, 1] - swath.window[-nSWATH, 2] < 0)
    swath.window[overlapped, 1] <- swath.window[overlapped - 1, 2]
  } else {
    nSWATH <- lenMS2/lenMS1
    swath.window.list <- (mzRange[2] - mzRange[1])/nSWATH * seq(0,nSWATH) + mzRange[1]
    swath.window.from <- swath.window.list[-(nSWATH + 1)]
    swath.window.to <- swath.window.list[-1]
    swath.window <- cbind(swath.window.from, swath.window.to)
  }
  infoSWATH$nSWATH <- nSWATH
  infoSWATH$window <- swath.window
  return(infoSWATH)
})

# get the index of the given data located in a range defined by row of matrix
setGeneric('getIntervalIndex', function(dataPoint, dataMatrix) {
  isDone <- FALSE
  idx <- 1
  while (!isDone) {
    if ((dataPoint >= dataMatrix[idx, 1] & dataPoint < dataMatrix[idx, 2]) | idx == length(dataMatrix[,2])) {
      isDone <- TRUE
    } else {
      idx <- idx + 1
    }
  }
  return(idx)
})

# smoothing peaks with loess method
# setGeneric('smoothLoess', function(preData, span = 0.3, degree = 1) {
setGeneric('smoothLoess', function(preData, span = 0.3, degree = 1L) {
  pre.inf <- data.frame(idx = seq(nrow(preData)), intensity = preData[,2])
  pre.fit <- loess(intensity~idx, data = pre.inf, span = span, degree = degree)
  pre.predict <- predict(pre.fit, data.frame(idx = seq(nrow(preData))))
  cbind(mz = preData[,1], intensity = pre.predict)
})

# match library data with mz range
setGeneric('searchLib', function(mz, lib) {
  idx.matched <- NULL
  nms <- names(lib)
  
  for (idx in seq_along(lib)) {
    mz.lib <- lib[[idx]]
    if (mz >= mz.lib[1] & mz <= mz.lib[2]) {
      idx.matched <- c(idx.matched, idx)
    }
  }
  lib.matched <- idx.matched
  if (!is.null(idx.matched)) {
    names(lib.matched) <- nms[lib.matched]
  }
  
  return(lib.matched)
})
#
# setGeneric('searchSpectrums', function(spec.exp, spec.lib.mzRange, spec.lib.mz) {
#   match.spectrum.peak <- lapply(seq_along(spec.lib.mzRange), function(idx.lib) {
#     lib.mzRange <- spec.lib.mzRange[[idx.lib]]
#     spec.matched <- c()
#     for (idx in seq_along(spec.exp)) {
#       spec.exp.single <- spec.exp[[idx]]
#       match.result <- c()
#       for (nr in seq(nrow(spec.exp.single))) {
#         mz <- spec.exp.single[nr, 'mz']
#         if (mz >= lib.mzRange[1] & mz <= lib.mzRange[2]) {
#           match.result <- rbind(match.result, spec.exp.single[nr, ])
#         }
#       }
#
#       if (is.null(match.result)) {
#         spec.matched <- rbind(spec.matched, matrix(c(spec.lib.mz[idx.lib], 0), ncol = 2))
#       } else {
#         # select the one with closest mz compared to spectrum mz in lib
#         spec.matched <- rbind(spec.matched, match.result[which.min(abs(match.result[,'mz'] - spec.lib.mz[idx.lib])), ])
#       }
#     }
#
#     if (!is.null(spec.matched)) {
#       colnames(spec.matched) <- c('mz', 'intensity')
#     }
#     spec.matched
#   })
#
#   return(match.spectrum.peak)
# })

setGeneric('getFeatureSimilarity', function(spectrums, features, isFWHM) {
  # using FWHM: Full Width at Half Maximum
  samplenames <- names(features[[1]])
  # browser()
  score <- lapply(seq_along(spectrums), function(idx.ft) {
    sc.smp <- lapply(seq_along(features[[idx.ft]]), function(idx.smp) {
      if (isFWHM) {
        feature <- features[[idx.ft]][[idx.smp]]
        # FWHM: Full Width at Half Maximum
        if (!all(feature[, 2] == 0)) {
          range.feature.FWHM <- range(which(feature[,2] > max(feature[,2])/2))
          idx.feature.FWHM <- range.feature.FWHM[1]:range.feature.FWHM[2]
          feature.FWHM <- feature[idx.feature.FWHM, , drop = FALSE]
          names.match <- names(spectrums[[idx.ft]][[idx.smp]])
        }
        sc.match <- lapply(spectrums[[idx.ft]][[idx.smp]], function(spec.matches) {
          if (!all(feature[, 2] == 0)) {
            # spec.matches <- spectrums[[idx.ft]][[idx.smp]][[nm.lib]]
            score.pearson <- NULL
            for (idx.spec in seq_along(spec.matches)) {
              spec <- spec.matches[[idx.spec]]
              spec.FWHM <- spec[idx.feature.FWHM,, drop = FALSE]
              sc <- pearsonSimilarity(spec.FWHM[, 'intensity'], feature.FWHM[, 'intensity'])
              if (is.na(sc)) {
                sc <- 0
              }
              score.pearson[idx.spec] <- sc
            }
            score.pearson
          } else {
            rep(0, length(spec.matches))
          }
        })
      } else {
        feature <- features[[idx.ft]][[idx.smp]]
        # if (!all(feature[, 2] == 0)) {
        #   range.feature.FWHM <- range(which(feature[,2]>max(feature[,2])/2))
        #   idx.feature.FWHM <- range.feature.FWHM[1]:range.feature.FWHM[2]
        #   feature.FWHM <- feature[idx.feature.FWHM, , drop=FALSE]
        #   names.match <- names(spectrums[[idx.ft]][[idx.smp]])
        # }
        sc.match <- lapply(spectrums[[idx.ft]][[idx.smp]], function(spec.matches) {
          if (!all(feature[, 2] == 0)) {
            # spec.matches <- spectrums[[idx.ft]][[idx.smp]][[nm.lib]]
            score.pearson <- NULL
            for (idx.spec in seq_along(spec.matches)) {
              spec <- spec.matches[[idx.spec]]
              # spec.FWHM <- spec[idx.feature.FWHM,, drop=FALSE]
              sc <- pearsonSimilarity(spec[, 'intensity'], feature[, 'intensity'])
              if (is.na(sc)) {
                sc <- 0
              }
              score.pearson[idx.spec] <- sc
            }
            score.pearson
          } else {
            rep(0, length(spec.matches))
          }
        })
      }
      
      sc.match
    })
    names(sc.smp) <- samplenames
    sc.smp
  })
})


setGeneric('getSpectrumSimilarity', function(spectrums, lib.spectrum, features, score.feature, cutoff.ftSim = 0.6) {
  samplenames <- names(features[[1]])
  # browser()
  score <- lapply(seq_along(spectrums), function(idx.ft) {
    sc.smp <- lapply(seq_along(features[[idx.ft]]), function(idx.smp) {
      feature <- features[[idx.ft]][[idx.smp]]
      if (!is.null(feature)) {
        idx.max.int <- unname(which.max(feature[,'intensity']))
        names.match <- names(spectrums[[idx.ft]][[idx.smp]])
        sc.match <- lapply(names.match, function(nm.lib) {
          spec.matches <- spectrums[[idx.ft]][[idx.smp]][[nm.lib]]
          spec.exp <- c()
          for (idx.spec in seq_along(spec.matches)) {
            spec.match <- spec.matches[[idx.spec]][idx.max.int, ]
            #             if (score.feature[[idx.ft]][[idx.smp]][[nm.lib]][idx.spec] <= cutoff.ftSim) {
            #               spec.match['intensity'] <- 0
            #             }
            spec.exp <- rbind(spec.exp, spec.match)
          }
          spec.lib <- lib.spectrum[[nm.lib]]
          score.cos <- cosineSimilarity(spec.exp[, 'intensity'], spec.lib[, 'LibraryIntensity'])
        })
        names(sc.match) <- names.match
        sc.match
      }
    })
    names(sc.smp) <- samplenames
    sc.smp
  })
})

setGeneric('getAccurateMSscore', function(match.feature, mz.feature, mz.lib.list, ppm) {
  samplenames <- names(match.feature[[1]])
  score <- lapply(seq_along(match.feature), function(idx.ft) {
    sc.smp <- lapply(seq_along(match.feature[[idx.ft]]), function(idx.smp) {
      match.smp <- match.feature[[idx.ft]][[idx.smp]]
      if (!is.null(match.smp)) {
        names.match <- names(match.smp)
        mz.exp <- mz.feature[[idx.ft]][[idx.smp]]
        sc.match <- lapply(names.match, function(nm.lib) {
          mz.lib <- mz.lib.list[[nm.lib]]
          score.accuMS <- accurateMSscore(mz.exp, mz.lib, ppm)
        })
        names(sc.match) <- names.match
        sc.match
      }
    })
    names(sc.smp) <- samplenames
    sc.smp
  })
})

setGeneric('pearsonSimilarity', function(int.ms2, int.ms1) {
  if (length(int.ms2) < 3) {
    return(NA)
  }
  score.pearson <- cor(int.ms2, int.ms1, method = 'pearson')
  return(score.pearson)
})

setGeneric('cosineSimilarity', function(int.exp, int.lib) {
  score.cos <- sum(int.exp * int.lib) / sqrt(sum(int.exp ^ 2) * sum(int.lib ^ 2))
  # browser()
  return(ifelse(is.nan(score.cos), 0, score.cos))
})

setGeneric('accurateMSscore', function(mz.exp, mz.lib, ppm, spliter = 400) {
  if (mz.exp < spliter) {
    dm <- ppm * 1e-6 * spliter
  } else {
    dm <- ppm * 1e-6 * mz.lib
  }
  score.accuMS <- exp(-(mz.exp - mz.lib) ^ 2 / (2*dm ^ 2))
  return(score.accuMS)
})

