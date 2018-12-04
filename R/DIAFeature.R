setMethod('getFeatures', 'diaFeatures', function(obj, files, polarity = c('positive', 'negative'), adduct.list = c('HILIC', 'RPLC'), int.filter.ms1 = NULL, int.filter.ms1.field = c('maxo', 'into', 'intb'), ppm.pd = 15, sn = 6, peakwidth = c(5,30), mzdiff = 0.01, minfrac = 0.5, is.plot.eic.feature = TRUE, nSlaves = 6) {
  is.xcms3 <- FALSE
  if (packageVersion("xcms") >= '1.50.1') {
    is.xcms3 <- TRUE
    require(BiocParallel)
    param <- SnowParam(workers = nSlaves, type = 'SOCK')
  }
  
  int.filter.ms1.field <- match.arg(int.filter.ms1.field)
  fn.adduct.list <- switch(polarity,
                           'positive' = paste0(adduct.list, '_POS.csv'),
                           'negative' = paste0(adduct.list, '_NEG.csv')
  )
  rules.camera <- read.csv(file.path(system.file(package = getPackageName(), 'rules'), fn.adduct.list))
  if (is.xcms3) {
    xset <- xcmsSet(files, method = 'centWave', ppm = ppm.pd, snthr = sn, peakwidth = peakwidth, mzdiff = mzdiff, BPPARAM = param)
  } else {
    xset <- xcmsSet(files, method = 'centWave', ppm = ppm.pd, snthr = sn, peakwidth = peakwidth, mzdiff = mzdiff, nSlaves = nSlaves)
  }
  require(CAMERA)
  if (length(files) == 1) {
    # for one sample only, no peak alignment
    xa <- xsAnnotate(xset, polarity = polarity)
    xa <- groupFWHM(xa)
    xa <- findIsotopes(xa)
    xa <- findAdducts(xa, rules = rules.camera, polarity = polarity)
    peaklist.filtered <- peaklist.anno <- getPeaklist(xa, intval = int.filter.ms1.field)
    nr.isotope <- which(grepl('\\]\\[M\\+\\d+', peaklist.anno[, 'isotopes']))
    write.csv(peaklist.anno, 'PeakTable-annotated_origin.csv',
              row.names = FALSE)
    
    nr.low.int <- which(xset@peaks[, int.filter.ms1.field] < int.filter.ms1)
    if (length(nr.low.int) > 0 | length(nr.isotope) > 0) {
      nr.remove <- unique(c(nr.isotope, nr.low.int))
      xset@peaks <- xset@peaks[-nr.remove, , drop = FALSE]
      peaklist.filtered <- peaklist.anno[-nr.remove, , drop = FALSE]
    }
    peaklist.filtered <- cbind('ft.idx' = seq(nrow(peaklist.filtered)), peaklist.filtered)
    write.csv(peaklist.filtered, 'PeakTable-annotated_filtered.csv',
              row.names = FALSE)
    
    tmp <- lapply(seq(nrow(xset@peaks)), function(idx.r) {
      r <- xset@peaks[idx.r, ]
      r.gp <- c(r[1:6], 1, 1)
      names(r.gp) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks", "10")
      r.gp
    })
    # browser()
    xset@groups <- do.call(rbind,tmp)
    xset@groupidx <- lapply(seq(nrow(xset@groups)), function(i) {
      i
    })
    
    sclassv <- as.character(sampclass(xset))
    for (i in 1:length( xset@filepaths)) {
      cat( i, xset@filepaths[ i], "\t [",  sclassv[ i], "]", sep = "")
      cat( "  --> ", length( which( xset@peaks[, "sample"] == i)), " deisotoped features. \n")
    }
    
    
    # obj@peaks <- xset@peaks[-nr.isotope, , drop = FALSE]
    obj@peaks <- peaklist.filtered
    obj@peakgroup$raw <- lapply(seq(nrow(obj@peaks)),
                                function(nr) {
                                  obj@peaks[nr, , drop = FALSE]
                                })
    obj@peakgroup$corrected <- obj@peakgroup$raw
    obj@samplenames <- gsub('(?i)\\.mzxml$' , '', basename(files))
    ftname <- sapply(seq_along(obj@peakgroup$raw), function(idx) {
      paste("M", sprintf('%.0f', obj@peakgroup$raw[[idx]][, 'mz']),
            "T", sprintf('%.0f', obj@peakgroup$raw[[idx]][, 'rt']),
            sep = "")
    })
    is.duplicated <- TRUE
    while (is.duplicated) {
      i <- 2
      idx.duplicated <- which(duplicated(ftname))
      if (length(idx.duplicated) > 0) {
        ftname[idx.duplicated] <- paste0(ftname[idx.duplicated], '_', i)
        i <- i + 1
      } else {
        is.duplicated <- FALSE
      }
    }
    
    obj@featurenames <- matrix(ftname, nrow = 1)
    obj@groupnames <- plotFeatureEIC(xset, obj@featurenames, nr.isotope, is.plot.eic.feature = is.plot.eic.feature, extra = 30)[]
  } else {
    # for multiple samples, peak alignment
    # browser()
    
    xset <- group(xset, minfrac = minfrac)
    sclassv <- as.character(sampclass( xset))
    for (i in 1:length( xset@filepaths)) {
      cat( i, xset@filepaths[ i], "\t [",  sclassv[ i], "]", sep = "")
      cat( "  --> ", length( which( xset@peaks[, "sample"] == i)), " Features. \n")
    }
    
    pdf('retcor-obiwarp.pdf')
    xset2 <- retcor(xset, method = 'obiwarp', plottype = 'deviation', profStep = 0.1)
    dev.off()
    
    # xset2 <- group(xset2, bw = 10)
    xset2 <- group(xset2, bw = 10, mzwid = 0.015, minfrac = minfrac)
    
    groupmat <- xcms::groups( xset2)
    pt <- data.frame( cbind( groupmat, groupval( xset2, "medret", "into")), row.names = NULL)
    write.csv( pt, file = "PeakTable-raw.csv", row.names = T)
    gc()
    
    cat( nrow( groupmat), "aligned features. \n")
    
    xset3 <- fillPeaks(xset2)
    
    xa <- xsAnnotate(xset3, polarity= polarity, nSlaves = nSlaves)
    xa <- groupFWHM(xa)
    xa <- findIsotopes(xa)
    xa <- findAdducts(xa, rules = rules.camera, polarity = polarity)
    peaklist.filtered <- peaklist.anno <- getPeaklist(xa, intval = int.filter.ms1.field)
    nr.isotope <- which(grepl('\\]\\[M\\+\\d+', peaklist.anno[, 'isotopes']))
    # peaklist.anno <- cbind('ft.idx' = seq(nrow(peaklist.anno)), peaklist.anno)
    write.csv(peaklist.anno, 'PeakTable-annotated.csv', row.names = FALSE)
 
    ms1.int <- unname(apply(groupval(xset3, value = int.filter.ms1.field), 1, max))
    nr.low.int <- which(ms1.int < int.filter.ms1)
    if (length(nr.low.int) > 0 | length(nr.isotope) > 0) {
      nr.remove <- unique(c(nr.isotope, nr.low.int))
      xset3@groups <- xset3@groups[-nr.remove, , drop = FALSE]
      xset3@groupidx <- xset3@groupidx[-nr.remove]
      peaklist.filtered <- peaklist.anno[-nr.remove, , drop = FALSE]
    }
    peaklist.filtered <- cbind('ft.idx' = seq(nrow(peaklist.filtered)), peaklist.filtered)
    write.csv(peaklist.filtered, 'PeakTable-annotated_filtered.csv',
              row.names = FALSE)
    
    obj@peaks <- data.frame(xset3@peaks)
    
    grpIdx.pre <- xset3@groupidx
    grpIdx <- lapply(seq_along(grpIdx.pre), function(idx) {
      Index <- grpIdx.pre[[idx]][order(abs(obj@peaks[grpIdx.pre[[idx]], 'rt'] - median(obj@peaks[grpIdx.pre[[idx]], 'rt'])))]
      Index[match(seq_along(sampnames(xset)), obj@peaks[Index,'sample'])]
    })
    
    obj@peakgroup$corrected <- lapply(grpIdx, function(idx) {
      obj@peaks[idx, , drop = FALSE]
    })
    
    obj@rt$corrected <- xset3@rt$corrected
    obj@rt$raw <- xset3@rt$raw
    
    obj@peakgroup$raw <- lapply(obj@peakgroup$corrected, function(peak) {
      peak.raw.list <- lapply(1:nrow(peak), function(y) {
        rt.idx <- which.min(abs(peak[y, "rt"] - obj@rt$corrected[[peak[y, "sample"]]]))
        peak[y, "rt"] <- obj@rt$raw[[peak[y, "sample"]]][rt.idx]
        rtmin.idx <- which.min(abs(peak[y, "rtmin"] - obj@rt$corrected[[peak[y, "sample"]]]))
        peak[y, "rtmin"] <- obj@rt$raw[[peak[y, "sample"]]][rtmin.idx]
        rtmax.idx <- which.min(abs(peak[y, "rtmax"] - obj@rt$corrected[[peak[y, "sample"]]]))
        peak[y, "rtmax"] <- obj@rt$raw[[peak[y, "sample"]]][rtmax.idx]
        peak[y, ]
      })
      peak.raw <- do.call(rbind, peak.raw.list)
    })
    
    obj@samplenames <- sampnames(xset3)
    obj@featurenames <- sapply(seq_along(obj@peakgroup$raw), function(idx) {
      paste("M", sprintf('%.0f', obj@peakgroup$raw[[idx]][, 'mz']), "T", sprintf('%.0f', obj@peakgroup$raw[[idx]][, 'rt']), sep = "")
    })
    obj@groupnames <- plotFeatureEIC(xset3, obj@featurenames, is.plot.eic.feature = is.plot.eic.feature, extra = 30)
  }
  return(obj)
})

setGeneric('plotFeatureEIC', function(xset, featurenames.smp, is.plot.eic.feature = FALSE, checkEIC = TRUE, extra = 30) {
  rt.range.data <- range(xset@rt)
  
  rtmin.origin <- unname(apply(groupval(xset, value = 'rtmin'), 1, median))
  rtmax.origin <- unname(apply(groupval(xset, value = 'rtmax'), 1, median))
  
  rtmin <- rtmin.origin - extra
  rtmax <- rtmax.origin + extra
  
  rtmin[which(rtmin < rt.range.data[1])] <- rt.range.data[1]
  rtmax[which(rtmax > rt.range.data[2])] <- rt.range.data[2]
  
  rt.range.eic <- cbind( rtmin, rtmax)
  eics <- getEIC(xset, rtrange = rt.range.eic, sampleidx = seq(length(xset@filepaths)), groupidx = seq(nrow(rt.range.eic)))
  
  if (is.plot.eic.feature) {
    cat('Ploting EICs ... ')
    smpnames <- names(eics@eic)
    clr <- rainbow(length(smpnames), end = 0.85)
    rgbvec <- pmin(col2rgb(clr) + 153, 255)
    clr_gray <- rgb(rgbvec[1, ], rgbvec[2, ], rgbvec[3, ], max = 255)
    
    
    dirEics <- file.path('featureEICs', c('all', smpnames))
    names(dirEics) <- c('all', smpnames)
    sapply(dirEics, function(fp) {
      if (!file.exists(fp)) {
        dir.create(fp, recursive = TRUE)
      }
    })
    
    fn.all <- file.path(dirEics['all'], paste0(eics@groupnames, '.png'))
    # browser()
    
    if (length(smpnames) > 1) {
      sapply(seq_along(eics@groupnames), function(idx.ft) {
        
        max.int <- max(unname(sapply(smpnames, function(smpname) {
          max(eics@eic[[smpname]][[idx.ft]][, 2])
        })))
        png(file = fn.all[idx.ft], width = 960, height = 480)
        plot(0, 0, type = "n", xlim = rt.range.eic[idx.ft,], ylim = c(0, max.int),
             xlab = "Retention Time (seconds)", ylab = "Intensity",
        )
        sapply(seq_along(smpnames), function(idx.smp) {
          ft <- eics@eic[[smpnames[idx.smp]]][[idx.ft]]
          lines(ft, col = clr_gray[idx.smp])
          idx.ft.peak <- which(ft[, 1] >= rtmin.origin[idx.ft] & ft[, 1] <= rtmax.origin[idx.ft])
          lines(ft[idx.ft.peak, ], col = clr[idx.smp])
          if (checkEIC) {
            abline(v = c(rtmin.origin[idx.ft], rtmax.origin[idx.ft]), col = clr[idx.smp])
          }
        })
        legend('topright', smpnames, col = clr, lty = 1)
        title(eics@groupnames[idx.ft])
        par(new = FALSE)
        dev.off()
      })
    }
    
    sapply(seq_along(smpnames), function(idx.smp) {
      sapply(seq_along(eics@groupnames), function(idx.ft) {
        png(file = file.path(dirEics[smpnames[idx.smp]], paste0(featurenames.smp[idx.smp,idx.ft], '.png')), width = 960, height = 480)
        ft <- eics@eic[[smpnames[idx.smp]]][[idx.ft]]
        plot(ft, col = 'black')
        lines(ft, col = clr_gray[1])
        idx.ft.peak <- which(ft[, 1] >= rtmin.origin[idx.ft] & ft[, 1] <= rtmax.origin[idx.ft])
        lines(ft[idx.ft.peak, ], col = clr[1])
        if (checkEIC) {
          abline(v = c(rtmin.origin[idx.ft], rtmax.origin[idx.ft]), col = clr[1])
        }
        title(featurenames.smp[idx.smp,idx.ft])
        dev.off()
      })
    })
  }
  return(eics@groupnames)
})
