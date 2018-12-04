setMethod('diaData', 'diaData', function(obj, xRaw) {
  
  # extracting MS1 scans
  obj@ms1$scantime <- xRaw@scantime
  obj@ms1$mzrange <- xRaw@mzrange
  obj@ms1$scan <- lapply(seq_along(obj@ms1$scantime) , function(idx) {
    getScan(xRaw, idx)
  })
  # browser()
  # extracting MS2 scans
  xRaw.ms2 <- msn2xcmsRaw(xRaw)
  obj@ms2$scantime <- xRaw.ms2@scantime
  obj@ms2$mzrange <- xRaw.ms2@mzrange
  obj@ms2$scan <- lapply(seq_along(obj@ms2$scantime) , function(idx) {
    getScan(xRaw.ms2, idx)
  })
  
  return(obj)
})

setMethod('checkMissingScans', 'diaData', function(obj, dData) {
  scan.num.ms1 <- sapply(dData, function(dD) {
    length(dD@ms1$scantime)
  })
  max.idx <- which.max(scan.num.ms1)
  max.scantime <- dData[[max.idx]]@ms1$scantime
  max.scan <- dData[[max.idx]]@ms1$scan
  # browser()
  if (length(dData) == 1) {
    scantime.ms1 <- obj@ms1$scantime
    scantime.ms2 <- obj@ms2$scantime
    is.miss.mid <- any(round(diff(scantime.ms1)/median(diff(scantime.ms1)), 0) > 1)
    is.miss.start <- scantime.ms1[1] > scantime.ms2[1]
    is.miss.end <- round(scantime.ms2[length(scantime.ms2)] - scantime.ms1[length(scantime.ms1)], 0) > round(median(diff(scantime.ms1)), 0)
    if (is.miss.start | is.miss.start | is.miss.end) {
      cat('Some of the MS1 scans may be missing. Trying to complete them ... \n')
      obj <- completeMS1Scans(obj, max.scantime, max.scan, 'single')
    }
  } else {
    if (length(obj@ms1$scantime) != max(scan.num.ms1)) {
      cat('Some of the MS1 scans may be missing. Trying to complete them ... \n')
      obj <- completeMS1Scans(obj, max.scantime, max.scan, 'multiple')
    }
  }
  
  # checking missing MS2 scans
  len.ms1 <- length(obj@ms1$scan)
  len.ms2 <- length(obj@ms2$scan)
  if(len.ms2 / len.ms1 != round(len.ms2 / len.ms1)) {
    cat('Some of the MS2 scans may be missing. Trying to complete them ... \n')
    obj <- completeMS2Scans(obj)
  }
  return(obj)
})

setMethod('completeMS1Scans', 'diaData', function(obj, max.scantime, max.scan, smpInfo) {
  switch(smpInfo,
         'multiple' = {
           scantime <- obj@ms1$scantime
           
           idx.exist <- sapply(scantime, function(sct) {
             which.min(abs(sct - max.scantime))
           })
           mtx.fill <- matrix(0, ncol = 2)
           colnames(mtx.fill) <- c('mz', 'intensity')
           max.scantime[idx.exist] <- scantime
           max.scan[idx.exist] <- obj@ms1$scan
           max.scan[-idx.exist] <- list(mtx.fill)
           
           obj@ms1$scantime <- max.scantime
           obj@ms1$scan <- max.scan
         },
         'single' = {
           scan.ms1 <- obj@ms1$scan
           scantime.ms1 <- obj@ms1$scantime
           scantime.ms2 <- obj@ms2$scantime
           
           mtx.fill <- matrix(0, ncol = 2)
           colnames(mtx.fill) <- c('mz', 'intensity')
           
           if(scantime.ms1[1] > scantime.ms2[1]) {
             scantime.ms1[2:(length(scantime.ms1)+1)] <- scantime.ms1
             scantime.ms1[1] <- median(diff(scantime.ms2)) * 1.6
             scan.ms1 <- append(list(mtx.fill), scan.ms1)
           }
           
           if(round(scantime.ms2[length(scantime.ms2)] - scantime.ms1[length(scantime.ms1)], 0) > round(median(diff(scantime.ms1)), 0)) {
             scantime.ms1[length(scantime.ms1+1)] <- scantime.ms2[length(scantime.ms2)] + median(diff(scantime.ms2)) * 2
             scan.ms1 <- append(scan.ms1, list(mtx.fill))
           }
           
           if(any(round(diff(scantime.ms1)/median(diff(scantime.ms1)), 0) > 1)) {
             md1 <- median(diff(scantime.ms1))
             dfmd1 <- round(diff(scantime.ms1) / md1, 0)
             idx.miss <- which(dfmd1 > 1)
             
             scantime.ms1.theo <- seq(sum(dfmd1), from = 0, by = md1)
             idx.exist <- sapply(scantime.ms1, function(t){
               which.min(abs(t - scantime.ms1.theo))
             })
             scantime.ms1.theo[idx.exist] <-scantime.ms1
             scantime.ms1 <- scantime.ms1.theo
             
             scan.ms1.empty <- rep(list(mtx.fill), sum(dfmd1))
             scan.ms1.empty[idx.exist] <- scan.ms1
             scan.ms1 <- scan.ms1.empty
           }
           
           obj@ms1$scantime <- scantime.ms1
           obj@ms1$scan <- scan.ms1
         })
  return(obj)
})

setMethod('completeMS2Scans', 'diaData', function(obj) {
  # browser()
  scantime.ms1 <- obj@ms1$scantime
  scantime.ms2 <- obj@ms2$scantime
  idx.remove <- which(scantime.ms2 <= scantime.ms1[1])
  scan.ms2 <- obj@ms2$scan[-idx.remove]
  scantime.ms2 <- scantime.ms2[-idx.remove]
  
  ms2.intervial <- median(diff(scantime.ms2))
  ms2.fold <- ceiling(length(scantime.ms2) / length(scantime.ms1))
  
  scantime.ms2.theoritical.list <- lapply(1:length(scantime.ms1), function(y){
    scantime.ms12mix.1group <- seq(from = scantime.ms1[y], length.out = ms2.fold + 1, by = ms2.intervial)
    scantime.ms2.1group <- scantime.ms12mix.1group[-1]
  })
  scantime.ms2.theoritical <- do.call(c, scantime.ms2.theoritical.list)
  
  names(scantime.ms2.theoritical) <- c(1 : length(scantime.ms2.theoritical))
  idx.ms2.theoritical <- sapply(scantime.ms2, function(y){
    tmp <- y - scantime.ms2.theoritical
    which.min(tmp[tmp >= 0])
  })
  
  spec.non <- matrix(0, ncol = 2)
  colnames(spec.non) <- c('mz', 'intensity')
  
  names(scan.ms2) <- names(scantime.ms2)[idx.ms2.theoritical]
  scan.ms.complete <- list()
  scan.ms.complete <- sapply(1:length(scantime.ms2.theoritical), function(y){
    scan.ms.complete <- list(spec.non)
  })
  
  scantime.ms2.theoritical[idx.ms2.theoritical] <- scantime.ms2
  scan.ms.complete[idx.ms2.theoritical] <- scan.ms2
  obj@ms2$scantime <- unname(scantime.ms2.theoritical)
  obj@ms2$scan <- scan.ms.complete
  # browser()
  return(obj)
})

# setMethod('completeMS2Scans', 'diaData', function(obj) {
#   scantime.ms1 <- obj@ms1$scantime
#   scantime.ms2 <- obj@ms2$scantime
#   scan.ms2 <- obj@ms2$scan
#
#   ms2.intervial <- median(diff(scantime.ms2))
#   ms2.fold <- ceiling(length(scantime.ms2) / length(scantime.ms1))
#
#   scantime.ms2.theoritical.list <- lapply(1:length(scantime.ms1), function(y){
#     scantime.ms12mix.1group <- seq(from = scantime.ms1[y], length.out = ms2.fold + 1, by = ms2.intervial)
#     scantime.ms2.1group <- scantime.ms12mix.1group[-1]
#   })
#   scantime.ms2.theoritical <- do.call(c, scantime.ms2.theoritical.list)
#
#   names(scantime.ms2.theoritical) <- c(1 : length(scantime.ms2.theoritical))
#   idx.ms2.theoritical <- sapply(scantime.ms2, function(y){
#     tmp <- y - scantime.ms2
#     which.min(tmp[tmp >= 0])
#   })
#
#   spec.non <- matrix(0, ncol = 2)
#   colnames(spec.non) <- c('mz', 'intensity')
#
#   names(scan.ms2) <- names(scantime.ms2)[idx.ms2.theoritical]
#   scan.ms.complete <- list()
#   scan.ms.complete <- sapply(1:length(scantime.ms2.theoritical), function(y){
#     scan.ms.complete <- list(spec.non)
#   })
#
#   scantime.ms2.theoritical[idx.ms2.theoritical] <- scantime.ms2
#   scan.ms.complete[idx.ms2.theoritical] <- scan.ms2
#   obj@ms2$scantime <- unname(scantime.ms2.theoritical)
#   obj@ms2$scan <- scan.ms.complete
#   return(obj)
# })
