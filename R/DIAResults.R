setMethod('getDIAResult', 'DIAResult', function(obj, dFeatures, mtpSpec, cutoff = 0.8, minfrac.vote = 0.5) {
  nm.smp <- dFeatures@samplenames
  num.smp <- length(nm.smp)
  sc.na <- matrix(rep(NA, num.smp), ncol = num.smp)
  colnames(sc.na) <- nm.smp
  # ft.name <- dFeatures@featurenames

  info.anno <- read.csv('PeakTable-annotated_filtered.csv',
                        stringsAsFactors = FALSE)
  # browser()
  result.detail <- lapply(seq_along(mtpSpec@featureMatch), function(ft.idx) {
    # ft.idx <- 532
    # cat(ft.idx, '\n')
    sc.ft <- mtpSpec@featureScore[[ft.idx]]
    # sc.spec <- mtpSpec@spectrumScore[[ft.idx]]
    # sc.accuMS <- mtpSpec@accurateMSScore[[ft.idx]]
    ft.info <- dFeatures@peakgroup$raw[[ft.idx]]
    rownames(ft.info) <- NULL
    anno.col <- c('ft.idx', 'isotopes', 'adduct', 'pcgroup')
    anno.info <- do.call(rbind, replicate(nrow(ft.info), info.anno[ft.idx, anno.col], simplify = FALSE))
    
    ft.info <- cbind(ft.info, anno.info)
    
    ft.info.ncol <- ncol(ft.info)
    ft.info.colname <- colnames(ft.info)
    ft.group.name <- dFeatures@groupnames[ft.idx]
    
    # p.id <- table(unname(unlist(sapply(sc.ft, names)))) / num.smp
    idx.match.ft.spec <- lapply(mtpSpec@featureScore[[ft.idx]], function(sc.ids) {
      lapply(sc.ids, function(sc) {
        which(sc >= cutoff)
      })
    })
    # ms2.mz <- lapply(mtpSpec@spectrumMatch[[ft.idx]][[1]]
    
    match.ft.spec.list <- lapply(nm.smp, function(smp) {
      smp.idx <- which(nm.smp == smp)
      # browser()
      if (length(idx.match.ft.spec[[smp]][[1]]) > 0) {
        score.match.ft.spec <- lapply(seq_along(mtpSpec@featureScore[[ft.idx]][[smp]]), function(idx.sc){
          sc <- mtpSpec@featureScore[[ft.idx]][[smp]][[idx.sc]]
          sc[idx.match.ft.spec[[smp]][[1]]]
        })
        
        mz.match.ft.spec <- lapply(seq_along(mtpSpec@spectrumMatch[[ft.idx]][[smp]]), function(idx.specs){
          sapply(idx.match.ft.spec[[smp]][[idx.specs]], function(idx.spec) {
            median(mtpSpec@spectrumMatch[[ft.idx]][[smp]][[idx.specs]][[idx.spec]][,1])
          })
        })
        
        idx.max <- which.max(mtpSpec@featureMatchSmooth[[ft.idx]][[smp]][,'intensity'])
        
        int.match.ft.spec <- lapply(seq_along(mtpSpec@spectrumMatch[[ft.idx]][[smp]]), function(idx.specs){
          sapply(idx.match.ft.spec[[smp]][[idx.specs]], function(idx.spec) {
            unname(mtpSpec@spectrumMatch[[ft.idx]][[smp]][[idx.specs]][[idx.spec]][idx.max,2])
          })
        })
        
        ms2.info.mz <- cbind(do.call(cbind, mz.match.ft.spec),do.call(cbind, int.match.ft.spec), do.call(cbind,score.match.ft.spec))
        colnames(ms2.info.mz) <- c('mz.spec', 'int.spec', 'score.PPC')
        
        ms2.info.mz.nrow <- nrow(ms2.info.mz)
        
        # rep(ft.info[smp.idx, ], 3)
        
        ms1.info.rep <- do.call(rbind, replicate(ms2.info.mz.nrow, ft.info[smp.idx, ], simplify = FALSE))
        
        gpname.rep <- rep(ft.group.name, ms2.info.mz.nrow)
        ft.idx.rep <- rep(ft.idx, ms2.info.mz.nrow)
        # browser()
        tmp <- data.frame('ft.idx' = ft.idx.rep,
                          'nm.gp' = gpname.rep,
                          ms2.info.mz,
                          ms1.info.rep, stringsAsFactors = FALSE)
      } else {
        ms2.info.mz <- matrix(rep(NA, 3), ncol = 3)
        colnames(ms2.info.mz) <- c('mz.spec', 'int.spec', 'score.PPC')
        # browser()
        
        ms1.info.rep <- ft.info
        tmp <- data.frame('ft.idx' = ft.idx,
                          'nm.gp' = ft.group.name,
                          ms2.info.mz,
                          as.matrix(ft.info[smp.idx,]),
                          stringsAsFactors = FALSE)
        
      }
    })
    id.matrix <- do.call(rbind, match.ft.spec.list)
    rownames(id.matrix) <- NULL
    id.matrix
  })
  result.detail.table <- do.call(rbind, result.detail)
  nr.remove <- which(is.na(result.detail.table[, 'mz.spec']))
  result.detail.final <- result.detail.table[-nr.remove, , drop = FALSE]
  # browser()
  ft.idx.unique <- unique(result.detail.final[, 'ft.idx'])
  result.list <- lapply(ft.idx.unique, function(ft.idx) {
    spec <- lapply(seq_along(nm.smp), function(smp.idx) {
      # nr <- which(info.anno[, 'ft.idx'] == ft.idx)
      nr.origin <- which(result.detail.final[, 'ft.idx'] == ft.idx &
                           result.detail.final[, 'sample'] == smp.idx)
      spec.smp <- result.detail.final[nr.origin, c('mz.spec', 'int.spec'), drop = FALSE]
      colnames(spec.smp) <- c('mz', 'intensity')
      spec.smp
    })
    spec.concensus <- GetConcensusSpec(spec, minfrac.vote = minfrac.vote)
    if (nrow(spec.concensus) == 0) {
      return(NULL)
    }
    nr <- which(info.anno[, 'ft.idx'] == ft.idx)
    # nr.origin <- which(result.detail.final[, 'ft.idx'] == ft.idx)
    info.precursor <- do.call(c, info.anno[nr, c('mz', 'rt', make.names(nm.smp))])
    ft.rep <- matrix(rep(c(ft.idx, info.precursor), nrow(spec.concensus)),
                     ncol = 3+length(nm.smp), byrow = TRUE)
    spec.combind <- cbind(ft.rep, spec.concensus)
    colnames(spec.combind) <- c('ft.idx', 'mz.precursor', 'rt', make.names(nm.smp), 'mz.product', 'int.product')
    spec.combind
  })
  result.list <- result.list[which(!sapply(result.list, is.null))]
  result <- do.call(rbind, result.list)
  obj@resultOrigin <- result.detail.table
  obj@result <- result
  return(obj)
})

setGeneric('RemoveRingEffect', function(spec, mz.diff.thr = 0.3, int.rel.thr = 0.2) {
  nr.ring <- nrow(spec) + 1
  mz <- spec[, 1]
  
  mz.diff <- diff(mz)
  idx.mzdiff <- which(mz.diff <= mz.diff.thr)
  if (length(idx.mzdiff) == 0) {
    return(spec)
  }
  nr.ring.possible <- unique(c(idx.mzdiff, idx.mzdiff + 1))
  
  while (TRUE) {
    
    idx.int.max <- which.max(spec[nr.ring.possible, 2])
    nr.int.max <- nr.ring.possible[idx.int.max]
    int.thr <- spec[nr.int.max, 2] * int.rel.thr
    
    mz.diff <- abs(mz[nr.ring.possible[-idx.int.max]] - mz[nr.int.max])
    int <- spec[nr.ring.possible[-idx.int.max], 2]
    nr.ring <- append(nr.ring, nr.ring.possible[which(mz.diff <= mz.diff.thr & int <= int.thr)])
    nr.ring.possible <- nr.ring.possible[!nr.ring.possible %in% c(nr.ring, nr.int.max)]
    if (length(nr.ring.possible) == 0) {
      break
    }
  }
  
  return(spec[-nr.ring, , drop = FALSE])
})

setGeneric('GetConcensusSpec', function(ms2.all, minfrac.vote = 0.5){
  # merge fragments within the same m/z bin across spectra
  ms2.merged <- MergeFragments(ms2.all)
  # quanlity control: removing ring effect and low intensity spectra
  ms2.merged <- RemoveRingEffect(ms2.merged)
  # get concensus spectra with peak voting
  ms2.final <- VoteSpectra(ms2.merged, ms2.all, minfrac.vote = minfrac.vote)
})

setGeneric('MergeFragments', function(ms2.all, ppm = 30, mz.bin.min = 0.004) {
  spec.all <- do.call(rbind, ms2.all)
  spec.all <- spec.all[order(spec.all[, 'mz']), , drop = FALSE]
  
  idx.left <- seq(nrow(spec.all))
  
  spec.merged <- {}
  while (length(idx.left) > 0 ) {
    idx <- tail(idx.left, 1)
    mz <- spec.all[idx, 'mz']
    mz.range <- c(-1, 0) * max(prod(mz, ppm, 1e-6), mz.bin.min) + mz
    idx.range <- idx.left[spec.all[idx.left, 'mz'] >= mz.range[1] &
                            spec.all[idx.left, 'mz'] <= mz.range[2]]
    spec.tmp <- sapply(c('mz', 'intensity'), function(x) {
      quantile(spec.all[idx.left[idx.range], x], 0.5)
    })
    spec.merged <- rbind(spec.merged, spec.tmp)
    idx.left <- idx.left[-idx.range]
  }
  colnames(spec.merged) <- c('mz', 'intensity')
  rownames(spec.merged) <- NULL
  spec.merged <- spec.merged[order(spec.merged[, 'mz']), , drop = FALSE]
  
  return(spec.merged)
})

setGeneric('RemoveRingEffect', function(spec, mz.diff.thr = 0.3, int.rel.thr = 0.2) {
  spec <- spec[order(spec[, 'mz']), , drop = FALSE]
  nr.ring <- nrow(spec) + 1
  mz <- spec[, 'mz']
  
  mz.diff <- diff(mz)
  idx.mzdiff <- which(mz.diff <= mz.diff.thr)
  if (length(idx.mzdiff) == 0) {
    return(spec)
  }
  
  nr.ring.possible <- unique(c(idx.mzdiff, idx.mzdiff + 1))
  while (TRUE) {
    idx.int.max <- which.max(spec[nr.ring.possible, 2])
    nr.int.max <- nr.ring.possible[idx.int.max]
    int.thr <- spec[nr.int.max, 2] * int.rel.thr
    
    mz.diff <- abs(mz[nr.ring.possible[-idx.int.max]] - mz[nr.int.max])
    int <- spec[nr.ring.possible[-idx.int.max], 2]
    nr.ring <- append(nr.ring, nr.ring.possible[-idx.int.max][which(mz.diff <= mz.diff.thr & int <= int.thr)])
    nr.ring.possible <- nr.ring.possible[!nr.ring.possible %in% c(nr.ring, nr.int.max)]
    if (length(nr.ring.possible) == 0) {
      break
    }
  }
  
  return(spec[-nr.ring, , drop = FALSE])
})

setGeneric('VoteSpectra', function(ms2.merged, ms2.all,
                                   minfrac.vote = 0.5,
                                   ppm = 30,
                                   mz.bin.min = 0.004) {
  num.contained <- sapply(ms2.merged[, 'mz'], function(mz) {
    mz.range <- mz + max(prod(mz, ppm, 1e-6), mz.bin.min) * c(-1, 1)
    is.contained <- sapply(ms2.all, function(spec) {
      any(spec[, 'mz'] >= mz.range[1] & spec[, 'mz'] <= mz.range[2])
    })
    sum(is.contained)
  })
  
  ms2.merged[num.contained/length(ms2.all) >= minfrac.vote, , drop = FALSE]
})
