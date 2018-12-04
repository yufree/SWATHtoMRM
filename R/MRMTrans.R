setGeneric('genMRMTrans', function(dt.spec) {
  # browser()
  dt.pk <- read.csv('PeakTable-annotated_filtered.csv',
                        stringsAsFactors = FALSE)
  idx.ft <- dt.pk[, 'ft.idx']
  mz.diff.spare <- cbind('x' = c(17.0235, 18.0076, 43.9869),
                         'y' = c(17.0295, 18.0136, 43.9929))
  rownames(mz.diff.spare) <- c('NH3', 'H2O', 'CO2')
  
  FilterSpec <- function(dt.spec, mz.diff.spare, mz.diff.thr = 14.0126) {
    colnames(dt.spec)
    # browser()
    mz.diff <- dt.spec[, 'mz.precursor'] - dt.spec[, 'mz.product']
    nr.rm.thr <- which(mz.diff <= 14.0126)
    nr.rm.spare <- unlist(apply(mz.diff.spare, 1, function(mz.diff.spare.1) {
      which(mz.diff >= mz.diff.spare.1[1] & mz.diff <= mz.diff.spare.1[2])
    }))
    nr.rm.all <- unique(c(nr.rm.thr, nr.rm.spare))
    return(dt.spec[-nr.rm.all, , drop = FALSE])
  }
  
  dt.spec.filtered <- FilterSpec(dt.spec, mz.diff.spare)
  
  idx.ft.spec <- unique(dt.spec.filtered[, 'ft.idx'])

  trans.spec <- lapply(idx.ft.spec, function(idx) {
    info <- dt.spec.filtered[dt.spec.filtered[, 'ft.idx'] == idx, , drop = FALSE]
    nr.max <- which.max(info[, 'int.product'])
    mz.precursor <- info[nr.max, 'mz.precursor']
    mz.product <- info[nr.max, 'mz.product']
    rt <- info[nr.max, 'rt']
    c(idx, mz.precursor, mz.product, round(rt/60, 2), round(rt, 1),
      round(mz.precursor, 1), round(mz.product, 1)
    )
  })
  dt.trans <- do.call(rbind, trans.spec)
  colnames(dt.trans) <- c('ft.idx', 'mz.precursor', 'mz.product', 'rt(min)', 'rt(s)',
                          'Q1', 'Q3')
  return(dt.trans)
})

