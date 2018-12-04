setGeneric('catLineSeparator', function(info) {
  info.sep <- paste(rep('-', 38), collapse  = '')
  cat('', info.sep, info, info.sep, '', sep = '\n')
})

setGeneric('loadData', function(file, keep.name = FALSE, env){
  if (missing( env)) env <- new.env()
  b <- load(file, env = env)
  if (keep.name | length( b) > 1) {
    r <- lapply( b, function(b1) env[[ b1]])
    names( r) <- b
    r
  } else {
    env[[ b]]
  }
})

setGeneric('outputResult', function(dt.matrix, fn.out = {}, d.out = {}, rewrite = FALSE, row.names = TRUE) {
  if (is.null(d.out)) {
    d.out <- d.out0 <- 'Results'
    i <- 1
    d.out.norewrite <- d.out
    while (file.exists(d.out)) {
      i <- i + 1
      d.out <- paste(d.out0, i, sep = '.')
      d.out.norewrite <- ifelse(i == 2, d.out0, paste(d.out0, i - 1, sep = '.'))
    }
  } else{
    d.out0 <- d.out.norewrite <- d.out
    i <- 1
    while (file.exists(d.out)) {
      i <- i + 1
      d.out <- paste(d.out0, i, sep = '.')
      d.out.norewrite <- ifelse(i == 2, d.out0, paste(d.out0, i - 1, sep = '.'))
    }
  }
  
  if (rewrite) {
    d.out <- d.out.norewrite
  }
  
  if (!file.exists(d.out)) {
    dir.create(d.out)
  }
  
  write.csv(dt.matrix, file = file.path(d.out, fn.out), row.names = row.names)
})