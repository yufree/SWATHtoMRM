##' @title SWATHtoMRM
##' @description An R package to generate MRM transitions from SWATH-MS data.
##' @author Yandong Yin \email{yinyandong@sioc.ac.cn}
##' @param d.in Directory path to load the orgnized data files in mzXML format
##' @param d.out Directory path to save analysis results, if not set, the results
##'  will be saved in 'Results(.x)' folder
##' @param polarity Polarity setup when acquiring SWATH-MS data
##' @param adduct.list Adduct list for CAMERA annotation, which is inherent in
##'  the package
##' @param peakwidth Peak width range in seconds when detecting peaks with xcms
##' @param sn Signal to noise ratio threshold when detecting peaks with xcms
##' @param minfrac Minimal fraction of appearance of peaks when aligning peaks
##'  with xcms
##' @param minfrac.vote Minimal fraction of appearance of fragments when 
##'  generating consensus spectra
##' @param nSlaves Number of threads when detecting peaks with xcms
##' @param ppm.pd Tolerance in ppm when detecting peaks with xcms
##' @param cutoff Peak-Peak Correlation(PPC) threshold to determe if a fragment
##'  is belonging to a precursor ion
##' @param ppm.ms2.mtp Tolerance in ppm of the fragment m/z when extracting MS2
##'  EICs with the fragments in multicomplexed spectra
##' @param int.filter.ms1.field Field of MS1 peaks after peak detection to decide
##'  if the peak shold be kept
##' @param int.filter.ms1 Intensity threshold when filtering MS1 peaks
##' @param int.filter.ms2 Intensity threshold when filtering the fragments in
##' the multiplexed spectra
##' @param isFWHM Logical, if using the full width half maximum of MS1 and MS2
##'  peaks when calculating PPC
##' @param is.plot.eic.feature Logical, if plot the EICs of MS1 peaks
##' @param is.plot.eic.spec Logical, if plot MS2 EICs
##' @param is.ms1.only Logical, if only detect peaks without generating MRM
##'  transitions
##' @param rerun Logical, if re-run all data processing works
SWATHtoMRM <- function(d.in='.', d.out={},
                       polarity = c('positive', 'negative'),
                       adduct.list = c('HILIC', 'RPLC'),
                       peakwidth = c(8, 30),
                       sn = 20, minfrac = 1, minfrac.vote = 0.5,
                       nSlaves = 6, ppm.pd = 15,
                       cutoff = 0.8,
                       ppm.ms2.mtp = 15,
                       int.filter.ms1.field = c('maxo', 'into', 'intb'),
                       int.filter.ms1 = NULL,
                       int.filter.ms2 = 200,
                       isFWHM=FALSE,
                       is.plot.eic.feature = FALSE,
                       is.plot.eic.spec = FALSE,
                       is.ms1.only = FALSE,
                       rerun = FALSE) {
# d.in <- '~/Resilio Sync/Data/Tissue/'
  span <- 0.3
  mz.diff <- 0.01
  # browser()
  wd0 <- getwd()
  setwd(d.in)
  polarity <- match.arg(polarity)
  adduct.list <- match.arg(adduct.list)
  int.filter.ms1.field <- match.arg(int.filter.ms1.field)
  
  files <- list.files(d.in, recursive = TRUE, full.names = TRUE, pattern = '(?i)mzxml$')
  nSlaves <- min(parallel::detectCores() - 1, nSlaves, length(files))

  catLineSeparator('Detecting and aligning features ...')
  fn.skip <- 'diaFeature.RData'
  if ((!rerun) & file.exists(fn.skip)) {
    cat('using existing results:', fn.skip, '...\n')
    load(fn.skip)
  } else {
    dia.feature <- getFeatures(new('diaFeatures'), files,
                               polarity = polarity,
                               adduct.list = adduct.list,
                               ppm.pd = ppm.pd,
                               sn = sn,
                               peakwidth = peakwidth,
                               mzdiff = mz.diff,
                               minfrac = minfrac,
                               is.plot.eic.feature = is.plot.eic.feature,
                               int.filter.ms1 = int.filter.ms1,
                               int.filter.ms1.field = int.filter.ms1.field,
                               nSlaves = nSlaves)
    
    save(dia.feature, file = fn.skip)
  }
  
  if (is.ms1.only) {
    setwd(wd0)
    return()
  }
  
  # browser()
  catLineSeparator('Extracting multicomplexed spectra ...')
  fn.skip <- 'mtpSpecs.RData'
  if ((!rerun) & file.exists(fn.skip)) {
    cat('using existing results:', fn.skip, '...\n')
    load(fn.skip)
  } else {
    mtp.spec <- getMultiplexSpectrum(new('multiplexSpectrum'), files, dia.feature, ppm.ms2.mtp,
                                     int.filter.ms2 = int.filter.ms2, isFWHM = isFWHM,
                                     span = span, rerun = rerun)
    # mtp.spec <- getSimilarityScore(mtp.spec, ppm.ms1, ppm.ms2.mtp, dmz.tol.ppm = dmz.tol.ppm, dmz.tol.abs = dmz.tol.abs, fp.mod = fp.mod, int.filter = int.filter, isFWHM = isFWHM, is.prefilter.dmz = is.prefilter.dmz, span = span)
    
    save(mtp.spec, file = fn.skip)
  }
  # browser()
  if (is.plot.eic.spec) {
    plotSpecEIC(mtp.spec)
  }
  
  catLineSeparator('Generating consensus spectra and MRM transitions ...')
  
  dia.res <- getDIAResult(new('DIAResult'), dia.feature, mtp.spec, cutoff = cutoff, minfrac.vote = minfrac.vote)
  
  res.origin <- dia.res@resultOrigin
  idx.remove <- which(is.na(res.origin[, 'mz.spec']))
  res.origin <- res.origin[-idx.remove, , drop = FALSE]
  outputResult(res.origin, d.out = d.out, fn.out = 'result_origin.csv', rewrite = FALSE, row.names = FALSE)
  outputResult(dia.res@result, d.out = d.out, fn.out = 'result.csv', rewrite = TRUE, row.names = FALSE)
  mrm.trans <- genMRMTrans(dia.res@result)
  outputResult(mrm.trans, d.out = d.out, fn.out = 'MRMTransition.csv', rewrite = TRUE, row.names = FALSE)
  cat('All processing work done!!')
}