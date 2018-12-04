setClass('libSpectrum', slots = list(precursorMZ = 'list', rt = 'matrix', spectrums = 'list', metaName = 'matrix', name = 'character'))

setClass('multiplexSpectrum', slots = list(mtpSpectrums = 'list', featurePeaks = 'list', featureMatchSmooth = 'list', featureMZ = 'list', featureMatch = 'list', spectrumMatch = 'list', featureScore = 'list', spectrumScore = 'list', accurateMSScore = 'list'))

setClass('diaData', slots = list(ms1 = 'list', ms2 = 'list'))

setClass('diaFeatures', slots = list(peaks = 'data.frame', peakgroup = 'list', rt = 'list', samplenames = 'character', featurenames = 'matrix', groupnames = 'character'))

setClass('DIAResult', slots = list(result = 'matrix', resultOrigin = 'data.frame'))
