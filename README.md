# HRF

Calculate hemodynamic response function (HRF) from simultaneously recorded neuronal activities and hemodynamic changes.

 

# get_photometry

Read raw spectra, timestamps and corresponding wavelengths from the Oceanview output .txt file into Matlab.

 

# interleaving_fix_withBG

1. Read 400/488 nm interleaved recorded GCaMP spectra time series output from Oceanview.

2. Resolve 400 and 488 excited GCaMP spectra, and fix lost-frames if there's any.

3. Remove the autofluorescence of optical cable.

4. Calculate the amplitude time courses of GCaMP signal excited by 400 and 488 nm lasers respectively.

5. Extract hemoglobin changes time course from GCaMP signal excited by 400 nm laser.

 

# unmixing2

Read and segregate the mixed raw spectra data from the Oceanview output .txt file using user provided reference spectra.

 

# wavelet_TFmap

Generate time-frequency map from GCaMP signal time course.
