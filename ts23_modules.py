from pyraf import iraf
import os
import shutil
import numpy as np

#########################################################################################################
########################################### IRAF modules ################################################
#########################################################################################################

################################################################################
#give the readnoise and gain values
#check these values in the header of your spectra
readnoise = 3.06
gain      = 0.584

#if you want to call a function
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.onedspec()
iraf.twodspec()
iraf.apextract()
iraf.echelle()

#defining apertures
iraf.echelle.setParam('dispaxi', 1)       #horizontal dispersion in TS23

################## organizing data #######################
def edit_ccdproc(image, output):
    iraf.ccdproc.setParam('images', image)
#    iraf.ccdproc.setParam('output', image)    #leaving ouput empty will overwrite the image, be careful!!
    iraf.ccdproc.setParam('ccdtype', '')
    iraf.ccdproc.setParam('max_cac', 0)
    iraf.ccdproc.setParam('noproc', 'no')
    iraf.ccdproc.setParam('fixpix', 'yes')
    iraf.ccdproc.setParam('overscan', 'yes')
    iraf.ccdproc.setParam('trim', 'yes')
    iraf.ccdproc.setParam('zerocor', 'no')
    iraf.ccdproc.setParam('darkcor', 'no')
    iraf.ccdproc.setParam('flatcor', 'no')
    iraf.ccdproc.setParam('illumco', 'no')
    iraf.ccdproc.setParam('fringec', 'no')
    iraf.ccdproc.setParam('readcor', 'no')
    iraf.ccdproc.setParam('scancor', 'no')
    iraf.ccdproc.setParam('readaxi', 'line')	
    iraf.ccdproc.setParam('fixfile', 'bad_pixels.dat')	
    iraf.ccdproc.setParam('biassec', '[2053:2080,2:2047]')   #Ivan's prescription
    iraf.ccdproc.setParam('trimsec', '[10:2040,2:2047]')     #Ivan's prescription 
    iraf.ccdproc.setParam('zero', '')
    iraf.ccdproc.setParam('dark', '')
    iraf.ccdproc.setParam('flat', '')
    iraf.ccdproc.setParam('illum', '')
    iraf.ccdproc.setParam('fringe', '')
    iraf.ccdproc.setParam('minrepl', 1.)
    iraf.ccdproc.setParam('scantyp', 'shortscan')
    iraf.ccdproc.setParam('nscan', 1)
    iraf.ccdproc.setParam('interac', 'no')   #change to yes if you need to check the fit
    iraf.ccdproc.setParam('functio', 'legendre')
    iraf.ccdproc.setParam('order', 5)         # this order works very well
    iraf.ccdproc.setParam('sample', '*')
    iraf.ccdproc.setParam('naverag', 1)
    iraf.ccdproc.setParam('niterat', 5)
    iraf.ccdproc.setParam('low_rej', 3.)
    iraf.ccdproc.setParam('high_rej', 3.)
    iraf.ccdproc.setParam('grow', 0.)
    iraf.ccdproc.setParam('mode', 'q1')
    iraf.ccdproc(image=image)

def edit_zerocombine(input):
    iraf.zerocombine.setParam('output', 'zero.fits')
    iraf.zerocombine.setParam('combine', 'average')
    iraf.zerocombine.setParam('reject', 'minmax')
    iraf.zerocombine.setParam('ccdtype', '')
    iraf.zerocombine.setParam('process', 'no')
    iraf.zerocombine.setParam('delete', 'no')
    iraf.zerocombine.setParam('clobber', 'no')
    iraf.zerocombine.setParam('scale', 'none')
    iraf.zerocombine.setParam('statsec', '')
    iraf.zerocombine.setParam('nlow', 0)
    iraf.zerocombine.setParam('nhigh', 1)
    iraf.zerocombine.setParam('nkeep', 1)
    iraf.zerocombine.setParam('mclip', 'yes')
    iraf.zerocombine.setParam('lsigma', 3.)
    iraf.zerocombine.setParam('hsigma', 3.)
    iraf.zerocombine.setParam('rdnoise', readnoise)
    iraf.zerocombine.setParam('gain', gain)
    iraf.zerocombine.setParam('snoise', 0.)
    iraf.zerocombine.setParam('pclip', -0.5)
    iraf.zerocombine.setParam('blank', 0.)
    iraf.zerocombine.setParam('mode', 'q1')
    iraf.zerocombine(input=input)

#to apply zero correction
def edit_zerocorrection(image):
    iraf.ccdproc.setParam('images', image)
    iraf.ccdproc.setParam('ccdtype', '')
    iraf.ccdproc.setParam('max_cac', 0)
    iraf.ccdproc.setParam('noproc', 'no')
    iraf.ccdproc.setParam('fixpix', 'yes')
    iraf.ccdproc.setParam('overscan', 'yes')
    iraf.ccdproc.setParam('trim', 'yes')
    iraf.ccdproc.setParam('zerocor', 'yes')
    iraf.ccdproc.setParam('darkcor', 'no')
    iraf.ccdproc.setParam('flatcor', 'no')
    iraf.ccdproc.setParam('illumco', 'no')
    iraf.ccdproc.setParam('fringec', 'no')
    iraf.ccdproc.setParam('readcor', 'no')
    iraf.ccdproc.setParam('scancor', 'no')
    iraf.ccdproc.setParam('readaxi', 'line')    
    iraf.ccdproc.setParam('fixfile', 'bad_pixels.dat')  
    iraf.ccdproc.setParam('biassec', '[2053:2080,2:2047]')   #Ivan's prescription
    iraf.ccdproc.setParam('trimsec', '[10:2040,2:2047]')     #Ivan's prescription 
    iraf.ccdproc.setParam('zero', 'zero')
    iraf.ccdproc.setParam('dark', '')
    iraf.ccdproc.setParam('flat', '')
    iraf.ccdproc.setParam('illum', '')
    iraf.ccdproc.setParam('fringe', '')
    iraf.ccdproc.setParam('minrepl', 1.)
    iraf.ccdproc.setParam('scantyp', 'shortscan')
    iraf.ccdproc.setParam('nscan', 1)
    iraf.ccdproc.setParam('interac', 'no')   #change to yes if you need to check the fit
    iraf.ccdproc.setParam('functio', 'legendre')
    iraf.ccdproc.setParam('order', 5)         # this order works very well
    iraf.ccdproc.setParam('sample', '*')
    iraf.ccdproc.setParam('naverag', 1)
    iraf.ccdproc.setParam('niterat', 5)
    iraf.ccdproc.setParam('low_rej', 3.)
    iraf.ccdproc.setParam('high_rej', 3.)
    iraf.ccdproc.setParam('grow', 0.)
    iraf.ccdproc.setParam('mode', 'q1')
    iraf.ccdproc(image=image)

#to create combined flat
def edit_flatcombine(input):
    iraf.flatcombine.setParam('output', 'flat.fits')
    iraf.flatcombine.setParam('combine', 'average')
    iraf.flatcombine.setParam('reject', 'avsigclip')
    iraf.flatcombine.setParam('ccdtype', '')
    iraf.flatcombine.setParam('process', 'yes')
    iraf.flatcombine.setParam('subsets', 'yes')
    iraf.flatcombine.setParam('delete', 'no')
    iraf.flatcombine.setParam('clobber', 'no')
    iraf.flatcombine.setParam('scale', 'mode')
    iraf.flatcombine.setParam('statsec', '')
    iraf.flatcombine.setParam('nlow', 0)
    iraf.flatcombine.setParam('nhigh', 1)
    iraf.flatcombine.setParam('nkeep', 1)
    iraf.flatcombine.setParam('mclip', 'yes')
    iraf.flatcombine.setParam('lsigma', 3.)
    iraf.flatcombine.setParam('hsigma', 3.)
    iraf.flatcombine.setParam('rdnoise', readnoise)
    iraf.flatcombine.setParam('gain', gain)
    iraf.flatcombine.setParam('snoise', 0.)
    iraf.flatcombine.setParam('pclip', -0.5)
    iraf.flatcombine.setParam('blank', 0.)
    iraf.flatcombine.setParam('mode', 'q1')
    iraf.flatcombine(input=input)

#to find apertures
def edit_apall(inp):
    iraf.apall.setParam('input', inp)
    iraf.apall.setParam('interactive', 'yes')
    iraf.apall.setParam('find', 'yes')
    iraf.apall.setParam('recenter', 'yes')
    iraf.apall.setParam('resize', 'yes')
    iraf.apall.setParam('edit', 'yes')
    iraf.apall.setParam('trace', 'yes')
    iraf.apall.setParam('fittrace', 'yes')
    iraf.apall.setParam('extract', 'yes')
    iraf.apall.setParam('extras', 'no')
    iraf.apall.setParam('review', 'no')
    iraf.apall.setParam('line', 1070)
    iraf.apall.setParam('nsum', 10)
    iraf.apall.setParam('width', 15.)
    iraf.apall.setParam('radius', 15.)
    iraf.apall.setParam('thresho', 0.)
    iraf.apall.setParam('minsep', 17)
    iraf.apall.setParam('maxsep', 45)
    #iraf.apall.setParam('nfind', 5.)          #60 in the tutorial reduction
    iraf.apall.setParam('shift', 'no')
    iraf.apall.setParam('avglimi', 'yes')
    iraf.apall.setParam('t_order', 3)
    iraf.apall.setParam('t_niterate', 5)
    iraf.apall.setParam('readnoi', readnoise)
    iraf.apall.setParam('gain', gain)
    iraf.apall(input=inp)

#to create combined flat
def edit_apflatten(input, interac):
    iraf.apflatten.setParam('output', 'flat_nrm.fits')
    iraf.apflatten.setParam('apertur', '')
    iraf.apflatten.setParam('referen', 'trz.fits')
    iraf.apflatten.setParam('interac', interac)
    iraf.apflatten.setParam('find', 'no')
    iraf.apflatten.setParam('recente', 'no')
    iraf.apflatten.setParam('resize', 'no')
    iraf.apflatten.setParam('edit', 'yes')
    iraf.apflatten.setParam('trace', 'no')
    iraf.apflatten.setParam('fittrac', 'no')
    iraf.apflatten.setParam('flatten', 'yes')
    iraf.apflatten.setParam('fitspec', 'yes')
    iraf.apflatten.setParam('readnoi', readnoise)
    iraf.apflatten.setParam('gain', gain)
    iraf.apflatten.setParam('line', 'INDEF')
    iraf.apflatten.setParam('nsum', 10)
    iraf.apflatten(input=input)

#apply normalize flat
def norm_ccdproc(image):
    iraf.ccdproc.setParam('images', image)
#    iraf.ccdproc.setParam('output', image)    #leaving ouput empty will overwrite the image, be careful!!
    iraf.ccdproc.setParam('ccdtype', '')
    iraf.ccdproc.setParam('max_cac', 0)
    iraf.ccdproc.setParam('noproc', 'no')
    iraf.ccdproc.setParam('fixpix', 'yes')
    iraf.ccdproc.setParam('overscan', 'yes')
    iraf.ccdproc.setParam('trim', 'yes')
    iraf.ccdproc.setParam('zerocor', 'no')
    iraf.ccdproc.setParam('darkcor', 'no')
    iraf.ccdproc.setParam('flatcor', 'yes')
    iraf.ccdproc.setParam('illumco', 'no')
    iraf.ccdproc.setParam('fringec', 'no')
    iraf.ccdproc.setParam('readcor', 'no')
    iraf.ccdproc.setParam('scancor', 'no')
    iraf.ccdproc.setParam('readaxi', 'line')    
    iraf.ccdproc.setParam('fixfile', 'bad_pixels.dat')  
    iraf.ccdproc.setParam('biassec', '[2053:2080,2:2047]')   #Ivan's prescription
    iraf.ccdproc.setParam('trimsec', '[10:2040,2:2047]')     #Ivan's prescription 
    iraf.ccdproc.setParam('zero', '')
    iraf.ccdproc.setParam('dark', '')
    iraf.ccdproc.setParam('flat', 'flat_nrm.fits')
    iraf.ccdproc.setParam('illum', '')
    iraf.ccdproc.setParam('fringe', '')
    iraf.ccdproc.setParam('minrepl', 1.)
    iraf.ccdproc.setParam('scantyp', 'shortscan')
    iraf.ccdproc.setParam('nscan', 1)
    iraf.ccdproc.setParam('interac', 'no')   #change to yes if you need to check the fit
    iraf.ccdproc.setParam('functio', 'legendre')
    iraf.ccdproc.setParam('order', 5)         # this order works very well
    iraf.ccdproc.setParam('sample', '*')
    iraf.ccdproc.setParam('naverag', 1)
    iraf.ccdproc.setParam('niterat', 5)
    iraf.ccdproc.setParam('low_rej', 3.)
    iraf.ccdproc.setParam('high_rej', 3.)
    iraf.ccdproc.setParam('grow', 0.)
    iraf.ccdproc.setParam('mode', 'q1')
    iraf.ccdproc(image=image)



################## object apertures #######################
#to trace object apertures
def trace_objectap_apall(inp, interac):
    iraf.apall.setParam('input', inp)
    iraf.apall.setParam('output', '')
    iraf.apall.setParam('apertur', '')
    iraf.apall.setParam('format', 'echelle')
    iraf.apall.setParam('referen', 'trz.fits')
    iraf.apall.setParam('profile', '')
    iraf.apall.setParam('interac', interac)
    iraf.apall.setParam('find', 'no')
    iraf.apall.setParam('recente', 'yes')
    iraf.apall.setParam('resize', 'yes')
    iraf.apall.setParam('edit', 'yes')
    iraf.apall.setParam('trace', 'yes')
    iraf.apall.setParam('fittrac', 'yes')   
    iraf.apall.setParam('extract', 'no')
    iraf.apall.setParam('extras', 'no')
    iraf.apall.setParam('t_order', 3)
    iraf.apall.setParam('t_niter', 5)
    iraf.apall.setParam('readnoi', readnoise)
    iraf.apall.setParam('gain', gain)
    iraf.apall(input=inp)

#to remove scattered light
def edit_apscatter(inp, interac):
    iraf.apscatter.setParam('input', inp)
    iraf.apscatter.setParam('output', '')
    iraf.apscatter.setParam('apertur', '')
    iraf.apscatter.setParam('scatter', '')
    iraf.apscatter.setParam('referen', '')
    iraf.apscatter.setParam('interac', interac)
    iraf.apscatter.setParam('find', 'no')
    iraf.apscatter.setParam('recente', 'no')
    iraf.apscatter.setParam('resize', 'no')
    iraf.apscatter.setParam('edit', 'no')
    iraf.apscatter.setParam('trace', 'no')
    iraf.apscatter.setParam('fittrac', 'no')   
    iraf.apscatter.setParam('subtrac', 'yes')
    iraf.apscatter.setParam('smooth', 'yes')
    iraf.apscatter.setParam('fitscat', 'yes')
    iraf.apscatter.setParam('fitsmoo', 'yes')
    iraf.apscatter.setParam('line', 'INDEF')
    iraf.apscatter.setParam('nsum', 10)
    iraf.apscatter.setParam('buffer', 1.0)
    iraf.apscatter.setParam('apscat1', '')
    iraf.apscatter.setParam('apscat2', '')
    iraf.apscatter.setParam('mode', 'q1')
    iraf.apscatter(input=inp)

#to re-trace and widen apertures
def retrace_objectap_apall(inp, interac):
    iraf.apall.setParam('input', inp)
    iraf.apall.setParam('interactive', interac)
    iraf.apall.setParam('find', 'yes')
    iraf.apall.setParam('recenter', 'yes')
    iraf.apall.setParam('resize', 'yes')
    iraf.apall.setParam('edit', 'yes')
    iraf.apall.setParam('trace', 'yes')
    iraf.apall.setParam('fittrace', 'yes')
    iraf.apall.setParam('extract', 'no')
    iraf.apall.setParam('extras', 'no')
    iraf.apall.setParam('review', 'no')
    iraf.apall.setParam('line', 1070)
    iraf.apall.setParam('nsum', 10)
    iraf.apall.setParam('width', 15.)
    iraf.apall.setParam('radius', 15.)
    iraf.apall.setParam('thresho', 0.)
    iraf.apall.setParam('minsep', 17)
    iraf.apall.setParam('maxsep', 45)
    #iraf.apall.setParam('nfind', 5.)          #60 in the tutorial reduction
    iraf.apall.setParam('shift', 'no')
    iraf.apall.setParam('avglimi', 'yes')
    iraf.apall.setParam('t_order', 3)
    iraf.apall.setParam('t_niterate', 5)
    iraf.apall.setParam('ylevel', 0.02)
    iraf.apall.setParam('readnoi', readnoise)
    iraf.apall.setParam('gain', gain)
    iraf.apall(input=inp)



################## wv calibration #######################
#to trace apertures of thar
def ap_thar_apall(inp, interac, clean):
    iraf.apall.setParam('input', inp)
    iraf.apall.setParam('output', '')
    iraf.apall.setParam('apertur', '')
    iraf.apall.setParam('format', 'echelle')
    iraf.apall.setParam('referen', 'trz.fits')
    iraf.apall.setParam('profile', '')
    iraf.apall.setParam('interac', interac)
    iraf.apall.setParam('find', 'no')
    iraf.apall.setParam('recente', 'no')
    iraf.apall.setParam('resize', 'no')
    iraf.apall.setParam('edit', 'no')
    iraf.apall.setParam('trace', 'no')
    iraf.apall.setParam('fittrac', 'no')   
    iraf.apall.setParam('extract', 'yes')
    iraf.apall.setParam('extras', 'no')
    iraf.apall.setParam('review', 'no')
    iraf.apall.setParam('line', 'INDEF')
    iraf.apall.setParam('nsum', 10)
    iraf.apall.setParam('width', 5.)
    iraf.apall.setParam('radius', 10.)
    iraf.apall.setParam('thresho', 0.)
    iraf.apall.setParam('minsep', 5)
    iraf.apall.setParam('maxsep', 100000.)
    #iraf.apall.setParam('nfind', 5.)          #60 in the tutorial reduction
    iraf.apall.setParam('shift', 'yes')
    iraf.apall.setParam('avglimi', 'no')
    iraf.apall.setParam('t_order', 2)
    iraf.apall.setParam('t_niterate', 0)
    iraf.apall.setParam('ylevel', 0.1)
    iraf.apall.setParam('clean', clean)
    iraf.apall.setParam('readnoi', readnoise)
    iraf.apall.setParam('gain', gain)
    iraf.apall(input=inp)

#to extract object spectra
def edit_doecslit(inp):
    iraf.doecslit.setParam('objects', inp)
    iraf.doecslit.setParam('apref', 'trz.fits')
    #iraf.doecslit.setParam('arcs', '@list_thar.txt')
    iraf.doecslit.setParam('arcs', 'thar.fits')
    iraf.doecslit(objects=inp)





################## useful python modules #######################
#check if file exist and then remove
def rm_file(input):
    for i in input:
        if os.path.exists(i):
            os.remove(i)

#remove and create not empy folder
def rem_cre_Folder(directory):
  try:
    if os.path.exists(directory):
      shutil.rmtree(directory)
    if not os.path.exists(directory):
      os.makedirs(directory)
  except OSError:
    print ('Error: Removing and creating directory. ' +  directory)

#copying files to path
def copy_files(input, directory):
    for i in input:
        shutil.copy(i, directory)

#moving files to path
def move_files(input, directory):
    for i in input:
        shutil.move(i, directory)


#Interquartile range (IQR)
def interquartile_range_filtering(data, k=1.5):
    """
    Interquartile range (IQR) is used to find outliers in data. By default, outliers
    are observations that fall below Quartile1 - k*(IQR) or above Quartile3 + k*(IQR).

    * k = 1.5 represents +/-2.698 * sigma (or standard dev) of a gaussian\
    distribution, which includes the 99.3% of the data.
    """
    # First and third quartile (the second is the median)
    q1 = np.percentile(data, 25) # 25% of the data (left to right)
    q3 = np.percentile(data, 75) # 75%
    # Interquartile range
    iqr = q3 - q1
    sfilter = np.logical_and(data > q1 - k * iqr, data < q3 + k * iqr)
    return data[sfilter], sfilter

#sigma clipping
def sigma_clipping(data, sig=3, meanfunc=np.mean):
    """
    Identify outliers considering the mean (if meanfunc=np.mean) or median (if meanfunc=np.median) value and 3 sigma (3*stdev),
    iterating until convergence.
    """
    last_total = len(data)

    # First iteration
    stdev = np.std(data)
    diff = data - meanfunc(data)
    sfilter = np.abs(diff) < sig*stdev
    current_total = len(data[sfilter])
    # Continue iterating until convergence (no more points are removed)
    while last_total > current_total:
        #print current_total, stdev
        last_total = current_total

        stdev = np.std(data[sfilter])
        diff = data - meanfunc(data[sfilter])
        sfilter = np.abs(diff) < sig*stdev

        current_total = len(data[sfilter])

    return data[sfilter], sfilter

#estimating snr
def estimate_snr(flux, num_points):
    """
    Estimate the Signal-to-Noise ratio for a given spectrum calculating the
    signal over standard deviation in blocks of N points and returning the average.
    """
    # Avoid negative values and outliers
    flux = flux[flux > 0.0]
    #flux, f = sigma_clipping(flux, sig=8, meanfunc=np.median)
    flux, f = interquartile_range_filtering(flux, k=1.5)
    if num_points == 1:
        snr = np.mean(flux) / np.std(flux)
    else:
        snr = []
        total_num_blocks = len(flux)-num_points
        for i in np.arange(total_num_blocks):
            values = flux[i:i+num_points]
            stdev = np.std(values)
            if stdev != 0:
                snr.append(np.mean(values) / stdev)
        snr = np.asarray(snr)
    #snr, s = sigma_clipping(snr, sig=8, meanfunc=np.median)
    snr, s = interquartile_range_filtering(snr, k=1.5)
    estimated_snr = np.mean(snr)
    return estimated_snr

# =====================================================================================

def DER_SNR(flux):
   
# =====================================================================================
   """
   DESCRIPTION This function computes the signal to noise ratio DER_SNR following the
               definition set forth by the Spectral Container Working Group of ST-ECF,
	       MAST and CADC. 

               signal = median(flux)      
               noise  = 1.482602 / sqrt(6) median(abs(2 flux_i - flux_i-2 - flux_i+2))
	       snr    = signal / noise
               values with padded zeros are skipped

   USAGE       snr = DER_SNR(flux)
   PARAMETERS  none
   INPUT       flux (the computation is unit independent)
   OUTPUT      the estimated signal-to-noise ratio [dimensionless]
   USES        numpy      
   NOTES       The DER_SNR algorithm is an unbiased estimator describing the spectrum 
	       as a whole as long as
               * the noise is uncorrelated in wavelength bins spaced two pixels apart
               * the noise is Normal distributed
               * for large wavelength regions, the signal over the scale of 5 or
	         more pixels can be approximated by a straight line
 
               For most spectra, these conditions are met.

   REFERENCES  * ST-ECF Newsletter, Issue #42:
               www.spacetelescope.org/about/further_information/newsletters/html/newsletter_42.html
               * Software:
	       www.stecf.org/software/ASTROsoft/DER_SNR/
   AUTHOR      Felix Stoehr, ST-ECF
               24.05.2007, fst, initial import
               01.01.2007, fst, added more help text
               28.04.2010, fst, return value is a float now instead of a numpy.float64
   """
   from numpy import array, where, median, abs 

   flux = array(flux)

   # Values that are exactly zero (padded) are skipped
   flux = array(flux[where(flux != 0.0)])
   n    = len(flux)      

   # For spectra shorter than this, no value can be returned
   if (n>4):
      signal = median(flux)

      noise  = 0.6052697 * median(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))

      return float(signal / noise)  

   else:

      return 0.0

# end DER_SNR -------------------------------------------------------------------------















