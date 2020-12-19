from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from glob import glob
from PyPDF2 import PdfFileReader, PdfFileWriter
from PyAstronomy import pyasl
from pyraf import iraf
import os
import shutil
from astropy.io import fits as pyfits
import numpy as np
from tqdm import trange, tqdm
import ts23_modules as ts23

#########################################################################################################
#modules
def plotting():
  #  plt.rcParams['figure.figsize'] = (14, 8)
    plt.rcParams['font.size'] = 17                                     #set the size of the font numbers       
    plt.rcParams['font.family'] = 'fantasy'                            #choose the style of the numbers
    plt.rcParams['axes.labelsize'] = 20                                #set the size of the word axes
    #plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    #plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']        #set the size of the label x
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']        #set the size of the label y
    plt.rcParams['xtick.major.size'] = 7                               #set the length of the x tick
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['xtick.major.width'] = 2                              #set the width of the x tick 
    #plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 7                               #set the length of the y tick
    plt.rcParams['ytick.minor.size'] = 4
    plt.rcParams['ytick.major.width'] = 2                              #set the width of the y tick
    #plt.rcParams['ytick.minor.width'] = 1.2
    #plt.rcParams['legend.frameon'] = True
    #plt.rcParams['legend.loc'] = 'upper left'
    plt.rcParams['axes.linewidth'] = 3

#merging pdfs
def merge_pdfs(file_pdf, output):
  pdf_writer = PdfFileWriter()
  
  for file_pdf in file_pdf:
    pdf_reader = PdfFileReader(file_pdf)
    for page in range(pdf_reader.getNumPages()):
      # Add each page to the writer object
      pdf_writer.addPage(pdf_reader.getPage(page))

  # Write out the merged PDF
  with open(output, 'wb') as out:
    pdf_writer.write(out)

################################################################################
#file names
path            = os.getcwd()    #path the pipeline 
reduced_spectra = 'reduced_spectra' 

#removing unnecessary data
for i in glob('apt*'):
    if os.path.isfile(path+'/'+i):
        os.remove(path+'/'+i)

file_remo = glob('ec.fits*')
for i in file_remo:
    if os.path.isfile(path+'/'+i):
        os.remove(path+'/'+i)

for i in glob('.png*'):
    if os.path.isfile(path+'/'+i):
        os.remove(path+'/'+i)

#########################################################################################################
#def making plots
def make_plot(input):
    for i in input:
        hdulist   = pyfits.open(i)
        hdu         = hdulist[0]
        object_id   = hdu.header['OBJECT']

        input_fits  = [i + '[*,' + str(n) +']' for n in range(1, 58, 1)]   #58
        output_fits = ['apt_' + str('{:02d}'.format(n)) + '_' + object_id  + '_' + i for n in range(1, 58, 1)]  #58
        
        for j in range(len(input_fits)):
            iraf.scopy(input_fits[j], output_fits[j])

    #making plots
    snr_d = []
    snr_e = []
    wv    = []
    cutted_spectra = glob('*ap*.fits')
    for i in trange(len(cutted_spectra), desc='making plots'):
        wvls, flxs   = pyasl.read1dFitsSpec(cutted_spectra[i])

        plotting()
        
        fig = plt.figure(figsize=(23,11))
        gs  = gridspec.GridSpec(1, 1)#, height_ratios=[5, 1.5])
        ax = plt.subplot(gs[0])
        ax.plot(wvls, flxs, linewidth = 2., color='blue', label = cutted_spectra[i])

        plt.legend(loc='best', numpoints=1, prop={'size':15}, shadow=True)
        plt.ylabel(r'$\rm{Flux}$', color='black', size = 30)
        plt.xlabel(r'$\rm{Wavelength}\ \rm{(\AA)}}$', color='black', size = 30)
        plt.grid()
        plt.savefig(cutted_spectra[i]+'.pdf')
        #plt.show()

        #estimating snr in each aperture
        snr_der  = ts23.DER_SNR(flxs)
        snr_est  = ts23.estimate_snr(flxs, 10)
        snr_d.append(snr_der) 
        snr_e.append(snr_est)
        wvls_at = np.mean(wvls)
        wv.append(wvls_at)
        #print (wvls_at)

    fig = plt.figure(figsize=(18,14))
    plt.scatter(wv, snr_d, c='k', label='derivated SNR', s=90)
    plt.scatter(wv, snr_e, c='r', label='estimated SNR', s=90)
    plt.legend(loc='best', numpoints=1, prop={'size':15}, shadow=True)
    plt.ylabel(r'$\rm{<SNR>}$', color='black', size = 35)
    plt.xlabel(r'$\rm{Wavelength}\ \rm{(\AA)}}$', color='black', size = 35)
    plt.grid()
    #plt.show()
    plt.savefig('SNR_'+input[0]+'.png')

    #Merging and organizing files
    pdf_files = np.sort(glob('*.pdf'))
    merge_pdfs(pdf_files, input[0]+'.pdf')
    
    #moving files
    shutil.move(input[0]+'.pdf', path+'/'+reduced_spectra)
    shutil.move('SNR_'+input[0]+'.png', path+'/'+reduced_spectra)

    if os.path.isfile(input[0]):
        os.remove(input[0])

    file_remo = glob('apt*')
    for i in file_remo:
        if os.path.isfile(path+'/'+i):
            os.remove(path+'/'+i)

##########################################################################################################
#the code begins here
list_reduc = glob(path+'/reduced_spectra/*ec.fits')
for i in list_reduc:
    shutil.copy(i, path)

#making plots
spec = glob('*ec.fits')
for i in spec:
        make_plot([i])

