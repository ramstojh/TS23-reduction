from pyraf import iraf
from glob import glob
from astropy.io import fits as pyfits
import shutil
import numpy as np
import ts23_modules as ts23
import os

#path the pipeline
path = os.getcwd()

#horizontal dispersion in TS23
iraf.echelle.setParam('dispaxi', 1)  

################################################################################
#creating file to store
reduced_spectra = 'reduced_spectra'
ts23.rem_cre_Folder(path+'/'+reduced_spectra)

print ('                ')
print ('########  (1) Creating a reference Th-Ar exposure  ########')
#shutil.copy('backup/trz/trz.fits', path) 
ts23.rm_file(glob(path+'/database/apthar'))
list_thar = ['thar.fits', 'thar.ec.fits']
ts23.rm_file(list_thar)
ts23.copy_files(glob(path+'/backup/thar/*.fits'), path)

lista = glob('*.fits')
for i in lista:
    hdulist = pyfits.open(i)
    #print (hdulist.info())         #all infos have only the primary information
    hdu = hdulist[0]
    #print(hdu.header)
    #print(hdu.header['OBJECT'])    #identify objects (flat, bias, star)
    if hdu.header['OBJECT'] == 'Th-Ar Blue Pass':
        iraf.imcopy.setParam('output', 'thar.fits')
        iraf.imcopy(i)

print ('                ')
print ('########  (2) Extracting calibration lamp  ########')
ts23.ap_thar_apall('thar.fits', 'yes', 'yes')

print ('                ')
resp = raw_input('Would you like to identify lines to make wavelenght calibration? (write yes/no):')

ts23.rm_file(glob(path+'/database/ecthar.ec'))
if resp == 'yes':

    iraf.ecidentify('thar.ec.fits')
else:
    print ('Be sure that the ecthar.ec is into the **data_input** file')
    ts23.copy_files(glob(path+'/data_input/ecthar.ec'), path+'/database')
    iraf.eciden('thar.ec.fits')

iraf.imdel('@list_thar.txt')

#write here the loop for one spectra
print ('                ')
print ('########  (3) Extracting objects  ########')
ts23.rm_file(glob('*ec.fits'))
ts23.edit_doecslit('@list_objects.txt')

print ('                ')
print ('########  (4) The extracted objects are into the **final_spectra** file  ########')
ts23.move_files(glob('*ec.fits'), path+'/reduced_spectra')
ts23.rm_file(glob(path+'/reduced_spectra/thar.ec.fits'))

print ('                ')
print ('                ')
print ('########################################')
print ('########  Well, that was easy!  ########')
print ('########################################')
print ('                ')
print ('                ')




