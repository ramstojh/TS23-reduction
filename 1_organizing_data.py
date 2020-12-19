import pandas as pd
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
#file names
data_input      = 'data_input'
data_to_reduce  = 'data_to_reduce'
database        = 'database'

#copy all spectra from data_to_reduce file to pipeline file
print ('                ')
print ('                ')
print ('########  (1) Organizing files  ########')
#removing unnecessary data
ts23.rm_file(glob('*.fits'))
ts23.rm_file(glob('*.txt'))
ts23.rm_file(glob('*.pdf'))

#creating files
ts23.rem_cre_Folder(path+'/'+database)
ts23.rem_cre_Folder(path+'/'+'backup')
ts23.rem_cre_Folder(path+'/'+'bias')
ts23.rem_cre_Folder(path+'/'+'flat')
ts23.rem_cre_Folder(path+'/'+'thar')
ts23.rem_cre_Folder(path+'/'+'trz')
ts23.rem_cre_Folder(path+'/'+'objects')

#copying spectra to path
ts23.copy_files(glob(path+'/'+data_to_reduce+'/'+'*.fits'), path)
lista = np.sort(glob('*.fits'))
np.savetxt('list_all.txt', lista, fmt="%s")

#overscan correction
print ('                ')
print ('########  (2) Making overscan correction  ########')
#be carefull it will rename all the *fits
for i in range(len(lista)):
    ts23.edit_ccdproc(lista[i], lista[i])

#creating empty list
bias    = []
flat    = []
thar    = []
objects = []

#organizing files
print ('                ')
print ('########  (3) Organizing ** .fits ** files  ########')
for i in lista:
    hdulist = pyfits.open(i)
    #print (hdulist.info())         #all infos have only the primary information
    hdu = hdulist[0]
    #print(hdu.header)
    #print(hdu.header['OBJECT'])    #identify objects (flat, bias, star)
    if hdu.header['OBJECT'] == 'zero':
        bias.append(i)
        np.savetxt('list_zero.txt', np.sort(bias), fmt="%s")
        #print('bias',i)
        #shutil.copy(i, 'bias/')
    elif hdu.header['OBJECT'] == 'bias':
        bias.append(i)
        np.savetxt('list_zero.txt', np.sort(bias), fmt="%s")
        #print('bias',i)
        #shutil.copy(i, 'bias/')
    elif hdu.header['OBJECT'] == 'flat':
        flat.append(i)
        np.savetxt('list_flat.txt', np.sort(flat), fmt="%s")
        #print('flat',i)
        #shutil.copy(i, 'flat/')
    elif hdu.header['OBJECT'] == 'Th-Ar Blue Pass':
        thar.append(i)
        np.savetxt('list_thar.txt', np.sort(thar), fmt="%s")
        #print('thar',i)
        #shutil.copy(i, 'thar/')
    else:
        objects.append(i)
        np.savetxt('list_objects.txt', np.sort(objects), fmt="%s")
        #print ('objects', i)
        #shutil.copy(i, 'objects/')

print ('                ')
print ('########  (4) Creating a combined zero image: zero.fits  ########')
ts23.edit_zerocombine('@list_zero.txt')       #zero combine

print ('                ')
print ('########  (5) Applying zero correction  ########')
ts23.edit_zerocorrection('@list_all.txt')     #zero correction

print ('                ')
print ('########  (6) Creating a combined flat: flat.fits  ########')
ts23.edit_flatcombine('@list_flat.txt')       #flat combine

print ('                ')
resp = raw_input('Do you want to define the apertures? (write yes/no):')

if resp == 'yes':
    #taking the image with the highest MEAN
    iraf.imstat('@list_objects.txt', Stdout="imstat_mean.txt")
    imst = pd.read_csv('imstat_mean.txt', comment='#', delim_whitespace=True, names = ['IMAGE', 'NPIX', 'MEAN', 'STDDEV', 'MIN', 'MAX'])
    idx = imst['MEAN'] == imst['MEAN'].max()
    image_max = imst['IMAGE'][idx].values.tolist()     #image with the highest MEAN
    iraf.imcopy.setParam('output', 'trz.fits')
    iraf.imcopy(image_max[0])
    ts23.edit_apall('trz.fits')        #Defining apertures
elif resp == 'no':
    shutil.copy(path+'/'+data_input+'/trz.fits', path)
    shutil.copy(path+'/'+data_input+'/aptrz', path+'/'+database)
    print ('You choosed NO, so do not forget to provide the apertures!!.')
    print ('The aptrz and trz.fits should be into the ** data_input ** file.')


#Normalize flat
print ('                ')
print ('########  (7) Normalizing the flat: flat_nrm  ########')
ts23.edit_apflatten('flat.fits', 'yes')        #we use o=1 and nit=5. Wide the apertures from the trace

#divide the object spectra by the flat_nrm
print ('                ')
print ('                ')
print ('########  (8) Applying the ** flat_nrm ** to the object spectra  ########')
ts23.norm_ccdproc('@list_objects.txt')

for i in lista:
    hdulist = pyfits.open(i)
    hdu = hdulist[0]
    #if hdu.header['OBJECT'] == 'zero':
    if hdu.header['OBJECT'] == 'bias':
        bias.append(i)
        shutil.move(i, 'bias/')
    elif hdu.header['OBJECT'] == 'flat':
        flat.append(i)
        shutil.move(i, 'flat/')
    elif hdu.header['OBJECT'] == 'Th-Ar Blue Pass':
        thar.append(i)
        shutil.move(i, 'thar/')
    else:
        objects.append(i)
        shutil.move(i, 'objects/')

#making a backup
files_flat = ['flat_nrm.fits', 'flat.fits']
ts23.move_files(files_flat, path+'/flat')
shutil.move('zero.fits', path+'/bias')
shutil.move('trz.fits', path+'/trz')
ts23.rm_file(glob('*.fits'))
files = ['bias', 'flat', 'objects', 'thar', 'trz']
ts23.move_files(files, path+'/backup')

print ('                ')
print ('--------------------------------------------------------------------------------')
print ('Check the ** .txt ** files to get information about flat, bias, lamps and stars.')
print ('--------------------------------------------------------------------------------')



