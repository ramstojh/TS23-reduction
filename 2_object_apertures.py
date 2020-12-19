from pyraf import iraf
from glob import glob
import shutil
import ts23_modules as ts23
import os

#path the pipeline
path = os.getcwd()

#horizontal dispersion in TS23
iraf.echelle.setParam('dispaxi', 1)  

#removing unnecessary files
ts23.rm_file(glob('*.fits'))
ts23.copy_files(glob(path+'/backup/objects/*.fits'), path)

#Tracing object apertures
print ('                ')
print ('                ')
print ('########  (1) Tracing the apertures of all spectra  ########')
ts23.trace_objectap_apall('@list_objects.txt', 'yes')   #again, we use o=3 and nit=5 for all

#removing scattered light
print ('                ')
print ('########  (2) Removing scattered light  ########')
print ('Answer yes to everything.')
ts23.edit_apscatter('@list_objects.txt', 'yes')    #o=12, nit=10 for the row, and o=6, nit=5 for the column  

print ('                ')
resp = raw_input('Are you happy with the removed scattered light? (write yes/no):')

while resp == 'no':
    iraf.imdel('@list_objects.txt')
    ts23.copy_files(glob(path+'/backup/objects/*.fits'), path)
    ts23.trace_objectap_apall('@list_objects.txt', 'yes')
    ts23.edit_apscatter('@list_objects.txt', 'yes')
    print ('                ')
    print ('                ')
    resp = raw_input('Are you happy with the removed scatterd light? (write yes/no):')
    if resp == 'yes':
        break
if resp == 'yes':
    #retracing and widing apertures
    print ('########  (3) Re-tracing and widen apertures to obtain the highest SNR  ########')
    ts23.retrace_objectap_apall('@list_objects.txt', 'yes')

print ('                ')
print ('                ')
print ('------------------------')
print ('The step 2 has finished!!')
print ('------------------------')







