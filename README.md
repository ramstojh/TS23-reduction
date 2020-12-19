# TS23-reduction

Useful semi-automathic python/IRAF scripts to reduce TS23 spectra. 

How to use?

1) Create a file and named it as data_input, an then put all your ts23 sectra (including bias, flats and only one calibration lamp) within this file. 
2) Run the 1_organizing_data.py scritp.py to check if all your data is ok.
3) Run the 2_object_apertures.py to choose the apertures manually and correct dispersion light.
4) Run 3_wv_calibration.py to perform wavelength calibration.
5) Run 4_plotting.py (it is optional) if you want to generate a pdf plot of the reduced spectra for each aperture.

## Importants Remarks
The scripts works only with TS23 spectra in its high resolution mode (R = 60 000).
