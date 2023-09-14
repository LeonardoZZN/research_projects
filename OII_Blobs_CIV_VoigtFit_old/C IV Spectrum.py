#!/usr/bin/env python
# coding: utf-8

import numpy as np
from astropy import constants as const
from astropy.io import fits
import matplotlib.pyplot as plt
import VoigtFit
import sys

target_name = 'OIIB146m-4'

#Input fits files
file_fcal = "../LRIS data/1D spectra/OIIB146m-4_w7_fcal.fits"
file_iraf = "../LRIS data/1D spectra/OIIB146m-4_w7_iraf.fits"
file_iraf_E = "../LRIS data/1D spectra/OIIB146m-4_w7_iraf_E.fits"

#Extract data from fits files
fits_fcal = fits.open(file_fcal)
fits_iraf = fits.open(file_iraf)
fits_iraf_E = fits.open(file_iraf_E)

wave = fits_fcal[0].data #Wavelength
flux = fits_iraf[0].data #Flux
error = fits_iraf_E[0].data #Flux uncertainty

#Add a mask on wavelengths from 3700A to 3900A
mask = np.where((wave>3700)&(wave<3900))


'''
# Plot the LRIS spectrum (in wavelength and velocity space)
fig1 = plt.figure(figsize=(6,8))
ax1  = fig1.add_subplot(211)
ax2  = fig1.add_subplot(212)

# fiting C IV: 1550.781, 1548.204 [wrest from Morton+03]
civ_wrest = np.array([1548.204, 1550.781])

# flux vs wavelength
ax1.step(wave[mask],flux[mask],color='k', where='mid')
ax1.step(wave[mask],error[mask],color='g', where='mid')
ax1.set_xlabel("Wavelength (A)")
ax1.set_ylabel("Flux")

# flux vs velocity
z_sys = 1.48260 # copied from z
vel   = (wave/(1.+z_sys) - civ_wrest[0])/civ_wrest[0] * 2.9979e5
ax2.step(vel[mask],flux[mask],color='k', where='mid')
ax2.step(vel[mask],error[mask],color='g', where='mid')
ax2.set_xlabel("velocity [km/s]")
ax2.set_ylabel("Flux")

fig1.show()
sys.exit()
'''

# -- The VoigtFit
#The spectral resolution of LRIS
R = 3820/10
res = const.c.to('km/s')/R
print(res.value)

#Set up the dataset
#z = 1.48260
#dataset = VoigtFit.DataSet(z)
#dataset.add_data(wave, flux, res.value, err=error, normalized=False)

#Load a dataset
dataset = VoigtFit.load_dataset('test_civ.dataset.hdf5')

#Define absorption lines
dataset.add_line('CIV_1550',10000)
dataset.add_line('CIV_1548',10000)

dataset.reset_components()

'''
FWHM_b = 600
FWHM_r_ = 500
b = FWHM_b/1.665
'''

# Set C IV velocity components
dataset.add_component_velocity('CIV', -4700.0, 75.0, 16.5, var_z=1, var_b=1, var_N=1)  
dataset.add_component_velocity('CIV', -460.0, 70.0, 15.0, var_z=1, var_b=1, var_N=1)  

# Set fitting order of the continuum
#   cheb_order = -1 : fit a linear function
#   cheb_order = 0  : fit a horizontal line
#   cheb_order = 1 means a straight line
#   cheb_order = 2 means a quadratic function
#   cheb_order = 3 means a cubic function
dataset.cheb_order = 1
dataset.prepare_dataset(norm=True, mask=True) # norm = True to fit continuum

## Save the normalized dataset when you have a satisfied Continuum
#  and masking, so you don't have to normalize and mask every time
#  i.e., you can load this normalized dataset in the future
#  [uncomment the next line if you want to save now]
#dataset.save('CIV_norm.dataset' % target_name)

popt, chi2 = dataset.fit()

dataset.plot_fit()

dataset.set_name('CIV')
dataset.save_fit_regions()

# save the best fit result
dataset.save_parameters(filename='test_civ.best_fit.txt')
# save the dataset (input spectra, initial guesses, best fit parameters,
# masketc.) as a hdf5 file
dataset.save(filename='test_civ.dataset.hdf5')

