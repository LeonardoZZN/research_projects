# clm:  7/02/2019
#       - Add fitted continuum level
#       - Add reduced chi2
# clm:  6/28/2019
#       - Fit y = [m(wave - w0) + y0] + a * exp(-[wave-waveB(1+z)]^2/2sig^2) + dr * a * exp(-[wave-waveR(1+z)]^2/2sig^2)
#
# Parameters:
#    wb:  rest wavelength of short wavelength member of doublet
#    wr:  rest wavelength of long wavelength member of doublet
#    dr:  doublet ratio
#    w0:  wavelength where continuum level is y0
# Fitted Quantities:
#    z:      redshift
#    sigma:  standard deviation of Gaussian line profile
#    a:      amplitude of shorter wavelength member of doublet
#    m,y0:   continuum slope and normalization at wavelength w0

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.io import ascii
from astropy.io import fits

### Parameters for the doubelt ###
# Vacuum setting
# Blue line
wb = 3727.09
# Red line
wr = 3729.88

# Doublet ratio = F(red) / F(blue)
# LDL
dr = 1.5
# HDL
#dr = 0.3

#func = lambda x, z, sigma, a:  a * np.exp(-(x - wb * (1.+ z))**2 / (2. * sigma **2)) + dr * a * np.exp(-(x - wr * (1.+ z))**2 / (2. * sigma **2))

func = lambda x, z, sigma, a, m, y0:  a * np.exp(-(x - wb * (1.+ z))**2 / (2. * sigma **2)) + dr * a * np.exp(-(x - wr * (1.+ z))**2 / (2. * sigma **2)) + m * (x - w0) + y0




#### Input Spectrum and Initial Guess ####
#input = 'test.txt'
#
#input = 'test_apcal_m5-width8.txt'
#guess = np.array([1.4727, 4.40, 200, 0, 0])

#input = 'test_apcal_m4-width6.txt'
#guess = np.array([1.4833, 4.40, 200, 0, 0])

#input = 'test_apcal_s2-width5.txt'
#input = 'test_apcal_s2-width5.txt'


input = 'OIIB146m-5_w15_fcal.fits'
#guess = np.array([1.4767, 4.40, 200, 0, 0])  # redshift, sigma (A), amplitude, continuum slope, continuum normalization
guess = np.array([1.473, 1.00, 6e-18, 0, 0])  # redshift, sigma (A), amplitude, continuum slope, continuum normalization






# Continuum through point (w0,y0) with slope m
w0 = 9210.




## Read the spectrum ##
# Observed Wavelength (A)
# Flux_lambda (arbitrary)
#list = ascii.read('apcal_m5-width8.txt')
"""
list = ascii.read(input)
wave = np.array(list['wave(A)'])
flux = np.array(list['flux'])
err = np.array(list['err'])
"""
hdulist = fits.open(input)
wave = np.array(hdulist[0].data)
flux = np.array(hdulist[1].data)
err = np.array(hdulist[2].data)

nspec = wave.size


if True:  # sanity check the input
    fig1 = plt.figure()
    plt.plot(wave,flux, label='flux')
    plt.plot(wave, err, label='err')
    plt.legend()
    fig1.show()
    print('Zoom in around 9200 to 9250 Angstrom to see the [OII] doublet.')



if True:  # sanity check the initial guess
    yguess = func(wave,guess[0],guess[1],guess[2],guess[3],guess[4])    
    fig1 = plt.figure()
    plt.plot(wave,flux, label='flux')
    plt.plot(wave, err, label='err')
    plt.plot(wave, yguess, label='guess')
    plt.legend()
    fig1.show()
    print('Adjust the guess to be reasonable.')






# Fit with LM
p, pcov = curve_fit(func, wave, flux, guess, err)
perr = np.sqrt(np.diag(pcov))

print ('Fit Completed')

# Calculate 'Goodness of Fit'
model = func(wave,p[0],p[1],p[2],p[3],p[4])
term = ((flux-model) / err)**2  
dof = wave.size - guess.size
chi2 = np.sum(term) / dof



# Plot the results
yguess = func(wave,guess[0],guess[1],guess[2],guess[3],guess[4])
yfit = func(wave,p[0],p[1],p[2],p[3],p[4])

yblue = lambda x, z, sigma, a, m, y0:       a * np.exp(-(x - wb * (1.+ z))**2 / (2. * sigma **2)) + m * (x - w0) + y0
yred  = lambda x, z, sigma, a, m, y0:  dr * a * np.exp(-(x - wr * (1.+ z))**2 / (2. * sigma **2)) + m * (x - w0) + y0

fig = plt.figure()
plt.plot(wave,flux,'k')
plt.plot(wave,yfit,'r')
plt.plot(wave,yblue(wave,p[0],p[1],p[2],p[3],p[4]), 'b:')
plt.plot(wave,yred(wave,p[0],p[1],p[2],p[3],p[4]), 'r:')
#
#plt.plot(wave,yguess,'k:')

dy = np.max(flux) - np.min(flux)
dx = np.max(wave) - np.min(wave)
epsx = 0.01 * dx
epsy = 0.1 * dy

# formatting the fitted values for printing
zlab = (f'z =  {p[0]:.5f} ( {perr[0]:.5f} )')
siglab = (f'sigma =  {p[1]:.5f} ( {perr[1]:.5f} )')
alab = (f'A =  {p[2]:.5f} ( {perr[2]:.5f} )')
contlab = (f'm =  {p[3]:.5f} ( {perr[3]:.5f} ); y0 =  {p[4]:.5f} ( {perr[4]:.5f} )')
goodlab = (f'chi2 =  {chi2:.5f} ')

plt.text(np.min(wave) + epsx, np.max(flux), zlab, fontsize = 12)
plt.text(np.min(wave) + epsx, np.max(flux) -  epsy, siglab, fontsize = 12)
plt.text(np.min(wave) + epsx, np.max(flux) - 2 * epsy, alab, fontsize = 12)
plt.text(np.min(wave) + epsx, np.max(flux) - 3 * epsy, contlab, fontsize = 12)
plt.text(np.min(wave) + epsx, np.max(flux) - 4 * epsy, goodlab, fontsize = 12)

plt.title(input + ' doublet ratio = ' +  str(dr))

plt.xlabel("Wavelength (A)")
plt.ylabel("Flux")
fig.show()

#temp = input[:-4]
#name = temp[11:]
name = input[:-4]
outfile = "plot_" + name + "_dr_" + str(dr) + ".ps"
plt.savefig(outfile,facecolor="w",bbox_inches="tight",pad_inches=0.1,format='eps') 
