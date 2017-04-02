# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:33:29 2017

@author: jph
"""

# Import Modules
import numpy as np
import matplotlib.pyplot as plt
#import math
#import lmfit
import camb
#from camb import model, initialpower

inFile = "COM_PowerSpect_CMB-TT-hiL-binned_R2.csv"
# Load Data (with bin values)
x_min, x_max, y, y_err = np.loadtxt(inFile, delimiter=',', 
                        usecols=(1, 2, 3, 4),
                         unpack=True, skiprows=3)

# plotted x value                          
x = x_min + (x_max - x_min)/2.

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#def function():
#This function sets up the cosmological parameters
pars.set_cosmology(H0=None, ombh2=0.02222, omch2=0.1197, omk=0, tau=0.078, cosmomc_theta = 0.010408)
pars.InitPower.set_params(ns=0.965, As=2.2E-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0)

# get results
results = camb.get_results(pars)

#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
#for name in powers: print(name)

#plot the total lensed CMB power spectra versus unlensed, and fractional difference
#All power spectra are l(l+1)C_l/2pi
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
#spectrum = []
#for name in powers : 
#    spectrum.append((7.4311)*(1E12)*powers[name][:,2])
#print(totCL.shape)

ls = np.arange(totCL.shape[0])

# unit conversion
spectrum = (7.4311)*(1E12)*totCL[:,0]


# Plot the data and camb models
plt.title(r'Figure 3: $TT$ Plot with Varying $\Omega_b h^2$')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$T_0^2 \ell (\ell + 1) C_\ell^{TT} / 2\pi$')
plt.errorbar(x, y, yerr=y_err, fmt='.', color='k', label=r'High $\ell$ Planck Data')
#plt.show()
#plt.plot(ls, spectrum[0], label='ult')
#plt.plot(ls, spectrum[1], label='uls')
#plt.plot(ls, spectrum[2], label='tens')
#plt.plot(ls, spectrum[3], label='lenspot')
#plt.plot(ls, spectrum[4], label='tot')
#plt.plot(ls, spectrum[5], label='lens')
plt.plot(ls[2:], spectrum[2:], label=r'Best Fit')

# repeat for 1.5*ombh2
hi_pars = camb.CAMBparams()
hi_pars.set_cosmology(H0=None, ombh2=0.03333, omch2=0.1197, omk=0, tau=0.078, cosmomc_theta = 0.010408)
hi_pars.InitPower.set_params(ns=0.965, As=2.2E-9)
hi_pars.set_for_lmax(2500, lens_potential_accuracy=0)
hi_results = camb.get_results(hi_pars)
hi_powers = hi_results.get_cmb_power_spectra(hi_pars)
hi_totCL=hi_powers['total']
hi_ls = np.arange(hi_totCL.shape[0])
hi_spectrum = (7.4311)*(1E12)*hi_totCL[:,0]
plt.plot(hi_ls[2:], hi_spectrum[2:], label=r'1.5$\Omega_b h^2$ Fit')

# repeat for 0.5*ombh2
lo_pars = camb.CAMBparams()
lo_pars.set_cosmology(H0=None, ombh2=0.01111, omch2=0.1197, omk=0, tau=0.078, cosmomc_theta = 0.010408)
lo_pars.InitPower.set_params(ns=0.965, As=2.2E-9)
lo_pars.set_for_lmax(2500, lens_potential_accuracy=0)
lo_results = camb.get_results(lo_pars)
lo_powers = lo_results.get_cmb_power_spectra(lo_pars)
lo_totCL=lo_powers['total']
lo_ls = np.arange(lo_totCL.shape[0])
lo_spectrum = (7.4311)*(1E12)*lo_totCL[:,0]
plt.plot(lo_ls[2:], lo_spectrum[2:], label=r'0.5$\Omega_b h^2$ Fit')

#plt.xscale('log')
#plt.show()

# lowL plot
inLowFile = "COM_PowerSpect_CMB-TT-loL-full_R2.csv"
# Load Data (with bin values)
lx, ly, y_up, y_down = np.loadtxt(inLowFile, delimiter=',', 
                        usecols=(0, 1, 2, 3),
                         unpack=True, skiprows=3)

plt.errorbar(lx, ly, yerr=[y_down, y_up], fmt='.', label=r'Low $\ell$ Planck Data')

plt.legend()
#lls = ls[2:len(x)+2]
lspectrum = spectrum[2:len(x)+2]
hi_lspectrum = hi_spectrum[2:len(x)+2]
lo_lspectrum = lo_spectrum[2:len(x)+2]
#plt.plot(lls, lspectrum, label=r'Low $\ell$ fit')
#plt.title(r'$TT$ at low $\ell$')
#plt.xlabel(r'$10^\ell$')
#plt.ylabel(r'$T_0^2 \ell (\ell + 1) C_\ell^{TT} / 2\pi$')
#plt.xscale('log')

# calculate chisquare
interpolated = []
hi_interpolated = []
lo_interpolated = []
for i in range(len(x) - 1) :
    interpolated.append(spectrum[int(x[i])] + 
        (spectrum[int(x[i+1])] - spectrum[int(x[i])])/(x_max[i] - x_min[i])/2)
    hi_interpolated.append(hi_spectrum[int(x[i])] + 
        (hi_spectrum[int(x[i+1])] - hi_spectrum[int(x[i])])/(x_max[i] - x_min[i])/2)
    lo_interpolated.append(lo_spectrum[int(x[i])] + 
        (lo_spectrum[int(x[i+1])] - lo_spectrum[int(x[i])])/(x_max[i] - x_min[i])/2)


# divide by two because we are increasing going up the interpolation slop by .5

chisq = 0
hi_chisq = 0
lo_chisq = 0
for i in range(len(interpolated)) : 
    chisq += (interpolated[i]-y[i])**2/y_err[i]**2
    hi_chisq += (hi_interpolated[i]-y[i])**2/y_err[i]**2
    lo_chisq += (lo_interpolated[i]-y[i])**2/y_err[i]**2
#chisq /= len(y) - 1
# for small l
lchisq = 0
hi_lchisq = 0
lo_lchisq = 0
for i in range(len(lx)) : 
    if(spectrum[int(lx[i])] > y[i]) :
        lchisq += (lspectrum[int(lx[i])] - ly[i])**2/y_up[i]**2
    else : 
        lchisq += (lspectrum[int(lx[i])] - ly[i])**2/y_down[i]**2
for i in range(len(lx)) : 
    if(hi_spectrum[int(lx[i])] > y[i]) :
        hi_lchisq += (hi_lspectrum[int(lx[i])] - ly[i])**2/y_up[i]**2
    else : 
        hi_lchisq += (hi_lspectrum[int(lx[i])] - ly[i])**2/y_down[i]**2
for i in range(len(lx)) : 
    if(lo_spectrum[int(lx[i])] > y[i]) :
        lo_lchisq += (lo_lspectrum[int(lx[i])] - ly[i])**2/y_up[i]**2
    else : 
        lo_lchisq += (lo_lspectrum[int(lx[i])] - ly[i])**2/y_down[i]**2
#lchisq /= len(ly) - 1
chisq = (chisq + lchisq)/(len(ly) + len(y) - 1)
hi_chisq = (hi_chisq + hi_lchisq)/(len(ly) + len(y) - 1)
lo_chisq = (lo_chisq + lo_lchisq)/(len(ly) + len(y) - 1)

# print results 
print(r'best reduced $\chi^2$ =', chisq)
print(r'1.5$\Omega_b h^2$ reduced $\chi^2$ =', hi_chisq)
print(r'0.5$\Omega_b h^2$ reduced $\chi^2$ =', lo_chisq)
print(results.get_derived_params())