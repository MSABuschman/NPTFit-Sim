import ps_mc

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

#User given index, flux breaks, and log-normilization
#n = np.array([16.99, -0.90])
n = np.array([17.79, -0.99])
#S = np.array([24.24])
S = np.array([14.29])
#A = -2.84
A = -2.34
#Path to relevant template
temp = "../fermi_data/template_gce.npy"
exp = "../fermi_data/fermidata_exposure.npy"
#Name of output file
name = "example"

#S is in terms of counts, converted to flux with mean exposure.
EXP = np.load("../fermi_data/fermidata_exposure.npy")
mean_exp = np.mean(EXP)
F = S / mean_exp

#Also account for distribution from counts to flux for norm. factor
cor_term = np.log10(mean_exp)
A = A + cor_term

# Define parameters that specify the Fermi-LAT PSF at 2 GeV
fcore = 0.748988248179
score = 0.428653790656
gcore = 7.82363229341
stail = 0.715962650769
gtail = 3.61883748683
spe = 0.00456544262478

# Define the full PSF in terms of two King functions
def king_fn(x, sigma, gamma):
    return 1./(2.*np.pi*sigma**2.)*(1.-1./gamma)*(1.+(x**2./(2.*gamma*sigma**2.)))**(-gamma)

def Fermi_PSF(r):
    return fcore*king_fn(r/spe,score,gcore) + (1-fcore)*king_fn(r/spe,stail,gtail)

#Lambda function to pass user defined PSF
psf_r = lambda r: r * Fermi_PSF(r)

#Run the simulation
ps_mc.run(n,F,A,temp,exp,psf_r,name)
