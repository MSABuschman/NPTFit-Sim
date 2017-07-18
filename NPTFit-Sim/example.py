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
temp = "../../fermi_data/template_gce.npy"
exp = "../../fermi_data/fermidata_exposure.npy"
#Name of output file
name = "example"

#S is in terms of counts, converted to flux with mean exposure.
EXP = np.load(exp)
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
psf_r = lambda r: np.sin(r) * Fermi_PSF(r)

#Run the simulation
ps_mc.run(n,F,A,temp,exp,psf_r,name)
sim = np.load("example.npy").astype(np.int32)

from NPTFit import nptfit
from NPTFit import create_mask as cm
from NPTFit import dnds_analysis
from NPTFit import psf_correction as pc

n = nptfit.NPTF(tag='GCE_Example')
n.load_data(sim,EXP)

gce = np.load("../../fermi_data/template_gce.npy")
# Remove exposure correction for PS
rescale = EXP / np.mean(EXP)
n.add_template(gce / rescale, 'gce_np', units='PS')

#Add the Non-poisson model
n.add_non_poiss_model('gce_np',
                      ['$A_\mathrm{gce}^\mathrm{ps}$','$n_1^\mathrm{gce}$','$n_2^\mathrm{gce}$','$S_b^{(1), \mathrm{gce}}$'],
                      [[-6,1],[2.05,30],[-2,1.95],[0.05,40]],
                      [True,False,False,False])

# Modify the relevant parameters in pc_inst and then make or load the PSF
pc_inst = pc.PSFCorrection(delay_compute=True)
pc_inst.psf_r_func = lambda r: Fermi_PSF(r)
pc_inst.sample_psf_max = 10.*spe*(score+stail)/2.
pc_inst.psf_samples = 10000
pc_inst.psf_tag = 'Fermi_PSF_2GeV'
pc_inst.make_or_load_psf_corr()

# Extract f_ary and df_rho_div_f_ary as usual
f_ary = pc_inst.f_ary
df_rho_div_f_ary = pc_inst.df_rho_div_f_ary

#pc_inst = pc.PSFCorrection(psf_sigma_deg=0.1812)
#f_ary, df_rho_div_f_ary = pc_inst.f_ary, pc_inst.df_rho_div_f_ary

n.configure_for_scan(f_ary, df_rho_div_f_ary, nexp=1)

#Run the scan
n.perform_scan(nlive=750)

n.load_scan()

an = dnds_analysis.Analysis(n)

an.make_triangle()

plt.figure(figsize=[6,5])

an.plot_source_count_median('gce_np',smin=0.01,smax=1000,nsteps=1000,color='firebrick',spow=0,label='GCE PS')
an.plot_source_count_band('gce_np',smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],color='firebrick',alpha=0.15,spow=0)
an.plot_source_count_band('gce_np',smin=0.01,smax=1000,nsteps=1000,qs=[0.025,0.5,0.975],color='firebrick',alpha=0.1,spow=0)

plt.xlim([5e-11,5e-9])
plt.ylim([2e5,2e9])
plt.yscale('log')
plt.xscale('log')
plt.xlim(F[-1]*1e-4,F[0]*1e4)
plt.tick_params(axis='x', length=5, width=2, labelsize=18)
plt.tick_params(axis='y', length=5, width=2, labelsize=18)
plt.ylabel('$dN/dF$ [counts$^{-1}$cm$^2$ s deg$^{-2}$]', fontsize=18)
plt.xlabel('$F$  [counts cm$^{-2}$ s$^{-1}$]', fontsize=18)
plt.title('Galactic Center NPTF', y=1.02)
plt.legend(fancybox=True)
plt.show()
