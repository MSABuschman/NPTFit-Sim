###############################################################################
# ps_mc.py
###############################################################################
#
# Program that does point source Monte Carlo based off of user defined source
# count distribution, template, exposure map, and user defined point source
# function. Creates .npy file with results of simulation.
#
###############################################################################

import int_src_cnt_dist as iscd
import create_flux as cf
import place_source as ps
import src_dist as sd

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def run(n,F,A,temp,exp,psf_r,name):
    """ Brings together serveral programs to run point source Monte Carlo by
        reading in template, source count distribution parameters, exposure 
        map, and the user defined PSF.

            :param n: numpy array of index values
            :param F: numpy array of flux break values
            :param A: float of log10 norm for the SCD
            :param temp: String for name of template
            :param exp: String for name of exposure map
            :param psf_r: User defined point spread function
            :param name: String for the name of output .npy file
    """
    #int. SCD to find mean couts, Poisson draws for # of sources in template
    num_src = iscd.run(n,F,A,temp)

    #Draws fluxes for each source from the CSD
    flux_arr = cf.run(num_src,n,F)

    #Make array of source positions based on template
    pos_arr = ps.run(num_src,temp)

    #Array with distances from point to all other pixels
    dist_arr = sd.run(pos_arr,temp)

    #Load the exposure map
    EXP_map = np.load(exp)

    #Takes first source, norm. its PSF, mulitply by flux and exposure, then
    #do a Poisson draw. Add to hold array wich is a running version of the final
    #map simulation.
    PSF_val = psf_r(dist_arr[0])
    norm = np.max(np.cumsum(PSF_val))
    hold = (PSF_val / norm) * (flux_arr[0] * EXP_map[int(pos_arr[0])])
    hold = np.random.poisson(hold)

    #Norms PSF, multiply by flux and exposure, do Poisson draw. Add result for
    #each source to a running map.
    i = 1
    while i < len(dist_arr):
        #Calc. PSF based on the distance array for each source
        PSF_val = psf_r(dist_arr[i])
        #Find the integrated value of the PSF for normilization
        norm = np.max(np.cumsum(PSF_val))
        #Norm the PSF, multiply by flux and exposure at pixel
        tmp = (PSF_val / norm) * (flux_arr[i] * EXP_map[int(pos_arr[i])])
        #Do Poisson draw for every pixel on map to get counts, add to running
        #map of the simulated sky
        hold = hold + np.random.poisson(tmp)
        i += 1

    #Save the file as an .npy file
    np.save(str(name) + ".npy",hold)

    hp.mollview(hold, title="PS Monte Carlo")
    plt.show()
