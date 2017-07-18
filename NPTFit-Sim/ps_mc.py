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
    dist_arr = sd.run(pos_arr,flux_arr,temp,exp,psf_r)

    #Save the file as an .npy file
    np.save(str(name) + ".npy",dist_arr)
    print "Done simulation."
