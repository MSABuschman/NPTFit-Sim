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
import make_map as mm

import numpy as np

def run(n,F,A,temp,exp,psf_r,name,rsamp=1):
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
            :param rsamp: set to 1 for rejection sampling, otherwise samples
                pixels in template like a pdf.

            :returns: saves numpy array of simulated map
    """
    # Int. SCD to find mean couts, Poisson draws for # of sources in template
    num_src = iscd.run(n,F,A,temp)

    # Draws fluxes for each source from the CSD
    flux_arr = cf.run(num_src,n,F)

    # Generate simulated counts map
    map_arr = np.asarray(mm.run(num_src,flux_arr,temp,exp,psf_r,rsamp))

    #Save the file as an .npy file
    np.save(str(name) + ".npy",map_arr.astype(np.int))
    print "Done simulation."
