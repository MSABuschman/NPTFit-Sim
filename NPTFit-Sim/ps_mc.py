###############################################################################
# ps_mc.py
###############################################################################
#
# Program that does point source Monte Carlo based off of user defined source
# count distribution, template, exposure map, and user defined point source
# function. Can save result of simulation to .npy file.
#
###############################################################################

import int_src_cnt_dist as iscd
import create_flux as cf
import make_map as mm

import numpy as np

def run(n,F,A,temp,exp,psf_r,Ebins,Edep,Eparam,name="map",save=False):
    """ Brings together serveral programs to run point source Monte Carlo by
        reading in template, source count distribution parameters, exposure 
        map, and the user defined PSF.

            :param n: numpy array of index values
            :param F: numpy array of flux break values
            :param A: float of log10 norm for the SCD
            :param temp: HEALPix numpy array of template
            :param exp: HEALPix numpy array of exposure map
            :param psf_r: user defined point spread function
            :param Ebins: numpy array of energy bin edges
            :param Edep: user defined energy dependence
            :param name: string for the name of output .npy file
            :param save: option to save map to .npy file

            :returns: HEALPix format numpy array of simulated maps
    """
    # Int. SCD to find mean couts, Poisson draws for # of sources in template
    num_src = iscd.run(n,F,A,temp)

    # Draws fluxes for each source from the SCD
    flux_arr = cf.run(num_src,n,F)

    # Generate simulated counts map
    maps = mm.run(num_src,flux_arr,temp,exp,psf_r,Ebins,Edep,Eparam)
    map_arr = np.asarray(maps[:])

    # Save the file as an .npy file
    if save:
        np.save(str(name) + ".npy",map_arr.astype(np.int32))

    print "Done simulation."

    return map_arr.astype(np.int32)
