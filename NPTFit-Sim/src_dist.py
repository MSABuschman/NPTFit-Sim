###############################################################################
# src_dist.py
###############################################################################
#
# For a list of source postions on Healpix map, returns a 2-D array of distances
# of every pixel from each source. 
#
###############################################################################

import healpy as hp
import numpy as np

import rej_samp as rs

# Global variables
th_arr = None # The theta coords of all pixels in the map
ph_arr = None # The phi coords of all pixels in the map
pix_dist = None # Distance of each pixel from the source of interest
map_arr = None # Array that is used to build the simulated counts map

def ang_dist(th,ph):
    """ This calculates the angular distance from a given Healpix pixel to all
        other pixels in degrees.

            :param th: source theta position in radians
            :param ph: source phi position in radians
    """
    # Load global variables
    global th_arr, ph_arr, pix_dist

    # Make an empty numpy array to hold pixel distances from source location.
    pix_dist = np.zeros(len(th_arr))

    # Loops over all pixels and calculates distances from source position
    i = 0
    while i < len(th_arr):
        # Calc. the angular distance from the two points with some trig.
        cos_term = np.cos(th) * np.cos(th_arr[i])
        sin_term = np.sin(th) * np.sin(th_arr[i]) * np.cos( ph - ph_arr[i])
        # Ensure we remain in domain of acos, add. occasionaly goes out of domain
        cs = cos_term + sin_term
        if cs >= 1.0:
            cs = 1.0
        elif cs <= -1.0:
            cs = -1.0
        # Make sure no source is greater than pi radians away
        dist = abs(np.arccos(cs))
        if dist <= np.pi:
            pix_dist[i] = dist
        else:
            pix_dist[i] = 2*np.pi - dist
        i += 1
    return None

def create_map(flux,EXP_map,psf_r):
    """ Reads in array of distances from the source and using the user defined
        PSF constructs a simulated counts map by doing Poisson draws.

            :param flux: value of flux for source
            :param EXP_map: Healpix map of exposure for whole sky
            :param psf_r: user defined point spread function
    """
    global pix_dist, map_arr
    # Calc. PSF based on the distance array for each source
    PSF_val = np.abs(psf_r(pix_dist))
    # Find the integrated value of the PSF for normilization
    norm = np.max(np.cumsum(PSF_val))
    # Norm the PSF, multiply by flux and exposure at pixel
    hold = (PSF_val / norm) * (flux * EXP_map)
    # Add to running sum of counts on sky
    map_arr = map_arr + hold
    return None

def run(N,flux_arr,temp,EXP_map,psf_r):
    """ Implements the creation of PSF map for each source. Returns a 2D array
        of PSF maps for each source.

            :param N: number of sources
            :param flux_arr: array source fluxes
            :param temp: numpy array for template
            :param EXP_map: numpy array of exposure map
            :param psf_r: user defined point spread function
    """
    # Load in the global variables
    global th_arr, ph_arr, pix_dist, map_arr

    # Determine the NSIDE of the Healpix map from length of template array
    NSIDE = hp.npix2nside(len(temp))    

    # Find all theta and phi values for pixels in template map
    th_arr, ph_arr = hp.pix2ang(NSIDE,range(len(temp)))

    # Create array to use when sim
    map_arr = np.zeros(len(EXP_map))

    print "Simulating counts map."

    # For each source, create a counts map. Keeps running total of counts from
    # all sources.
    i = 0
    while i < N:
        # Determine a source position, in terms of theta and phi, using
        # rejection sampling. 
        th,ph = rs.reject(temp)
        # Find the angular distance for all pixels from this position
        ang_dist(th,ph)
        # Generate simulated counts map for source
        create_map(flux_arr[i],EXP_map,psf_r)
        #print "Done source " + str(i + 1) + " out of " + str(N)
        i += 1
    # Do Poisson draw for every pixel on map to get counts, add to running
    # map of the simulated sky
    map_arr = np.random.poisson(map_arr)
    return map_arr
