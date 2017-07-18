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
import matplotlib.pyplot as plt

#Global variables
th_arr = None #The theta coords of all pixels in the map
ph_arr = None #The phi coords of all pixels in the map
NSIDE = None #The NSIDE of the Healpix map, found in find_all_coords
pix_dist = None #Distance of each pixel from the source of interest
map_arr = None #Array that is used to build the simulated counts map

def find_all_coords(temp):
    """ This goes through whole template and determines the theta and phi for
        each point in the Healpix grid.

            :param temp: String of template file name
    """
    #Load in the three global variables
    global th_arr, ph_arr, NSIDE

    #Determine the NSIDE of the Healpix map from length of template array
    NSIDE = hp.npix2nside(len(temp))

    #Make two empty arrays to fill with theta and phi values
    th_arr = np.zeros(len(temp))
    ph_arr = np.zeros(len(temp))

    #Loop over whole template and calc. theta and phi for each pixel
    i = 0
    while i < len(temp):
        th_arr[i],ph_arr[i] = hp.pix2ang(NSIDE,i)
        i += 1
    return None

def ang_dist(src_hp_index):
    """ This calculates the angular distance from a given Healpix pixel to all
        other pixels in degrees.

            :param src_hp_index: integer of Healpix index
    """
    #Load global variables
    global th_arr, ph_arr, NSIDE, pix_dist

    #Calculate theta and phi for pixel of interest
    th,ph = hp.pix2ang(NSIDE,src_hp_index)

    #Make an empty numpy array to hold pixel distances from source location.
    pix_dist = np.zeros(len(th_arr))

    #Loops over all pixels and calculates distances from source position
    i = 0
    while i < len(th_arr):
        #Calc. the angular distance from the two points with some trig.
        cos_term = np.cos(th) * np.cos(th_arr[i])
        sin_term = np.sin(th) * np.sin(th_arr[i]) * np.cos( ph - ph_arr[i])
        #Ensure we remain in domain of acos, add. occasionaly goes out of domain
        cs = cos_term + sin_term
        if cs >= 1.0:
            cs = 1.0
        elif cs <= -1.0:
            cs = -1.0
        #Convert from radians to degrees
        dist = abs(np.arccos(cs)) * 180./np.pi
        if dist <= 90:
            pix_dist[i] = dist
        else:
            pix_dist[i] = 0
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
    #Calc. PSF based on the distance array for each source
    PSF_val = np.abs(psf_r(pix_dist))
    #Find the integrated value of the PSF for normilization
    norm = np.max(np.cumsum(PSF_val))
    #Norm the PSF, multiply by flux and exposure at pixel
    hold = (PSF_val / norm) * (flux * EXP_map)
    #Do Poisson draw for every pixel on map to get counts, add to running
    #map of the simulated sky
    hold = np.random.poisson(hold)
    map_arr = map_arr + hold
    return None

def run(pos_arr,flux_arr,temp,exp,psf_r):
    """ Implements the creation of PSF map for each source. Returns a 2D array
        of PSF maps for each source.

            :param pos_arr: array of Healpix indexes
            :param temp: String of templates file name
    """
    #Load in the global variables
    global th_arr, ph_arr, NSIDE, pix_dist, map_arr

    #Load the exposure map
    EXP_map = np.load(exp)

    map_arr = np.zeros(len(EXP_map))

    #Make sure the position array is of type int (sould be Healpix index)
    pos_arr = pos_arr.astype(int)

    #Load in the correct template
    temp = np.load(temp)

    #Calculate theta and phi for all pixels
    find_all_coords(temp)

    #Loop, starting at second source, over whole position array, calculating
    #distance arrays for all sources. These 1-D arrays are stacked vertically
    #and returned.

    print "Simulating counts map."

    i = 0
    while i < len(pos_arr):
        #Find the angular distance for all pixels from first point
        ang_dist(pos_arr[i])
        #Stack the arrays on top of eachother
        #two_d_arr = np.vstack((two_d_arr,pix_dist))
        create_map(flux_arr[i],EXP_map,psf_r)
        #print "Done source " + str(i + 1) + " out of " + str(len(pos_arr))
        i += 1
    return map_arr
