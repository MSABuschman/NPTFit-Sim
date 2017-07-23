###############################################################################
# make_map.pyx
###############################################################################
#
# Given a number of sources, array of fluxes for each source, template & 
# exposure maps, and user defined PSF, simulates and returns a counts map. The
# user has the option to draw source positions weighted by the template pixels
# or using a rejection sampling routine. 
#
###############################################################################

import healpy as hp
import numpy as np
cimport numpy as np
cimport cython

import rej_samp as rs
import place_source as ps

# Call in cython functions
cdef extern from "math.h":
    double cos(double x) nogil
    double sin(double x) nogil
    double acos(double x) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double[::1] ang_dist(double th, double ph, double[::1] th_arr,
                          double[::1] ph_arr):
    """ This calculates the angular distance from a given Healpix pixel to all
        other pixels in degrees.

            :param th: source theta position in radians
            :param ph: source phi position in radians
            :param th_arr: array of theta positions of all pixels in map
            :param ph_arr: array of phi positions of all pixels in map

            :returns: array with distances from source to each pixel in radians
    """

    # Make an empty numpy array to hold pixel distances from source location.
    cdef double[::1] pix_dist = np.zeros(len(th_arr))

    cdef double cos_term, sin_term, cs, dist

    # Loops over all pixels and calculates distances from source position
    cdef int i = 0
    while i < len(th_arr):
        # Calc. the angular distance from the two points with some trig.
        cos_term = cos(th) * cos(th_arr[i])
        sin_term = sin(th) * sin(th_arr[i]) * cos( ph - ph_arr[i])
        # Ensure we remain in domain of acos, float addition may go just out of
        # the domain.
        cs = cos_term + sin_term
        if cs >= 1.0:
            cs = 1.0
        elif cs <= -1.0:
            cs = -1.0
        # Make sure no source is greater than pi radians away
        dist = abs(acos(cs))
        if dist <= np.pi:
            pix_dist[i] = dist
        else:
            pix_dist[i] = 2*np.pi - dist
        i += 1
    return pix_dist

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double[::1] create_map(double flux, double[::1] EXP_map,psf_r,
                            double[::1] pix_dist):
    """ Reads in array of distances from the source and using the user defined
        PSF constructs a simulated counts map by doing Poisson draws.

            :param flux: value of flux for source
            :param EXP_map: Healpix map of exposure for whole sky
            :param psf_r: user defined point spread function
            :param pix_dist: array of pixel distances from source location

            :returns: array of counts from a source prior ro Pois. draw
    """
    # Calc. PSF based on the distance array for each source
    cdef double [::1] PSF_val = psf_r(np.asarray(pix_dist))
    # Multiply PSF by Jacobian factor
    PSF_val *= np.asarray(pix_dist)

    # Find the integrated value of the PSF for normilization
    cdef double norm = np.sum(PSF_val)

    # Norm the PSF, multiply by flux and exposure at pixel to get counts
    cdef double[::1] hold = (np.asarray(PSF_val) / norm)
    hold *= (flux * np.asarray(EXP_map))

    return hold

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef long[::1] sum_map(int N,double[::1] flux_arr,double[::1] temp,
                       double[::1] EXP_map,psf_r, int rsamp):
    """ Implements the creation of PSF map for each source. Returns a 2D array
        of PSF maps for each source.

            :param N: number of sources
            :param flux_arr: array source fluxes
            :param temp: numpy array for template
            :param EXP_map: numpy array of exposure map
            :param psf_r: user defined point spread function
            :param rsamp: set True for rejection sampling

            :returns: array of simulated counts map
    """
    # Load in the global variables
    #global th_arr, ph_arr, pix_dist, map_arr

    # Determine the NSIDE of the Healpix map from length of template array
    cdef int NSIDE = hp.npix2nside(len(temp))
    cdef double[::1] th_arr, ph_arr, pix_dist
    cdef long[::1] pix_loc   
    cdef int ploc
    # Find all theta and phi values for pixels in template map
    th_arr, ph_arr = hp.pix2ang(NSIDE,range(len(temp)))

    # Create array to use when sim
    cdef double[::1] map_arr = np.zeros(len(EXP_map))

    print "Simulating counts map ..."

    # For each source, create a counts map. Keeps running total of counts from
    # all sources.
    cdef int i = 0
    # If no rejection sampleing, sample positions from pixel values in template.
    if rsamp == 0:
        print "Sampling template pixels for source positions ..."
        # Treat template as pdf, draw N weighted pixel positions
        pix_loc = np.asarray(ps.run(N,temp))
        while i < N:
            # Determine the source position based on given pixel
            ploc = pix_loc[i]
            th,ph = np.asarray(hp.pix2ang(NSIDE,ploc))
            # Find the angular distance to all pixels, from this pixel
            pix_dist = np.asarray(ang_dist(th,ph,th_arr,ph_arr))
            # Generate simulated counts map for source
            map_arr += np.asarray(create_map(flux_arr[i],EXP_map,psf_r,pix_dist))
            print "Done source " + str(i + 1) + " out of " + str(N)
            i += 1
    # Otherwise, do rejection sampling method.
    else:
        print "Using rejection sampling for source positions ..."
        while i < N:
            # Determine a source position, in terms of theta and phi, using
            # rejection sampling. 
            th,ph = np.asarray(rs.run(temp))
            # Find the angular distance for all pixels from this position
            pix_dist = np.asarray(ang_dist(th,ph,th_arr,ph_arr))
            # Generate simulated counts map for source
            map_arr += np.asarray(create_map(flux_arr[i],EXP_map,psf_r,pix_dist))
            #print "Done source " + str(i + 1) + " out of " + str(N)
            i += 1
    # Do Poisson draw for every pixel on map to get counts, add to running
    # map of the simulated sky
    cdef long[::1] r_map_arr = np.random.poisson(map_arr)
    return r_map_arr

def run(N,flux_arr,temp,EXP_map,psf_r,rsamp):
    """ Python wrapper for simulating counts map from template.

            :param N: number of sources
            :param flux_arr: array source fluxes
            :param temp: numpy array for template
            :param EXP_map: numpy array of exposure map
            :param psf_r: user defined point spread function
            :param rsamp: set to 1 for rejection sampeling

            :returns: array of simulated counts map
    """
    return sum_map(N,flux_arr,temp,EXP_map,psf_r,rsamp)
