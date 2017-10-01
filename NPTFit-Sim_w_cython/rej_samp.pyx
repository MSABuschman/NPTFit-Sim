################################################################################
# rej_samp.pyx
################################################################################
#
# Given a template as a numpy array, returns theta and phi position for a source
# using rejection sampeling.
#
################################################################################

import numpy as np
import healpy as hp

cimport numpy as np
cimport cython

DTYPE = np.float

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double[::1] coords():
    """ Returns an array of random theta and phi values (in radians).
    """

    # Create an empty array to hold position
    cdef double [::1] crds = np.zeros(2,dtype=DTYPE)
    # Choose random float from 0 <= x < 1 twice
    cdef double th_ran = np.random.random()
    cdef double ph_ran = np.random.random()
    # Calculate random theta and phi positions
    crds[0] = 2*np.arcsin(np.sqrt(th_ran))
    crds[1] = 2*np.pi*ph_ran 
    return crds

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef double[::1] run(double[::1] temp):
    """ Returns a source position from a give template in terms of theta and phi
        (in radians) using rejection sampling. 

            :params temp: numpy array corresponding to template
    """

    # Determine the NSIDE of the template
    cdef int NSIDE = hp.npix2nside(len(temp))
    cdef double[::1] crds = np.zeros(2)
    cdef int pos = 0
    cdef double rnd = 0.0
    # Make max pixel value in template 1.0 for rejection sampling
    temp = temp / np.max(temp)
    cdef int i = 0
    while i < 1:
        # Grab random source position in terms of theta and phi
        crds = coords()
        # Find coresponding Healpix pixel
        pos = hp.pixelfunc.ang2pix(NSIDE,crds[0],crds[1])
        # Choose random float from 0 <= x < 1 twice
        rnd = np.random.random()
        # If the random number is less that template value, accept coords
        if rnd <= temp[pos]:
            i += 1
    # Returns coords
    return crds
