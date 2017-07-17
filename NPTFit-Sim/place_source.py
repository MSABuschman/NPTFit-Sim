###############################################################################
# place_source.py
###############################################################################
#
# Determines N source positions based on the user-input tempalte. This is done
# using iversion sampling method.
#
###############################################################################

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def run(N,temp):
    """
        Creates and returns array of positions for each source based on the 
        template given by user. Returns array of length N with Healpix indexes.

            :params N: number of source positions to draw
            :params temp: String for the file name of template
    """

    print "Determining locations for sources."

    #Load in the Healpix template map
    temp = np.load(temp)

    #Determine the length of the Healpix array 
    L = len(temp)

    #Find the integrated value of template PDF for normilization
    MAX = np.max(np.cumsum(temp))

    #For each source, draw a weighted locations from Healpix map
    i = 0
    #Array to hold Healpix index for each source
    hold = np.zeros(N)
    while i < N:
        #Randomly select source position. Weighted based on template.
        v = np.random.choice(L,1,p=(temp/MAX))
        hold[i] = v[0]
        i += 1
    #Return the array
    return hold
