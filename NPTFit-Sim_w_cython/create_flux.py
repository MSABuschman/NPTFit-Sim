###############################################################################
# create_flux.py
###############################################################################
#
# Script that takes a given number of sources and draws a flux for each source 
# based on a user defined source count distribution.
#
###############################################################################

import numpy as np
import pdf_sampler

import matplotlib.pyplot as plt

## Global Variables ##
i = 1 #index used for facilitating recursion
hold = 1 #hold running values of coefficent for dN/dF

def dNdF(n,F,f):
    """ Defines the source count distribution (SCD).

            :param n: Python array of index
            :param F: Python array of flux breaks
            :param f: input flux
    """
    #load in the global variables to aid in recursion
    global i,hold

    #Make sure arrays have correct data type
    n = n.astype(float)
    F = F.astype(float)

    #Assure that f is indeed a float
    f = float(f)

    #Determine if we have reached the end of the SCD
    if f < F[len(F) - 1]:
        co = 1
        k = 0
        while k < len(F) - 1:
            co *= np.power(F[k + 1]/F[k],-1*n[k+1])
            k += 1
        x = np.power( f / F[len(F) - 1], -1*n[len(n) - 1])
        return co * x
    #See if we are above the upper most break
    elif f > F[0]:
        return np.power( f / F[0], -1*n[0])
    #Deal with all breaks in the middle
    elif f > F[i]:
        return hold * np.power( f / F[i-1], -1*n[i])
    #Update the index, coefficent(hold), and continue recursion.
    else:
        hold *= np.power( F[i] / F[i-1], -1*n[i])
        i += 1
        return dNdF(n,F,f) 

def run(N,n,F,lobo=-1,upbo=-1):
    """ Runs program and returns list of N flux values drawn from SCD using
        inversion sampeling.

            :param N: number of flux values to be drawn
            :param n: numpy array of index
            :param F: numpy array of flux breaks
            :param lobo: how far to go past lowest break
            :param upbo: how far to go past highest break 
    """
    # Read in the global variables
    global i, hold

    #Determine the upper and lower bounds, if not specified go 4 orders up/down
    if lobo == -1:
        lobo = F[0]*1e-4
    if upbo == -1:
        upbo = F[-1]*1e4
    #Use numpy to sample the distribution, also accounts for log-space w/ Jacob.
    print "Sampling the source count distribution."
    f = np.logspace(np.log10(lobo),np.log10(upbo),1e5)
    dv = f[1:] - f[:-1]
    dv = np.append(dv,dv[-1])
    #For each sample of dist, calc corr. dNdF value.
    j = 0
    #Empty Python list to hold dNdF values
    DNDF = []
    while j < len(f):
        DNDF.append(dNdF(n,F,f[j]))
        #Reset the two global variables
        i = 1
        hold = 1
        j += 1
    #Convert Python list to numpy array.
    DNDF = np.array(DNDF)
    #Draw N flux values from the SCD using inversion sampling
    pdf = pdf_sampler.PDFSampler(f, dv * DNDF)

    return pdf(N)
