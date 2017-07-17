###############################################################################
# int_src_cnt_dist.py
###############################################################################
#
# Program to integrate broken power law with n breaks for a given source count
# distribution. Returns number of sources after a Poisson draw.
#
###############################################################################

import numpy as np

## GLOBAL VARIABLES ##
coef = 1 #running coefficant
term_arr = None #array to store terms of integration
i = 0 #place holder for recursion steps

def recur(n,F):
    """ Function to recursively do the integration required for calculating
        the number of sources from SCD.

            :param n: numpy array of indexes
            :param F: numpy array of flux breaks
    """
    #Make sure arrays have correct data type
    n = n.astype(float)
    F = F.astype(float)
    #Read in global variables
    global coef, term_arr,i

    #Determine if we are at last term in recursion
    if i >= len(n):
        return term_arr
    #Otherwise, contiue
    else:
        pass

    #Calc. coef. common to all terms
    c = np.power( float(1 - n[i]), -1 )

    #Determine if first term
    if i == len(n) - 1:
        term_arr[i] = 1* c * F[i - 1]
    #Determine if last term
    elif i == 0:
        term_arr[i] = - 1 * c * F[i]
    #Otherwise, it must be middle terms
    else:
        coef = coef * np.power(F[i] / F[i - 1], -1*n[i])
        coef = coef * c
        f = np.power( F[i], 1 - n[i])/np.power(F[i - 1],-1*n[i])
        term_arr[i] = coef * f - F[i - 1]
    #index i to contiue recursion
    i += 1
    #returns int. terms when i is larger than len(n) - 1
    return recur(n,F)

def run(n,F,A,temp):
    """ Function to run calc. of number of sources for a template give a source
        count distribution.

            :param n: numpy array of indexes
            :param F: numpy array of Flux breaks
            :param A: log10 norm. coef. for the SCD
            :param temp: String with name of template     
    """
    #Read in the global variable
    global term_arr

    #Make empty numpy array of len(n)
    term_arr = np.zeros(len(n))

    #Start recursion function, and sum terms
    #Returns final integral function, i.e. number of sources
    runInt = recur(n,F)

    #Sum terms to calculate value of integral
    intgrl =  sum(runInt)

    #Combine normilization and integral. A is in log space, is converted to lin.
    coef = np.power(10.0,A) * intgrl

    #Load in the template map
    gce_map = np.load(str(temp))

    #Sum the pixels of the map
    gce_sum = np.sum(gce_map)

    #Multiply through to get total expected sources from template
    exp_num = coef * gce_sum

    pois_draw = np.random.poisson(exp_num)

    #Do a Poisson draw, where exp_num is the mean
    print "Number of sources from Poisson draw: " + str(pois_draw)
    return pois_draw
