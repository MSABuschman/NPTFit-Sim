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
coef = None # running coefficant
term_arr = None # array to store terms of integration
i = 0 # place holder for recursion steps

def recur(n,F):
    """ Function to recursively do the integration required for calculating
        the number of sources from SCD.

            :param n: numpy array of indexes
            :param F: numpy array of flux breaks
    """
    # Make sure arrays have correct data type
    n = n.astype(float)
    F = F.astype(float)

    # Read in global variables
    global coef, term_arr, i

    # Determine if we are at last term in recursion
    if i >= len(n):
        return term_arr
    # Otherwise, contiue
    else:
        pass

    # Calc. coef. common to all terms
    c = np.power( float(1 - n[i]), -1 )

    # Determine if first term (index above the highest break)
    if i == 0:
        term_arr[i] = - 1 * c * F[i]
    # Determine if last term (index below lowest break)
    elif i == len(n) - 1:
        term_arr[i] = coef[i] * c * F[i - 1]
    # Otherwise, it must be the middle terms
    else:
        coef[i:] *= np.power(F[i] / F[i - 1], -1*n[i])
        term_arr[i] = c * (coef[i-1] * F[i-1] - coef[i] * F[i])
    # Index i to contiue recursion
    i += 1
    # Returns int. terms when i is larger than len(n) - 1
    return recur(n,F)

def run(n,F,A,temp_map):
    """ Function to run calc. of number of sources for a template give a source
        count distribution.

            :param n: numpy array of indexes
            :param F: numpy array of Flux breaks
            :param A: log10 norm. coef. for the SCD
            :param temp_map: numpy array of template

            :returns: number of sources from Poisson draw   
    """
    # Read in the global variables
    global term_arr, coef, i

    # Make empty numpy array of len(n), and ones for coef
    term_arr = np.zeros(len(n))
    coef = np.ones(len(n))

    # Start recursion function, and sum terms
    # Returns final integral function, i.e. number of sources
    runInt = recur(n,F)

    # Sum terms to calculate value of integral
    intgrl = sum(runInt)

    # Combine normilization and integral. A is in log space, is converted to lin
    coef = np.power(10.0,A) * intgrl

    # Sum the pixels of the map
    temp_sum = np.sum(temp_map)

    # Multiply through to get total expected sources from template
    exp_num = coef * temp_sum

    # Do a Poisson draw, where exp_num is the mean
    pois_draw = np.random.poisson(exp_num)

    print "Number of sources from Poisson draw: " + str(pois_draw)

    # Reset Global variables in case of multiple calls of function
    coef = None
    term_arr = None 
    i = 0

    return pois_draw
