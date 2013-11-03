"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebyshev
polynomials to interpolate and approximate the data.

Method based on Judd, Maliar, Maliar, Valero 2013 (W.P)

Author: Chase Coleman and Spencer Lyon
"""
import sys
import numpy as np
import numpy.linalg as la
from numpy.polynomial import chebyshev
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations_with_replacement, permutations
import time as time
from dolo.numeric.interpolation import smolyak as smm
from smoly_cyutils import *


def find_A_n(n):
    """
    This method finds all of the unidimensional disjoint sets
    that we will use to construct the grid.  It improves on
    past algorithms by noting that A_{n} = S_{n}[evens] except for
    A_1 = {0} and A_2 = {-1, 1}. Additionally, A_{n} = A_{n+1}[odds]
    This prevents the calculation of these nodes repeatedly.  Thus
    we only need to calculate biggest of the S_n's to build the
    sequence of A_n's
    """

    # # Start w finding the biggest Sn(We will subsequently reduce it)
    Sn = find_S_n(n)
    A_chain = {}
    A_chain[1] = [0.]
    A_chain[2] = [-1., 1.]

    # Need a for loop to extract remaining elements
    for seq in xrange(n, 2, -1):
        num = Sn.size
        # Need odd indices in python because indexing starts at 0
        A_chain[seq] = tuple(Sn[range(1, num, 2)])
        # A_chain.append(list(Sn[range(1, num, 2)]))
        Sn = Sn[range(0, num, 2)]

    return A_chain


def find_S_n(n):
    """
    This method finds the element S_n for the Chebyshev Extrema
    """

    if n==1:
        return np.array([0.])

    # Apply the necessary transformation to get the nested sequence
    m_i = 2**(n-1) + 1

    # Create an array of values that will be passed in to calculate
    # the set of values
    comp_vals = np.arange(1., m_i + 1.)

    # Values are - cos(pi(j-1)/(n-1)) for j in [1, 2, ..., n]
    vals = -1. * np.cos(np.pi*(comp_vals - 1.)/(m_i-1.))
    vals[np.where(np.abs(vals) < 1e-14)] = 0.0

    return vals



def permute(a):
    """
    Creates all unique combinations of a list a that is passed in.
    Function is based off of a function written by John Lettman:
    TCHS Computer Information Systems.  My thanks to him.
    """

    a.sort() # Sort.

    ## Output the first input sorted.
    yield list(a)

    i = 0
    first = 0
    alen = len(a)

    ## "alen" could also be used for the reference to the last element.

    while(True):
        i = alen - 1

        while(True):
            i -= 1 # i--

            if(a[i] < a[(i + 1)]):
                j = alen - 1

                while(not (a[i] < a[j])): j -= 1 # j--

                a[i], a[j] = a[j], a[i] # swap(a[j], a[i])
                t = a[(i + 1):alen]
                t.reverse()
                a[(i + 1):alen] = t

                # Output current.
                yield list(a)

                break # next.

            if(i == first):
                a.reverse()

                yield list(a)
                return


def check_points(d, mu):
    if mu == 1:
        return 2*d + 1

    if mu == 2:
        return 1 + 4*d + 4*d*(d-1)/2.

    if mu == 3:
        return 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6.


# my_args = sys.argv[1:]

# if __name__ == '__main__':

#     for d, mu in [(20, 2), (10, 3), (20, 3)]:
#         s = Smoly_grid(d, mu)
#         print(s.build_sparse_grid().shape, check_points(d, mu))

# if 'prof' in my_args or 'profile' in my_args:
    import cProfile
    cProfile.run("s.build_sparse_grid()")
    cProfile.run("smm.smolyak_grids(d, mu)")





