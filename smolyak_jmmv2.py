"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebyshev
polynomials to interpolate and approximate the data.

Method based on Judd, Maliar, Maliar, Valero 2013 (W.P)

Author: Chase Coleman
"""

import numpy as np
import numpy.linalg as la
from numpy.polynomial import chebyshev
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
import time as time
import cProfile

class Smoly_JMMV(object):
    """
    This class currently takes a dimension and a degree of polynomial
    and builds the Smolyak Sparse grid.  We base this on the work by
    Judd Maliar Maliar and Valero.  Hope to obtain speed ups beyond
    what they achieved.

    Attributes:
    ------------------
    d : scalar : integer
        This is the dimension of grid that you are building

    mu : scalar : integer
        mu is a parameter that defines the fineness of grid that we
        want to build


    Methods:
    --------
    grid_build : This method builds the sparse grid

    """
    def __init__(self, d, mu):
        """
        Parameters
        ----------
        d : scalar : integer
            This is the dimension of grid that you are building

        mu : scalar : integer
            mu is a parameter that defines the fineness of grid that we
            want to build
        """
        self.d = d
        self.mu = mu

        if mu<1:
            raise ValueError('The parameter mu needs to be > 1.')

        if d<=1:
            raise ValueError('You are trying to build a one dimensional\
                             grid.')

    def _find_A_n(self, n):
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
        Sn = self._find_S_n(n)
        A_chain = []

        # Need a for loop to extract remaining elements
        for seq in xrange(2, n):
            num = Sn.size
            # Need odd indices in python because indexing starts at 0
            A_chain.append(Sn[range(1, num, 2)])
            Sn = Sn[np.arange(0, num, 2)]

        A_chain.append([-1, 1])
        A_chain.append([0])
        A_chain.reverse()

        return A_chain

    def _find_S_n(self, n):
        """
        This method finds the element S_n for the Chebyshev Extrema
        """

        if n==1:
            return np.array([0])

        # Apply the necessary transformation to get the nested sequence
        m_i = 2**(n-1) + 1

        # Create an array of values that will be passed in to calculate
        # the set of values
        comp_vals = np.arange(1., m_i + 1.)

        # Values are - cos(pi(j-1)/(n-1)) for j in [1, 2, ..., n]
        vals = -1. * np.cos(np.pi*(comp_vals - 1.)/(m_i-1.))
        vals[np.where(np.abs(vals) < 1e-14)] = 0.0

        return vals

    def build_sparse_grid(self):
        """
        This method builds a grid for the object
        """
        d = self.d
        mu = self.mu

        # Need to capture up to value mu + 1 so in python need mu+2
        possible_values = range(1, mu + 2)

        An = self._find_A_n(mu + 1)

        points = []
        inds = []

        # check for all indices in range we want
        for max in range(mu+1):
            inds.extend([v for v in product(range(1, mu+2), repeat=d) if sum(v)==max])

        for el in inds:
            temp = [An[i-1] for i in el]
            # Save these indices that we iterate through gecause
            # we need them for the chebyshev polynomial combination
            # inds.append(el)
            points.extend(list(product(*temp)))

        grid = np.array(points)

        self.grid = grid

        return self.grid


    def create_phi_grid(self):
        """
        Makes a matrix as presented in section 3 of JMMV paper.
        """
        d = self.d
        mu = self.mu
        indix = self.inds
        achain = self.Achain

        phi_chain = []
        phi_chain.append(1)


        for el in indix:
            print(chebyshev.chebval(achain[0], [1]))





    def cheb_poly_seq(self):
        """
        Makes a sequence of chebyshev polynomials evaluated at x
        """



def check_points(d, mu):
    if abs(mu - 1) < 1e-14:
        return 2*d + 1

    if abs(mu - 2) < 1e-14:
        return 1 + 4*d + 4*d*(d-1)/2.

    if abs(mu - 3) < 1e-14:
        return 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6.

d = 12
mu = 2
s = Smoly_JMMV(d, mu)
cProfile.run("s.build_sparse_grid()")



