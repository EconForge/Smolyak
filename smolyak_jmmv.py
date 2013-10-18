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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
import time as time


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

    def build_grid(self):
        """
        This method builds a grid for the object
        """
        return None

    def _find_A_n(n):
        """
        This method finds all of the unidimensional disjoint sets
        that we will use to construct the grid.  It improves on
        past algorithms by noting that A_{n} = S_{n}[evens] except for
        A_1 = {0} and A_2 = {-1, 1}. Additionally, A_{n} = A_{n+1}[odds]
        This prevents the calculation of these nodes repeatedly.  Thus
        we only need to calculate biggest of the S_n's to build the
        sequence of A_n's
        """
        # n = self.n

        # # Start w finding the biggest Sn(We will subsequently reduce it)
        # Sn = _find_S_n(n)
        A_chain = []

        A_chain.append([0])
        A_chain.append([-1, 1])

        # Need a for loop to extract remaining elements
        for seq in xrange(2, n):
            temp = _find_S_n(n)
            num = Sn.size
            # Need odd indices in python because indexing starts at 0
            A_chain.append(temp[np.arange(1, num, 2)])

        return A_chain

    def _find_S_n(n):
        """
        This method finds the element S_n for the Chebyshev Extrema
        """
        # n = self.n

        if n==1:
            return np.array([0])

        # Apply the necessary transformation to get the nested sequence
        # m_i = 2**(i-1) + 1
        m_i = 2**(n-1) + 1

        # Create an array of values that will be passed in to calculate
        # the set of values
        comp_vals = np.arange(1., m_i + 1.)

        # Values are - cos(pi(j-1)/(n-1)) for j in [1, 2, ..., n]
        vals = -1. * np.cos(np.pi*(comp_vals - 1.)/(m_i-1.))
        vals[np.where(np.abs(vals) < 1e-14)] = 0.0

        return vals



