"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebyshev
polynomials to interpolate and approximate the data.

Method based on Malin, Krueger, Kubler 2010 and
                Judd, Maliar, Maliar, Valero 2013 (W.P)

Author: Chase Coleman
"""

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
import time as time

class Smolyak(object):
    """
    This class currently takes a dimension and a degree of polynomial
    and builds the Smolyak Sparse grid.

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

    def _m_seq(self, vals):
        """
        Pass in a value and it will calculate the ith element of the
        sequence of m's as described in paper.
        """

        # Create vec of ones in same shape as vals.
        # Everywhere vals=1 we want to put a 1, so saves a step.
        ret_vals = np.ones(vals.shape)

        # Everywhere else it should be 2**(i-1) + 1
        ret_vals[np.where(vals > 1)] = 2**(vals[np.where(vals > 1)] - 1) + 1

        return ret_vals

    def _g_func(self, n):
        """
        This function returns the set G.  It corresponds to the
        Chebyshev extrema.
        """

        if abs(n - 1) < 1e-10:
            g_set = np.array([0])
        else:
            g_set = np.cos(np.pi*(np.arange(n))/(n-1.))
            # make sure they are 0 where I want them 0
            g_set[np.where(abs(g_set) < 1e-14)] = 0.

        return g_set

    def find_points(self):
        """
        For integers between 1 and d find the sets of d integers that
        add to mu + d.

        i.e. mu=1, d=3
        {1, 1, 2}, {1, 2, 1}, {2, 1, 1}
        """

        d = self.d
        mu = self.mu

        # We want all possible values up to dimension + 1
        poss_vals = np.arange(mu + 1) + 1

        # Create empty list to fill with the points
        points = []

        # Use this for loop to create the list of points
        for element in product(poss_vals, repeat=d):
            # We only want the ones that sum to d + mu
            if sum(element) == d + mu:

                # Transform the element i to m(i)
                elements = self._m_seq(np.array(element))

                # Now take these points and make the set g
                temp_points = map(self._g_func, elements)

                # Put those points into the list of points
                points.append(list(product(*temp_points)))

        return points


    def build_grid(self):
        """
        This method builds a d dimensional grid to be used for
        interpolation.
        """
        d = self.d
        mu = self.mu

        # Call the find_points function to find the elements that we
        # need to build the grid
        points = self.find_points()
        point_set = set()

        # Put each of the points into a set so that we don't have
        # duplicate points
        for elem in points:
            for temp in elem:
                point_set.add(temp)

        # Now put them in an array.  Each row is a d-dimensional vector
        # that represents the point.
        point_array = np.array(list(point_set))

        self.grid = point_array

#----------------------------------------------------------------------#
# Example
#----------------------------------------------------------------------#

def check_points(d, mu):
    if abs(mu - 1) < 1e-14:
        return 2*d + 1

    if abs(mu - 2) < 1e-14:
        return 1 + 4*d + 4*d*(d-1)/2.

    if abs(mu - 3) < 1e-14:
        return 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6.


d = 3
mu = 10
# start = time.time()
smolyak_grid = Smolyak(d, mu)
smolyak_grid.build_grid()
# elapsed = start - time.time()
# point_array = smolyak_grid.grid

# print(smolyak_grid.grid.shape)
# print(check_points(d, mu), "Took %.5f" %elapsed)


# Plot to show the where points go for 3-D case

# fig = plt.figure(1)
# fig.suptitle("Smolyak Grid on 3 Dimensions with $\mu =$ %d"  %smolyak_grid.mu)
# ax1 = fig.add_subplot(221)
# ax1.set_title("X-Y Grid")
# ax1.plot(point_array[:, 0], point_array[:, 1], 'ko')
# ax1.set_ylim((-1.5, 1.5))
# ax1.set_xlim((-1.5, 1.5))
# ax2 = fig.add_subplot(222)
# ax2.set_title("Y-Z Grid")
# ax2.plot(point_array[:, 1], point_array[:, 2], 'ko')
# ax2.set_ylim((-1.5, 1.5))
# ax2.set_xlim((-1.5, 1.5))
# ax3 = fig.add_subplot(223)
# ax3.set_title("X-Z Grid")
# ax3.plot(point_array[:, 0], point_array[:, 2], 'ko')
# ax3.set_ylim((-1.5, 1.5))
# ax3.set_xlim((-1.5, 1.5))
# ax4 = fig.add_subplot(224, projection='3d')
# ax4.set_title("3-D Projection")
# surf = ax4.scatter(point_array[:, 0], point_array[:, 1], point_array[:, 2])
# ax4.set_ylim((-1.5, 1.5))
# ax4.set_xlim((-1.5, 1.5))
# ax4.set_zlim((-1.5, 1.5))
# plt.show()

