"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebyshev
polynomials to interpolate and approximate the data.

Method based on Judd, Maliar, Maliar, Valero 2013 (W.P)

Authors: Chase Coleman and Spencer Lyon

"""
import sys
from operator import mul
from itertools import product, combinations_with_replacement
from itertools import chain
import numpy as np
from scipy.linalg import lu
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from util import *


def a_chain(n):
    """
    This method finds all of the unidimensional disjoint sets
    that we will use to construct the grid.  It improves on
    past algorithms by noting that A_{n} = S_{n}[evens] except for
    A_1 = {0} and A_2 = {-1, 1}. Additionally, A_{n} = A_{n+1}[odds]
    This prevents the calculation of these nodes repeatedly.  Thus
    we only need to calculate biggest of the S_n's to build the
    sequence of A_n's

    Parameters
    ----------
    n : scalar : integer
      This is the number of disjoint sets from Sn that this should make

    Returns
    -------
    A_chain : dictionary : lists
      This is a dictionary of the disjoint sets that are made.  They are
      indexed by the integer corresponding
    """

    # # Start w finding the biggest Sn(We will subsequently reduce it)
    Sn = _s_n(n)

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

def _s_n(n):
    """
    This method finds the element S_n for the Chebyshev Extrema
    """

    if n == 1:
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

def smol_inds(d, mu):
    """
    This method finds all of the indices that satisfy the requirement
    that d \leq \sum_{i=1}^d \leq d + \mu.  Once we have these, then
    they can be used to build both the grid and the polynomial
    basis.

    Notes
    =====
    This method sets the attribute smol_inds
    """
    if type(mu)==int:
        max_mu = mu
    else:
        max_mu = int(np.max(mu))

    # Need to capture up to value mu + 1 so in python need mu+2
    possible_values = range(1, max_mu + 2)

    # find all (i1, i2, ... id) such that their sum is in range
    # we want; this will cut down on later iterations
    poss_inds = [el for el in combinations_with_replacement(possible_values, d)
                  if d < sum(el) <= d+max_mu]

    if type(mu)==int:
        true_inds = [[el for el in permute(list(val))] for val in poss_inds]
    else:
        true_inds = [[el for el in permute(list(val)) if all(el <= mu+1)] for val in poss_inds]

    # Add the d dimension 1 array so that we don't repeat it a bunch
    # of times
    true_inds.extend([[[1]*d]])

    tinds = list(chain.from_iterable(true_inds))

    return tinds

def build_grid(d, mu):
    """
    This method builds a grid for the object

    Notes
    =====
    This method sets the attribute grid
    """
    # Get An chain
    An = a_chain(mu + 1)

    points = []

    # Need to get the correct indices

    tinds = smol_inds(d, mu)

    for el in tinds:
        temp = [An[i] for i in el]
        # Save these indices that we iterate through because
        # we need them for the chebyshev polynomial combination
        # inds.append(el)
        points.extend(list(product(*temp)))

    grid = pd.lib.to_object_array_tuples(points).astype(float)

    return grid, tinds

def phi_chain(n):
    """
    Finds the disjoint sets of aphi's that will be used to compute
    which functions we need to calculate
    """

    # First create a dictionary
    aphi_chain = {}

    aphi_chain[1] = [1]
    aphi_chain[2] = [2, 3]

    curr_val = 4
    for i in xrange(3, n+1):
        end_val = 2**(i-1) + 1
        temp = range(curr_val, end_val+1)
        aphi_chain[i] = temp
        curr_val = end_val+1

    return aphi_chain

def poly_inds(d, mu):
    """
    This function builds the indices of the basis polynomials that
    will be used to interpolate.
    """
    inds = smol_inds(d, mu)
    aphi = phi_chain(mu + 1)

    # Bring in polynomials
    # cheb_dict = self.calc_chebvals()

    base_polys = []

    for el in inds:
        temp = [aphi[i] for i in el]
        # Save these indices that we iterate through because
        # we need them for the chebyshev polynomial combination
        # inds.append(el)
        base_polys.extend(list(product(*temp)))

    return base_polys

def build_B(d, mu, grid):
    """
    This function builds the matrix B that will be used to calc
    the interpolation coefficients for a given set of data.

    Notes
    =====
    This method sets the attributes B, B_L, B_U
    """
    Ts = cheby2n(grid.T, m_i(mu + 1))
    base_polys = poly_inds(d, mu)
    n = len(grid)
    B = np.empty((n, n), order='F')
    for ind, comb in enumerate(base_polys):
        B[:, ind] = reduce(mul, [Ts[comb[i] - 1, i, :]
                           for i in range(d)])

    return B

def cheby_eval(x, n):
    past_val = np.ones_like(x)
    curr_val = x

    if n == 1:
        return past_val

    if n == 2:
        return curr_val

    for i in xrange(3, n + 1):
        temp = 2*x*curr_val - past_val
        past_val = curr_val
        curr_val = temp

    return curr_val

def cheby2n(x, n):
    # computes the chebychev polynomials of the first kind
    dim = x.shape
    results = np.zeros((n+1, ) + dim)
    results[0, ...] = np.ones(dim)
    results[1, ...] = x
    for i in range(2,n+1):
        results[i,...] = 2 * x * results[i-1,...] - results[i-2,...]
    return results


def num_grid_points(d, mu):
    if mu == 1:
        return 2*d + 1

    if mu == 2:
        return 1 + 4*d + 4*d*(d-1)/2.

    if mu == 3:
        return 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6.

def m_i(i):
    if i < 0:
        raise ValueError('i must be positive')
    elif i == 0:
        return 0
    elif i == 1:
        return 1
    else:
        return 2**(i - 1) + 1


class SmolyakGrid(object):
    """
    This class currently takes a dimension and a degree of polynomial
    and builds the Smolyak Sparse grid.  We base this on the work by
    Judd Maliar Maliar and Valero.  Hope to obtain speed ups beyond
    what they achieved.

    Parameters
    ----------
    d : scalar : integer
        This is the dimension of grid that you are building

    mu : scalar : integer
        mu is a parameter that defines the fineness of grid that we
        want to build


    Attributes
    ----------
    d : scalar : integer
        This is the dimension of grid that you are building

    mu : scalar : integer
        mu is a parameter that defines the fineness of grid that we
        want to build

    grid : array, float, ndim=2
        This is the sparse grid that we need to build

    inds : list : list
        This is a lists of lists that contains all of the indices

    B : np.ndarray : floats
        This is the B matrix that is used to do lagrange interpolation

    B_L : np.ndarray : floats
        Lower triangle matrix of the decomposition of B

    B_U : np.ndarray : floats
        Upper triangle matrix of the decomposition of B

    Examples
    --------
    >>> s = SmolyakGrid(3, 2)
    >>> s
    Smolyak Grid:
        d: 3
        mu: 2
        npoints: 25

    """

    def __init__(self, d, mu):
        self.d = d
        self.mu = mu

        if mu < 1:
            raise ValueError('The parameter mu needs to be > 1.')

        if d <= 1:
            raise ValueError('You are trying to build a one dimensional\
                             grid.')

        self.grid, self.inds = build_grid(self.d, self.mu)
        # self.inds = smol_inds(self.d, self.mu)
        self.B = build_B(self.d, self.mu, self.grid)
        # Compute LU decomposition
        l, u = lu(self.B, True)  # pass permute_l as true. See scipy docs
        self.B_L = l
        self.B_U = u



    def __repr__(self):
        msg = "Smolyak Grid:\n\td: {0} \n\tmu: {1} \n\tnpoints: {2}"
        return msg.format(self.d, self.mu, self.grid.shape[0])

    def __str__(self):
        return str(self.__repr__)

    def plot_grid(self):
        """
        This is a simple function that beautifully plots the grid for
        the 2 and 3d cases

        """
        grid = self.grid
        if grid.shape[1] == 2:
            xs = grid[:, 0]
            ys = grid[:, 1]
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(xs, ys)
            ax.grid(True, linestyle='--', color='0.75')
            plt.show()
        elif grid.shape[1] == 3:
            xs = grid[:, 0]
            ys = grid[:, 1]
            zs = grid[:, 2]
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(xs, ys, zs)
            ax.grid(True, linestyle='--', color='0.75')
            plt.show()
        else:
            raise ValueError('Can only plot 2 or 3 dimensional problems')
