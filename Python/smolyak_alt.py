"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebyshev
polynomials to interpolate and approximate the data.

Method based on Judd, Maliar, Maliar, Valero 2013 (W.P)

Author: Chase Coleman and Spencer Lyon
"""
from operator import mul
import sys
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations_with_replacement, permutations
from itertools import chain
import time as time
from dolo.numeric.interpolation import smolyak as smm
import pandas as pd


# class Smoly_grid(object):
#     """
#     This class currently takes a dimension and a degree of polynomial
#     and builds the Smolyak Sparse grid.  We base this on the work by
#     Judd Maliar Maliar and Valero.  Hope to obtain speed ups beyond
#     what they achieved.

#     Attributes:
#     ------------------
#     d : scalar : integer
#         This is the dimension of grid that you are building

#     mu : scalar : integer
#         mu is a parameter that defines the fineness of grid that we
#         want to build


#     Methods:
#     --------
#     grid_build : This method builds the sparse grid

#     """
#     def __init__(self, d, mu):
#         """
#         Parameters
#         ----------
#         d : scalar : integer
#             This is the dimension of grid that you are building

#         mu : scalar : integer
#             mu is a parameter that defines the fineness of grid that we
#             want to build
#         """
#         self.d = d
#         self.mu = mu

#         if mu<1:
#             raise ValueError('The parameter mu needs to be > 1.')

#         if d<=1:
#             raise ValueError('You are trying to build a one dimensional\
#                              grid.')

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

    # # Start w finding the biggest Sn(We will subsequently reduce it)
    Sn = _find_S_n(n)

    # Save to be used later in evaluations
    # self.Sn = Sn

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

def _find_S_n(n):
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

def _find_smol_inds(d, mu):
    """
    This method finds all of the indices that satisfy the requirement
    that d \leq \sum_{i=1}^d \leq d + \mu.  Once we have these, then
    they can be used to build both the grid and the polynomial
    basis.
    """

    # Need to capture up to value mu + 1 so in python need mu+2
    possible_values = range(1, mu + 2)

    # find all (i1, i2, ... id) such that their sum is in range
    # we want; this will cut down on later iterations
    poss_inds = [list(el) for el in combinations_with_replacement(possible_values, d) \
                  if d<sum(el)<=d+mu]

    true_inds = []

    true_inds = [[el for el in permute(list(val))] for val in poss_inds]


    # Add the d dimension 1 array so that we don't repeat it a bunch
    # of times
    true_inds.extend([[[1]*d]])

    tinds = list(chain.from_iterable(true_inds))

    return tinds, poss_inds

def build_sparse_grid(d, mu):
    """
    This method builds a grid for the object
    """

    # Get An chain
    An = _find_A_n(mu + 1)

    points = []

    # Need to get the correct indices

    tinds = _find_smol_inds(d, mu)[0]

    for el in tinds:
        temp = [An[i] for i in el]
        # Save these indices that we iterate through because
        # we need them for the chebyshev polynomial combination
        # inds.append(el)
        points.extend(list(product(*temp)))

    grid = pd.lib.to_object_array_tuples(points).astype(float)

    return grid

def calc_chebvals(n):
    """
    We will eventually need the chebyshev polynomials evaluated at
    every poiOnt that we calculate in the grid.  We will find them
    and store them in a dictionary of dictionaries.  The outer dict
    will say which chebyshev polynomial it is and the inner dict
    will have each point evaluated in it.

    TODO: Speed this function up.  Costly to build
    """

    pts = _find_S_n(n + 1)  # TODO: save this when we call _find_A_n
    max_poly = pts.size

    cheb_dict = {}
    for phi_n in xrange(1, max_poly + 1):
        cheb_dict[phi_n] = {pt: t_pt for (pt, t_pt) in
                            zip(pts, cheby_eval(pts, phi_n))}

    return cheb_dict, pts

def _find_aphi(n):
    """
    Finds the disjoint sets of aphi's that will be used to compute
    which functions we need to calculate
    """
    # Need to find which is biggest function that is going to be
    # called is
    max_ind = 2**(n-1) + 1


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

def build_basis_polynomials(self):
    """
    This function builds the base polynomials that will be used to
    interpolate.
    """
    d = self.d
    mu = self.mu

    smol_inds = self.smol_inds
    aphi = self._find_aphi(mu + 1)

    # Bring in polynomials
    cheb_dict = self.calc_chebvals()

    base_polys = []

    for el in smol_inds:
        temp = [aphi[i] for i in el]
        # Save these indices that we iterate through because
        # we need them for the chebyshev polynomial combination
        # inds.append(el)
        base_polys.extend(list(product(*temp)))

    self.base_polys = base_polys

    return base_polys

def _build_B(inds, basis, d, mu):

    cd, pts = calc_chebvals(mu)
    inds = _find_smol_inds(d, mu)[1]
    aphi = _find_aphi(mu + 1)

    bs = []

    for el in inds:
        # if sum(el) == len(el) * el[0]:
        temp = [aphi[i] for i in el]
        bs.extend(list(product(*temp)))






    for row in xrange(n):
        for col in xrange(n):
            pt = grid[row]
            f = basis[col]
            B[row, col] = reduce(mul, [cd[i][pt[k]] for k, i in enumerate(f)])

    return B



def plot_grid(self):
    import matplotlib.pyplot as plt
    grid = self.grid
    if grid.shape[1] == 2:
        xs = grid[:, 0]
        ys = grid[:, 1]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(xs, ys)
        ax.grid(True, linestyle='--',color='0.75')
        plt.show()
    elif grid.shape[1] == 3:
        xs = grid[:, 0]
        ys = grid[:, 1]
        zs = grid[:, 2]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xs, ys, zs)
        ax.grid(True, linestyle='--',color='0.75')
        plt.show()
    else:
        raise ValueError('Can only plot 2 or 3 dimensional problems')


def cheby_eval(x, n):
   past_val = np.ones_like(x)
   curr_val = x

   if n==1:
       return past_val

   if n==2:
       return curr_val

   for i in xrange(3, n+1):
       temp = 2*x*curr_val - past_val
       past_val = curr_val
       curr_val = temp

   return curr_val


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

                # yield list(a)
                return


def check_points(d, mu):
    if mu == 1:
        return 2*d + 1

    if mu == 2:
        return 1 + 4*d + 4*d*(d-1)/2.

    if mu == 3:
        return 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6.






