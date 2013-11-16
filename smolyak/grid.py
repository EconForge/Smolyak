"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebychev
polynomials to interpolate and approximate the data.

Method based on Judd, Maliar, Maliar, Valero 2013 (W.P)

Authors
=======

- Chase Coleman (ccoleman@stern.nyu.edu)
- Spencer Lyon (slyon@stern.nyu.edu)

References
==========
Judd, Kenneth L, Lilia Maliar, Serguei Maliar, and Rafael Valero. 2013.
    "Smolyak Method for Solving Dynamic Economic Models: Lagrange
    Interpolation, Anisotropic Grid and Adaptive Domain".

Krueger, Dirk, and Felix Kubler. 2004. "Computing Equilibrium in OLG
    Models with Stochastic Production." Journal of Economic Dynamics and
    Control 28 (7) (April): 1411-1436.

"""
from __future__ import division
from operator import mul
from itertools import product, combinations_with_replacement
from itertools import chain
import numpy as np
from scipy.linalg import lu
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from util import *

## --------------- ##
#- Building Blocks -#
## --------------- ##


def num_grid_points(d, mu):
    """
    Checks the number of grid points for a given d, mu combination.

    Parameters
    ----------
    d, mu : int
        The parameters d and mu that specify the grid

    Returns
    -------
    num : int
        The number of points that would be in a grid with params d, mu

    Notes
    -----
    This function is only defined for mu = 1, 2, or 3

    """
    if mu == 1:
        return 2*d + 1

    if mu == 2:
        return 1 + 4*d + 4*d*(d-1)/2.

    if mu == 3:
        return 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6.


def m_i(i):
    """
    Compute one plus the "total degree of the interpolating
    polynoimals" (Kruger & Kubler, 2004). This shows up many times in
    Smolyak's algorithm. It is defined as:

    .. math::

        m_i = \\begin\{cases\}
        1 \\quad & \\text\{if \} i = 1 \\\\
        2^\{i-1\} + 1 \\quad & \\text\{if \} i \\geq 2
        \\end\{cases\}

    Parameters
    ----------
    i : int
        The integer i which the total degree should be evaluated

    Returns
    -------
    num : int
        Return the value given by the expression above

    """
    if i < 0:
        raise ValueError('i must be positive')
    elif i == 0:
        return 0
    elif i == 1:
        return 1
    else:
        return 2**(i - 1) + 1


def cheby2n(x, n, kind=1.):
    """
    Computes the first :math:`n+1` Chebychev polynomials of the first
    kind evaluated at each point in :math:`x` .

    Parameters
    ----------
    x : float or array(float)
        A single point (float) or an array of points where each
        polynomial should be evaluated

    n : int
        The integer specifying which Chebychev polynomial is the last
        to be computed

    kind : float, optional(default=1.0)
        The "kind" of Chebychev polynomial to compute. Only accepts
        values 1 for first kind or 2 for second kind

    Returns
    -------
    results : array (float, ndim=x.ndim+1)
        The results of computation. This will be an :math:`(n+1 \\times
        dim \\dots)` where :math:`(dim \\dots)` is the shape of x. Each
        slice along the first dimension represents a new Chebychev
        polynomial. This dimension has length :math:`n+1` because it
        includes :math:`\\phi_0` which is equal to 1 :math:`\\forall x`

    """
    x = np.asarray(x)
    dim = x.shape
    results = np.zeros((n + 1, ) + dim)
    results[0, ...] = np.ones(dim)
    results[1, ...] = x * kind
    for i in range(2, n+1):
        results[i, ...] = 2 * x * results[i-1, ...] - results[i-2, ...]
    return results


def s_n(n):
    """
    Finds the set :math:`S_n` , which is the :math:`n` th Smolyak set of
    Chebychev extrema

    Parameters
    ----------
    n : int
        The index :math:`n` specifying which Smolyak set to compute

    Returns
    -------
    s_n : array (float, ndim=1)
        An array containing all the Chebychev extrema in the set
        :math:`S_n`

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


def a_chain(n):
    """
    Finds all of the unidimensional disjoint sets of Chebychev extrema
    that are used to construct the grid.  It improves on past algorithms
    by noting  that :math:`A_{n} = S_{n}` [evens] except for :math:`A_1
    = \{0\}`  and :math:`A_2 = \{-1, 1\}` . Additionally, :math:`A_{n} =
    A_{n+1}` [odds] This prevents the calculation of these nodes
    repeatedly. Thus we only need to calculate biggest of the S_n's to
    build the sequence of :math:`A_n` 's

    Parameters
    ----------
    n : int
      This is the number of disjoint sets from Sn that this should make

    Returns
    -------
    A_chain : dict (int -> list)
      This is a dictionary of the disjoint sets that are made.  They are
      indexed by the integer corresponding

    """

    # # Start w finding the biggest Sn(We will subsequently reduce it)
    Sn = s_n(n)

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


def phi_chain(n):
    """
    For each number in 1 to n, compute the Smolyak indices for the
    corresponding basis functions. This is the :math:`n` in
    :math:`\\phi_n`

    Parameters
    ----------
    n : int
        The last Smolyak index :math:`n` for which the basis polynomial
        indices should be found

    Returns
    -------
    aphi_chain : dict (int -> list)
        A dictionary whose keys are the Smolyak index :math:`n` and
        values are lists containing all basis polynomial subscripts for
        that Smolyak index

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

## ---------------------- ##
#- Construction Utilities -#
## ---------------------- ##


def smol_inds(d, mu):
    """
    Finds all of the indices that satisfy the requirement that
    :math:`d \leq \sum_{i=1}^d \leq d + \mu`.

    Parameters
    ----------
    d : int
        The number of dimensions in the grid

    mu : int or array (int, ndim=1)
        The parameter mu defining the density of the grid. If an array,
        there must be d elements and an anisotropic grid is formed

    Returns
    -------
    true_inds : array
        A 1-d Any array containing all d element arrays satisfying the
        constraint

    Notes
    -----
    This function is used directly by build_grid and poly_inds

    """
    if isinstance(mu, int):
        max_mu = mu
    else:
        if mu.size != d:
            raise ValueError("mu must have d elements. It has %i" % mu.size)
        max_mu = int(np.max(mu))

    # Need to capture up to value mu + 1 so in python need mu+2
    possible_values = range(1, max_mu + 2)

    # find all (i1, i2, ... id) such that their sum is in range
    # we want; this will cut down on later iterations
    poss_inds = [el for el in combinations_with_replacement(possible_values, d)
                 if d < sum(el) <= d+max_mu]

    if isinstance(mu, int):
        true_inds = [[el for el in permute(list(val))] for val in poss_inds]
    else:
        true_inds = [[el for el in permute(list(val)) if all(el <= mu+1)]
                     for val in poss_inds]

    # Add the d dimension 1 array so that we don't repeat it a bunch
    # of times
    true_inds.extend([[[1]*d]])

    tinds = list(chain.from_iterable(true_inds))

    return tinds


def poly_inds(d, mu, inds=None):
    """
    Build indices specifying all the Cartesian products of Chebychev
    polynomials needed to build Smolyak polynomial

    Parameters
    ----------
    d : int
        The number of dimensions in grid / polynomial

    mu : int
        The parameter mu defining the density of the grid

    inds : list (list (int)), optional (default=None)
        The Smolyak indices for parameters d and mu. Should be computed
        by calling `smol_inds(d, mu)`. If None is given, the indices
        are computed using this function call

    Returns
    -------
    phi_inds : array : (int, ndim=2)
        A two dimensional array of integers where each row specifies a
        new set of indices needed to define a Smolyak basis polynomial

    Notes
    -----
    This function uses smol_inds and phi_chain. The output of this
    function is used by build_B to construct the B matrix

    """
    if inds is None:
        inds = smol_inds(d, mu)

    if isinstance(mu, int):
        max_mu = mu
    else:
        max_mu = max(mu)

    aphi = phi_chain(max_mu + 1)

    base_polys = []

    for el in inds:
        temp = [aphi[i] for i in el]
        # Save these indices that we iterate through because
        # we need them for the chebychev polynomial combination
        # inds.append(el)
        base_polys.extend(list(product(*temp)))

    return base_polys


def build_grid(d, mu, inds=None):
    """
    Use disjoint Smolyak sets to construct Smolyak grid of degree d and
    density parameter :math:`mu`

    The return value is an :math:`n \\times d` Array, where :math:`n`
    is the number of points in the grid

    Parameters
    ----------
    d : int
        The number of dimensions in the grid

    mu : int
        The density parameter for the grid

    inds : list (list (int)), optional (default=None)
        The Smolyak indices for parameters d and mu. Should be computed
        by calling `smol_inds(d, mu)`. If None is given, the indices
        are computed using this function call

    Returns
    -------
    grid : array (float, ndim=2)
        The Smolyak grid for the given d, :math:`mu`

    """
    if inds is None:
        inds = smol_inds(d, mu)

    # Get An chain
    if isinstance(mu, int):
        An = a_chain(mu + 1)
    else:  # Anisotropic case
        An = a_chain(max(mu) + 1)

    points = []

    # Need to get the correct indices

    for el in inds:
        temp = [An[i] for i in el]
        # Save these indices that we iterate through because
        # we need them for the chebychev polynomial combination
        # inds.append(el)
        points.extend(list(product(*temp)))

    grid = pd.lib.to_object_array_tuples(points).astype(float)

    return grid


def build_B(d, mu, grid, inds=None):
    """
    Compute the matrix B from equation 22 in JMMV 2013
    Translation of dolo.numeric.interpolation.smolyak.SmolyakBasic

    Parameters
    ----------
    d : int
        The number of dimensions on the grid

    mu : int or array (int, ndim=1, legnth=d)
        The mu parameter used to define grid

    grid : array (float, dims=2)
        The smolyak grid returned by calling `build_grid(d, mu)`

    inds : list (list (int)), optional (default=None)
        The Smolyak indices for parameters d and mu. Should be computed
        by calling `smol_inds(d, mu)`. If None is given, the indices
        are computed using this function call

    Returns
    -------
    B : array (float, 2)
        The matrix B that represents the Smolyak polynomial
        corresponding to grid

    """
    if inds is None:
        inds = smol_inds(d, mu)

    if isinstance(mu, int):
        max_mu = mu
    else:
        max_mu = max(mu)

    Ts = cheby2n(grid.T, m_i(max_mu + 1))
    base_polys = poly_inds(d, mu, inds)
    n = len(grid)
    B = np.empty((n, n), order='F')
    for ind, comb in enumerate(base_polys):
        B[:, ind] = reduce(mul, [Ts[comb[i] - 1, i, :]
                           for i in range(d)])

    return B


## ------------------ ##
#- Class: SmolyakGrid -#
## ------------------ ##

class SmolyakGrid(object):
    """
    This class currently takes a dimension and a degree of polynomial
    and builds the Smolyak Sparse grid.  We base this on the work by
    Judd, Maliar, Maliar, and Valero (2013).

    Parameters
    ----------
    d : int
        The number of dimensions in the grid

    mu : int or array(int, ndim=1, length=d)
        The "density" parameter for the grid


    Attributes
    ----------
    d : int
        This is the dimension of grid that you are building

    mu : int
        mu is a parameter that defines the fineness of grid that we
        want to build

    grid : array (float, ndim=2)
        This is the sparse grid that we need to build

    inds : list (list (int))
        This is a lists of lists that contains all of the indices

    B : array (float, ndim=2)
        This is the B matrix that is used to do lagrange interpolation

    B_L : array (float, ndim=2)
        Lower triangle matrix of the decomposition of B

    B_U : array (float, ndim=2)
        Upper triangle matrix of the decomposition of B

    Examples
    --------
    >>> s = SmolyakGrid(3, 2)
    >>> s
    Smolyak Grid:
        d: 3
        mu: 2
        npoints: 25
        B: 0.65% non-zero
    >>> ag = SmolyakGrid(3, [1, 2, 3])
    >>> ag
    Anisotropic Smolyak Grid:
        d: 3
        mu: 1 x 2 x 3
        npoints: 51
        B: 0.68% non-zero

    """
    def __init__(self, d, mu):
        self.d = d

        if d <= 1:
            raise ValueError('Number of dimensions must be >= 2')

        if isinstance(mu, int):  # Isotropic case
            if mu < 1:
                raise ValueError('The parameter mu needs to be > 1.')

            self.mu = mu
            self.inds = smol_inds(d, mu)
            self.grid = build_grid(self.d, self.mu, self.inds)
            self.B = build_B(self.d, self.mu, self.grid, self.inds)

        else:  # Anisotropic case
            mu = np.asarray(mu)

            if any(mu < 1):
                raise ValueError('Each element in mu needs to be > 1.')

            if len(mu) != d:
                raise ValueError("For Anisotropic grid, mu must have len d ")

            self.mu = mu
            self.inds = smol_inds(d, mu)
            self.grid = build_grid(self.d, self.mu, self.inds)
            self.B = build_B(self.d, self.mu, self.grid, self.inds)

        # Compute LU decomposition of B
        l, u = lu(self.B, True)  # pass permute_l as true. See scipy docs
        self.B_L = l
        self.B_U = u

    def __repr__(self):
        npoints = self.grid.shape[0]
        nz_pts = np.count_nonzero(self.B)
        pct_nz = nz_pts / (npoints ** 2.)

        if isinstance(self.mu, int):
            msg = "Smolyak Grid:\n\td: {0} \n\tmu: {1} \n\tnpoints: {2}"
            msg += "\n\tB: {3:.2f}% non-zero"
            return msg.format(self.d, self.mu, self.grid.shape[0], pct_nz)
        else:  # Anisotropic grid
            msg = "Anisotropic Smolyak Grid:"
            msg += "\n\td: {0} \n\tmu: {1} \n\tnpoints: {2}"
            msg += "\n\tB: {3:.2f}% non-zero"
            mu_str = " x ".join(map(str, self.mu))
            return msg.format(self.d, mu_str, self.grid.shape[0], pct_nz)

    def __str__(self):
        return self.__repr__()

    def plot_grid(self):
        """
        Beautifully plots the grid for the 2d and 3d cases

        Parameters
        ----------
        None

        Returns
        -------
        None

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
