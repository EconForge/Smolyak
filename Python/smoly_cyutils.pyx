import cython
import numpy as np
cimport cython
cimport numpy as np


# cdef find_A_n(int n):
#     """
#     This method finds all of the unidimensional disjoint sets
#     that we will use to construct the grid.  It improves on
#     past algorithms by noting that A_{n} = S_{n}[evens] except for
#     A_1 = {0} and A_2 = {-1, 1}. Additionally, A_{n} = A_{n+1}[odds]
#     This prevents the calculation of these nodes repeatedly.  Thus
#     we only need to calculate biggest of the S_n's to build the
#     sequence of A_n's
#     """

#     # # Start w finding the biggest Sn(We will subsequently reduce it)
#     cdef double [::1] Sn
#     Sn = find_S_n(n)

#     cdef dict A_chain
#     A_chain = {}
#     A_chain[1] = [0.]
#     A_chain[2] = [-1., 1.]

#     # Need a for loop to extract remaining elements
#     cdef int seq, num
#     cdef tuple temp

#     for seq in range(n, 2, -1):
#         num = Sn.size
#         # Need odd indices in python because indexing starts at 0
#         temp = tuple(Sn[range(1, num, 2)])
#         A_chain[seq] = temp
#         # A_chain.append(list(Sn[range(1, num, 2)]))
#         Sn = Sn[range(0, num, 2)]

#     return A_chain


# def cyfind_S_n(int n):
#     """
#     This method finds the element S_n for the Chebyshev Extrema
#     """

#     if n==1:
#         cdef double [::1] ret_val
#         ret_val = np.array([0.])
#         return ret_val

#     cdef int m_i
#     cdef int [::1] comp_vals
#     cdef double [::1] vals
#     # Apply the necessary transformation to get the nested sequence
#     m_i = 2**(n-1) + 1

#     # Create an array of values that will be passed in to calculate
#     # the set of values
#     comp_vals = np.arange(1., m_i + 1.)

#     # Values are - cos(pi(j-1)/(n-1)) for j in [1, 2, ..., n]
#     vals = -1. * np.cos(np.pi*(comp_vals - 1.)/(m_i-1.))
#     # vals[np.where(np.abs(vals) < 1e-14)] = 0.0
#     return vals



cdef gen_points2(long[:, ::1] indices, dict An):
    from itertools import product
    cdef list points
    cdef long[::1] el
    cdef int i

    points = []
    for el in indices:
        points.extend(list(product(*[An[i] for i in el])))

    return points


def build_sparse_grid(int d, int mu):
    """
    This method builds a grid for the object
    """
    # Get An chain
    cdef double [:, ::1] An
    An = find_A_n(mu + 1)

    # Need to capture up to value mu + 1 so in python need mu+2
    cdef list possible_values
    possible_values = range(1, mu + 2)

    # find all (i1, i2, ... id) such that their sum is in range
    # we want; this will cut down on later iterations
    poss_inds = [el for el in combinations_with_replacement(possible_values, d) \
                  if d<sum(el)<=d+mu]

    true_inds = [[el for el in permute(list(val))][1:] for val in poss_inds]

    # Add the d dimension 1 array so that we don't repeat it a bunch
    # of times
    true_inds.extend([[1]*d])
    tinds = np.vstack(true_inds)

    points = []

    for el in tinds:
        temp = [An[i] for i in el]
        # Save these indices that we iterate through because
        # we need them for the chebyshev polynomial combination
        # inds.append(el)
        points.extend(list(product(*temp)))

    grid = np.array(points)
    self.grid = grid

    return grid
