"""
This file contains the interpolation routines for the grids that are
built using the grid.py file in the smolyak package...  Write more doc
soon.
"""

from __future__ import division
from operator import mul
from itertools import product, combinations_with_replacement
from itertools import chain
import numpy as np
import numpy.linalg as la
from scipy.linalg import lu
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from util import *
from grid import *


def find_theta(sg, f_on_grid):
    """
    Given a SmolyakGrid object and the value of the function on the
    points of the grid, this function will return the coefficients theta
    """
    return la.solve(sg.B_U, la.solve(sg.B_L, f_on_grid))


class SmolyakInterp(object):
    """
    This class is going to take several inputs.  It will need a
    SmolyakGrid object to be passed in and the values of the function
    evaluated at the grid points
    """
    def __init__(self, sg, f_on_grid):
        self.sg = sg
        self.f_on_grid = f_on_grid
        self.theta = find_theta(sg, f_on_grid)

    def interpolate(self, pts, with_deriv=False):
        """
        Basic interpolation.

        Parameters
        ==========
        si : SmolyakInterp
            Instance of SmolyakInterp type

        pts : array (float, ndim=2)
            A 2d array of points on which to evaluate the function. Each
            row is assumed to be a new d-dimensional point. Therefore, pts
            must have the same number of columns as `si.SGrid.d`

        Returns
        =======
        ipts : array(float, ndim=1)
            The interpolated values for each point in `pts`

        Notes
        =====
        This is a stripped down port of `dolo.SmolyakBasic.interpolate`

        TODO: There may be a better way to do this

        """
        n = pts.shape[0]
        d = pts.shape[1]
        sg = self.sg

        # Give ourselves the value of theta
        theta = self.theta
        # Move points to correct domain
        trans_points = sg.trans_points(pts)

        new_B = build_B(d, sg.mu, trans_points, sg.pinds)
        vals = new_B.dot(theta)
        return vals

