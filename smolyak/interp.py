"""
This file contains the interpolation routines for the grids that are
built using the grid.py file in the smolyak package...  Write more doc
soon.
"""

from __future__ import division
import numpy.linalg as la
from grid import SmolyakGrid


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
        self.f_on_grid = self.trans_vals(f_on_grid)
        self.theta = find_theta(sg, self.f_on_grid)

    def trans_vals(self, pts):
        """
        Takes the y-values that are being passed in and moves them into
        the correct domain of [-1, 1].  Note that pts should be 1 dim
        """
        low_bounds = np.min(pts)
        up_bounds = np.max(pts)

        self.centers = low_bounds + (up_bounds - low_bounds)/2
        self.radii = (up_bounds - low_bounds)/2

        trans_pts = (pts-self.centers)/self.radii

        return trans_pts

    def inv_trans_vals(self, pts):
        """
        Takes transformed y-values and moves them back into the original
        domain for the interpolation values.
        """
        inv_trans_vals = pts*self.radii + self.centers

        return inv_trans_vals

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
        trans_points = self.sg.trans_points(pts)

        new_B = build_B(d, sg.mu, trans_points, sg.pinds)
        vals = new_B.dot(theta)
        vals = self.inv_trans_vals(vals)
        return vals

