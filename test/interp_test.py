# TESTS: Leaving these here for now so that I don't have to do anything
# fancy
from smolyak.grid import SmolyakGrid
from smolyak.interp import SmolyakInterp
import numpy as np

func = lambda x, y: np.exp(x**2 - y**2)

func1 = lambda points: func(points[:, 0], points[:, 1])

def test_interp_2d(d, mu, f):
    sg = SmolyakGrid(d, mu, np.array([-1, -1]), np.array([1, 1]))

    f_on_grid = f(sg.grid)

    si = SmolyakInterp(sg, f_on_grid)

    np.random.seed(42)
    test_points = np.random.randn(100, 2)
    # Make sure it is bounded by -1, 1
    test_points = test_points/np.max(np.abs(test_points))

    true_vals = f(test_points)
    interp_vals = si.interpolate(test_points)

    mean_ad = np.mean(np.abs(interp_vals - true_vals))
    max_ad = np.max(np.abs(interp_vals - true_vals))
    min_ad = np.min(np.abs(interp_vals - true_vals))

    print("The mean abs diff is {}\nThe max abs diff is {}\nThe min abs diff is {}"
          .format(mean_ad, max_ad, min_ad))
    return

test_interp_2d(2, 3, func1)
