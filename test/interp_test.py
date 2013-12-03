# TESTS: Leaving these here for now so that I don't have to do anything
# fancy
from smolyak.grid import SmolyakGrid
from smolyak.interp import SmolyakInterp
import numpy as np

# func = lambda x, y: np.exp(x**2 - y**2)
func = lambda x, y: x**2 - y**2


func1 = lambda points: func(points[:, 0], points[:, 1])
func1_prime = lambda x: np.column_stack([2 * x[:, 0], -2*x[:, 1]])

func2 = lambda x: np.sum(x ** 2, axis=1)
func2_prime = lambda x: 2 * x

def test_interp_2d(d, mu, f):
    sg = SmolyakGrid(d, mu, np.array([-2, -2.]), np.array([2., 2.]))

    f_on_grid = f(sg.org_grid)

    si = SmolyakInterp(sg, f_on_grid)

    np.random.seed(42)
    test_points = np.random.randn(100, 2)
    # Make sure it is bounded by -2, 2
    test_points = 2*test_points/np.max(np.abs(test_points))

    true_vals = f(test_points)
    interp_vals = si.interpolate(test_points)

    mean_ad = np.mean(np.abs(interp_vals - true_vals))
    max_ad = np.max(np.abs(interp_vals - true_vals))
    min_ad = np.min(np.abs(interp_vals - true_vals))

    print("The mean abs diff is {}\nThe max abs diff is {}\nThe min abs diff is {}"
          .format(mean_ad, max_ad, min_ad))
    return

def test_interp_2d1(d, mu, f):
    sg = SmolyakGrid(d, mu, np.array([-1, -1.]), np.array([1., 1.]))

    f_on_grid = f(sg.org_grid)

    si = SmolyakInterp(sg, f_on_grid)

    np.random.seed(42)
    test_points = np.random.randn(100, 2)
    # Make sure it is bounded by -2, 2
    test_points = test_points/np.max(np.abs(test_points))

    true_vals = f(test_points)
    interp_vals = si.interpolate(test_points)

    mean_ad = np.mean(np.abs(interp_vals - true_vals))
    max_ad = np.max(np.abs(interp_vals - true_vals))
    min_ad = np.min(np.abs(interp_vals - true_vals))

    print("The mean abs diff is {}\nThe max abs diff is {}\nThe min abs diff is {}"
          .format(mean_ad, max_ad, min_ad))
    return


def test_interp2d_dolo(d, mu, f):
    from dolo.numeric.interpolation.smolyak import SmolyakGrid as SG_dolo

    bds = np.ones(d)

    sg = SG_dolo(bds * -2, bds * 2, mu + 1)

    f_on_grid = f(sg.grid.T)
    sg.set_values(f_on_grid)

    np.random.seed(42)
    test_points = np.random.randn(100, 2)
    # Make sure it is bounded by -2, 2
    test_points = 2*test_points/np.max(np.abs(test_points))

    true_vals = f(test_points)
    interp_vals = sg.interpolate(test_points.T)

    mean_ad = np.mean(np.abs(interp_vals - true_vals))
    max_ad = np.max(np.abs(interp_vals - true_vals))
    min_ad = np.min(np.abs(interp_vals - true_vals))

    print("The mean abs diff is {}\nThe max abs diff is {}\nThe min abs diff is {}"
          .format(mean_ad, max_ad, min_ad))
    return


def test_interp2d_derivs_dolo(d, mu, f, f_prime):
    from dolo.numeric.interpolation.smolyak import SmolyakGrid as SG_dolo
    bds = np.ones(d)

    sg = SG_dolo(bds * -2, bds * 2, mu + 1)

    f_on_grid = f(sg.grid.T)
    sg.set_values(f_on_grid)

    np.random.seed(42)
    test_points = np.random.randn(100, 2)
    # Make sure it is bounded by -2, 2
    test_points = 2*test_points/np.max(np.abs(test_points))

    true_vals = f(test_points)
    true_vals_prime = f_prime(test_points)
    i_vals, i_vals_prime = sg.interpolate(test_points.T, with_derivative=True)
    i_vals_prime = np.squeeze(i_vals_prime).T

    mean_ad = np.mean(np.abs(i_vals - true_vals))
    max_ad = np.max(np.abs(i_vals - true_vals))
    min_ad = np.min(np.abs(i_vals - true_vals))

    mean_ad_prime = np.mean(np.abs(i_vals_prime - true_vals_prime))
    max_ad_prime = np.max(np.abs(i_vals_prime - true_vals_prime))
    min_ad_prime = np.min(np.abs(i_vals_prime - true_vals_prime))

    msg = "mean abs diff is {}\nmax abs diff is {}\nmin abs diff is {}"

    print("Interpolation results\n" + "#" * 21)
    print(msg.format(mean_ad, max_ad, min_ad))

    print("Derivative results\n" + "#" * 18)
    print(msg.format(mean_ad_prime, max_ad_prime, min_ad_prime))
    return true_vals_prime, i_vals_prime


# test_interp_2d(2, 3, func1)
# test_interp2d_dolo(2, 3, func1)
# test_interp_2d1(2, 3, func1)
