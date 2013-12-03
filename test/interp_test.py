# TESTS: Leaving these here for now so that I don't have to do anything
# fancy


def func(x, y):
    return np.exp(x**2 - y**2)

def func1(points):
    return func(points[:, 0], points[:, 1])

sg = SmolyakGrid(2, 3, np.array([-1, -1]), np.array([1, 1]))

f_on_grid = func1(sg.grid)

si = SmolyakInterp(sg, f_on_grid)

test_points = np.random.randn(100, 2)
# Make sure it is bounded by -1, 1
test_points = test_points/np.max(np.abs(test_points))

true_vals = func1(test_points)
interp_vals = si.interpolate(test_points)

mean_ad = np.mean(np.abs(interp_vals - true_vals))
max_ad = np.max(np.abs(interp_vals - true_vals))
min_ad = np.min(np.abs(interp_vals - true_vals))

print("The mean abs diff is {} \nThe max abs diff is {}\n The min abs diff is{}"\
      .format(mean_ad, max_ad, min_ad))
