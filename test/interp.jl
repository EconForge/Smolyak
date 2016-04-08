using Smolyak


const verbose = true
const print_str = "mean abs diff: %.3e\nmax abs diff: %.3e\nmin abs diff: %.3e"
@eval print_func(x, y, z) = @printf($print_str, x, y, z)

# L2-norm squared
f1(x) = squeeze(sum(x .^ 2, 2), 2)
f1_p(x) = 2 * x

f2(x) = sin(x[:,1]) + cos(x[:, 2])
f2_p(x) = [cos(x[:, 1]) -sin(x[:, 2])]


function test_interp2d_derivs(d::Int, mu::Int, f::Function, f_prime::Function, bds::Real=1)
    sg = SmolyakGrid(d, mu, -bds, bds)

    f_on_grid = f(sg.grid)
    si = SmolyakInterp(sg, f_on_grid)

    srand(42)
    test_points = randn(100, d)

    # Make sure it is bounded by -bds, bds
    test_points = (bds - 0.05) * test_points/maximum(abs(test_points))

    true_vals = f(test_points)
    true_vals_prime = f_prime(test_points)
    i_vals = si(test_points)

    mean_ad = mean(abs(i_vals - true_vals))
    max_ad = maximum(abs(i_vals - true_vals))
    min_ad = minimum(abs(i_vals - true_vals))

    if verbose
        println("Interpolation results\n" * "#" ^ 21)
        print_func(mean_ad, max_ad, min_ad)
    end

    # return i_vals_prime
end

test_interp2d_derivs(15, 3, f2, f2_p, -2)

@time test_interp2d_derivs(15, 3, f2, f2_p, -2)
