using Distributions

require("../src/grid.jl")
require("../src/interp.jl")

# Define test function
func(x::Vector{Float64}, y::Vector{Float64}) = exp(- x.*x - y.*y)
func(g::Matrix{Float64}) = func(g[:, 1], g[:, 2])

# Construct Grid
sg = SmolyakGrid(2, 8)

# Evaluate the function on the grid and get coefficient vector
f_on_grid = func(sg.grid)

si = SmolyakInterp(sg, f_on_grid)

# Generate some random normal data
mvtnorm = MvNormal([0., 0.], eye(2).*0.05)
s = rand(mvtnorm, 50)'  # Transpose it so it has same dims as sg.grid

true_vals = func(s)
interp_vals = interpolate(si, s)

mean_abs_diff = mean(abs(interp_vals - true_vals))
println("mean_abs_diff = $mean_abs_diff")
