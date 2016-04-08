module Smolyak

using Iterators

export SmolyakGrid, build_B, build_B!, build_grid,
       SmolyakInterp, interpolate, find_theta, update_theta

include("util.jl")
include("grid.jl")
include("interp.jl")

end # module
