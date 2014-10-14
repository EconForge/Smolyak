module Smolyak

import Cartesian
# import PyPlot
import Base.show


export SmolyakGrid, build_B, build_grid, #=plot, =#  # grid.jl
       SmolyakInterp, interpolate, find_theta, update_theta  # interp.jl

include("util.jl")
include("grid.jl")
include("interp.jl")

end # module
