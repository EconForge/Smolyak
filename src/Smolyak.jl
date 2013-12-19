module Smolyak

import Cartesian
import PyPlot
import Base.show

include("util.jl")

export
    SmolyakGrid,
    build_B,
    build_grid,
    plot
include("grid.jl")

export
    SmolyakInterp,
    interpolate,
    find_theta,
    update_theta
include("interp.jl")

end # module
