"""
Uses the grid type from smolyak.jl to construct an interpolation type
"""
require("smolyak.jl")


type SmolyakInterp
    SGrid::SmolyakGrid
    u::Vector{Float64}
    l::Vector{Float64}
    f_on_grid::Matrix{Float64}
    theta::Vector{Float64}

    # Default Constructor (just in case)
    function SmolyakInterp(sg, u, l, f, theta)
        new(sg, u, l, f, theta)
    end

end


function find_theta{T :< Number}(sg::SmolyakGrd, f_on_grid::Matrix{T})
    # Use pre-computed LU decomp to get coefficients
   return sg.B_U \ (sg.B_L \ f_on_grid)
end


function SmolyakInterp(sg::SmolyakGrid,
                       u::Vector{Float64},
                       l::Vector{Float64},
                       f::Function)
    f_on_grid = f(sg.grid)
    theta = find_theta(sg, f_on_grid)
    SmolyakInterp(sg, u, l, )
end


# Convenience method to update coefficients.
function update_theta{T :< Number}(si::SmolyakInterp, f_on_grid::Matrix{T})
    si.theta = find_theta(si, f_on_grid)
end
