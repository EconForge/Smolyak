"""
Uses the grid type from grid.jl to construct an interpolation type

"""

type SmolyakInterp
    SGrid::SmolyakGrid
    f_on_grid::Vector{Float64}
    theta::Vector{Float64}

    function SmolyakInterp(sg::SmolyakGrid, fgrid::Vector{Float64})
        # Constructor with a grid and the value of function on the grid

        theta = find_theta(sg, fgrid)::Vector{Float64}
        return new(sg, fgrid, theta)
    end

    # Default Constructor (just in case)
    function SmolyakInterp(sg, f, theta)
        new(sg, f, theta)
    end
end


function find_theta{T <: Real}(sg::SmolyakGrid, f_on_grid::Vector{T})
    # Use pre-computed inverse of B to get coefficients
   return sg.B_fact \ f_on_grid
end


# Convenience method to update coefficients.
function update_theta{T <: Real}(si::SmolyakInterp, f_on_grid::Vector{T})
    si.theta = find_theta(si.SGrid, f_on_grid)
end


# Interpolation
function interpolate(si::SmolyakInterp, pts::Matrix{Float64})
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
    # We should do some type checking
    n = size(pts, 1)
    d = size(pts, 2)

    @assert d == si.SGrid.d # check number of dimensions

    cube_pts = dom2cube(pts, si.SGrid.lb, si.SGrid.ub)::Matrix{Float64}
    new_B = build_B(d, si.SGrid.mu, cube_pts, si.SGrid.pinds)::Matrix{Float64}
    ipts = (new_B * si.theta)::Vector{Float64}
    return ipts
end


# function interpolate(si::SmolyakInterp, pts::Matrix{Float64}, deriv::Bool=false)
#     if deriv
#         n = size(pts, 1)
#         d = size(pts, 2)
#         @assert d == si.SGrid.d # check number of dimensions

#         Us = cheby2n(pts, maximum(si.SGrid.mu) + 1, kind=2.0)
#         Us = cat(3, zeros(n, d, 1), Us)
#         # for i = 1:length(si.theta)

#     else
#         return interpolate(si, pts)
#     end

# end


