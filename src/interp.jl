type SmolyakInterp
    grid::SmolyakGrid
    f_on_grid::Vector{Float64}
    coefs::Vector{Float64}
end

SmolyakInterp(sg::SmolyakGrid, fgrid::AbstractVector) =
    SmolyakInterp(sg, fgrid, find_coefs(sg, fgrid))

# Use pre-computed LU fact of B to get coefficients
find_coefs(sg::SmolyakGrid, f_on_grid::AbstractVector) = sg.B_fact \ f_on_grid

# Convenience method to update coefficients.
update_coefs!(si::SmolyakInterp, f_on_grid) =
    (si.coefs = find_coefs(si.grid, f_on_grid))


# Interpolation
function interpolate(si::SmolyakInterp, pts::AbstractMatrix)
    # We should do some type checking
    n = size(pts, 1)
    d = size(pts, 2)

    @assert d == si.grid.d # check number of dimensions

    cube_pts = dom2cube(pts, si.grid.lb, si.grid.ub)
    new_B = build_B(d, si.grid.mu, cube_pts, si.grid.pinds)
    new_B * si.coefs
end

Base.call(si::SmolyakInterp, pts::AbstractMatrix) = interpolate(si, pts)
