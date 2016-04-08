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
function interpolate!(out::AbstractVector, new_B::AbstractMatrix,
                      si::SmolyakInterp, pts::AbstractMatrix)
    d = size(pts, 2)

    d != si.grid.d && error("pts should have $d columns")

    cube_pts = dom2cube(pts, si.grid.lb, si.grid.ub)
    build_B!(new_B, d, si.grid.mu, cube_pts, si.grid.pinds)
    A_mul_B!(out, new_B, si.coefs)
end

function interpolate!(out::AbstractVector, si::SmolyakInterp, pts::AbstractMatrix)
    new_B = ones(size(pts, 1), size(si.grid.pinds, 1))
    interpolate!(out, new_B, si, pts)
end

interpolate(si::SmolyakInterp, pts::AbstractMatrix) =
    interpolate!(Array(Float64, size(pts, 1)), si, pts)

Base.call(si::SmolyakInterp, pts::AbstractMatrix) = interpolate(si, pts)
