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


function interp_withderiv(si::SmolyakInterp, pts::AbstractMatrix)
    cube_pts = dom2cube(pts, si.grid.lb, si.grid.ub)
    new_B, der_B = build_B_with_deriv(d, si.grid.mu, cube_pts, si.grid.pinds)
    vals = new_B * si.coefs
    deriv = zeros(size(pts, 1), si.grid.d)

    # TODO: once we are confident that this works, add @inbounds/@simd
    for i in 1:length(si.coefs)
        for I in CartesianRange(size(deriv))
            deriv[I] += der_B[I, i] * si.coefs[i]
        end
    end

    dom2cube!(deriv, deriv, si.grid.lb, si.grid.ub)
    vals, deriv
end
