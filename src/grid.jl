"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebyshev
polynomials to interpolate and approximate the data.

Method based on Judd, Maliar, Maliar, Valero 2013 (W.P)

Date: Fri Oct 18 07:47:11 2013 -0700

Author: Chase Coleman and Spencer Lyon

"""

using Iterators
using Cartesian
import PyPlot
require("util")
import Base.show

## --------------- ##
#- Building Blocks -#
## --------------- ##


function num_grid_pts(d::Int, mu::Int)
    """
    Uses known polynomials for generating the number of grid points for
    a given d and mu. As of right now, this only works when mu <= 3.
    """

    if mu == 1
        return int(2d - 1)
    elseif mu == 2
        return int(1 + 4d + 4d*(d-1)/2)
    elseif mu == 3
        return int(1 + 8d + 12d*(d-1)/2 + 8d*(d-1)*(d-2)/6)
    else
        error("We only know the number of grid points for mu \\in [1, 2, 3]")
    end
end


function m_i(i::Int)
    """
    Compute m(i) for a given i. This is defined m_i(1) = 1,
    m_i(i) = 2^(i-1) + 1 for all i > 1.
    """
    if i < 0
        error("DomainError: i must be positive")
    elseif i == 0
        return 0
    elseif i == 1
        return 1
    else
        return 2^(i - 1) + 1
    end
end


function cheby2n{T<:Number}(x::Array{T, 1}, n::Int; kind::Real=1)
    """
    Evaluates the first n+1 Chebyshev polynomials of the 'kind' kind at x
    NOTE: This will only work when kind = 1 or kind = 2
    NOTE: We evaluate the first n+1, with the first dimension being all 1's
    """

    dim = length(x)
    results = Array(T, dim, n + 1)
    results[:, 1] = 1.
    results[:, 2] = kind * x
    for i=3:n+1
        results[:, i] = 2x .* results[:, i - 1] - results[:, i - 2]
    end
    return results
end


function cheby2n{T<:Number}(x::Array{T, 2}, n::Int; kind::Real=1)
    # Evaluates the first n+1 Chebyshev polynomials of the 'kind' kind at x
    # NOTE: This will only work when kind = 1 or kind = 2
    # NOTE: We evaluate the first n+1, with the first dimension being all 1's

    dim = size(x)
    results = Array(T, dim..., n + 1)
    results[:, :, 1] = 1.
    results[:, :, 2] = kind * x
    for i=3:n+1
        results[:, :, i] = 2x .* results[:, :, i - 1] - results[:, :, i - 2]
    end
    return results
end


function s_n(n::Int)
    """
    Compute the set s_n. All points are extrema of Chebyshev polynomials
    """

    if n < 1
        error("DomainError: n must be positive")
    end

    if n == 1
        return [0.]
    end

    m = m_i(n)
    j = [1:m]
    pts = cos(pi * (j - 1.) / (m - 1))::Array{Float64, 1}
    for i = 1:size(pts, 1)
        pts[i] = abs(pts[i]) < 1e-12 ? 0.0 : pts[i]
    end
    return pts
end


function a_n(n::Int)
    """
    Compute the nth disjoint set of smolyak points.
    """
    if n < 1
        error("DomainError: n must be positive")
    end

    if n == 1
        return [0.]
    elseif n == 2
        return [-1., 1.]
    end

    sn = s_n(n)
    return sn[[2:2:size(sn, 1)]]
end


function a_chain(n::Int)
    """
    Return Dict of all A_n from n to 1. Keys are n, values are arrays of
    points in A_n

    This function improves on other algorithms by noting that A_{n} =
    S_{n}[evens] except for A_1 = {0} and A_2 = {-1, 1}. Additionally,
    A_{n} = A_{n+1}[odds] This prevents the calculation of these nodes
    repeatedly. Thus we only need to calculate biggest of the S_n's to
    build the sequence of A_n's

    """

    sn = s_n(n)

    a = Dict{Int, Array{Float64, 1}}()
    sizehint(a, n)

    # These are constant and don't follow the pattern.
    a[1] = [0.]
    a[2] = [-1., 1.]

    for i=n:-1:3
        # push!(a, sn[[2:2:size(sn, 1)]])
        a[i] = sn[[2:2:size(sn, 1)]]
        sn = sn[[1:2:size(sn, 1)]]
    end

    return a
end


function phi_chain(n::Int)
    max_ind = m_i(n)
    phi = Dict{Int, Array{Int64, 1}}()
    phi[1] = [1]
    phi[2] = [2, 3]
    low_ind = 4  # current lower index
    for i = 3:n
        high_ind = m_i(i)
        phi[i] = [low_ind:high_ind]
        low_ind = high_ind + 1
    end

    return phi
end


## ---------------------- ##
#- Construction Utilities -#
## ---------------------- ##

# Isotropic inds
function smol_inds(d::Int, mu::Int)
    """
    This function finds all of the indices that satisfy the requirement
    that d \leq \sum_{i=1}^d \leq d + \mu.

    Parameters
    ----------
    d : int
        The number of dimensions in the grid

    mu : int
        The parameter mu defining the density of the grid

    Returns
    -------
    true_inds : Array
        A 1-d Any array containing all d element arrays satisfying the
        constraint

    Notes
    -----
    This function is used directly by build_grid and poly_inds

    """

    p_vals = [1:mu + 1]

    poss_inds = {}

    for el in Task(()-> comb_with_replacement(p_vals, d))
        if d < sum(el) <= d + mu
            push!(poss_inds, el)
        end
    end

    true_inds = Any[ones(Int64, d)]  # we will always have (1, 1, ...,  1)
    for val in poss_inds
        for el in Task(()->pmute(val))
            push!(true_inds, el)
        end
    end

    return true_inds
end


# Ansotropic inds
function smol_inds(d::Int, mu::Array{Int, 1})
    # Compute indices needed for anisotropic smolyak grid given number of
    # dimensions d and a vector of mu parameters mu

    if length(mu) != d
        error("ValueError: mu must have d elements.")
    end

    mu_max = maximum(mu)
    mup1 = mu + 1

    p_vals = [1:mu_max + 1]

    poss_inds = {}

    for el in Task(()-> comb_with_replacement(p_vals, d))
        if d < sum(el) <= d + mu_max
            push!(poss_inds, el)
        end
    end

    true_inds = Any[ones(Int64, d)]
    for val in poss_inds
        for el in Task(()->pmute(val))
            if all(el .<= mup1)
                push!(true_inds, el)
            end
        end
    end

    return true_inds
end


function poly_inds(d::Int, mu::Int)
    """
    Build indices specifying all the cartesian products of Chebychev
    polynomials needed to build smolyak polynomial

    Parameters
    ----------
    d : Int
        The number of dimensions in grid / polynomial

    mu : int
        The parameter mu defining the density of the grid

    Returns
    -------
    phi_inds : Array{Int, 2}
        A two dimensional array of integers where each row specifies a
        new set of indices needed to define a smolyak basis polynomial

    Notes
    -----
    This function uses smol_inds and phi_chain. The output of this
    function is used by build_B to construct the B matrix

    """

    true_inds = smol_inds(d, mu)
    phi_n = phi_chain(mu + 1)
    return vcat([ cartprod([phi_n[i] for i in el]) for el in true_inds]...)
end


# Build grid
function build_grid(d::Int, mu::Union(Int, Vector{Int}))
    """
    Use disjoint Smolyak sets to construct Smolyak grid of degree d and
    density parameter mu

    The return value is an n x d Array, where n is the number of points
    in the grid

    """

    true_inds = smol_inds(d, mu)
    An = a_chain(maximum(mu) + 1)
    return vcat([ cartprod([An[i] for i in el]) for el in true_inds]...)
end


# Build B-matrix
function build_B(d::Int, mu::Int, grid::Array{Float64, 2})
    """
    Compute the matrix B from equation 22 in JMMV 2013
    Translation of dolo.numeric.interpolation.smolyak.SmolyakBasic

    Parameters
    ----------
    d : Int
        The number of dimensions on the grid

    mu : Int
        The mu parameter used to define grid

    grid : Array{Float64, 2}
        The smolyak grid returned by calling `build_grid(d, mu)`

    Returns
    -------
    B : Array{Float64, 2}
        The matrix B that represents the Smolyak polynomial
        corresponding to grid

    """

    Ts = cheby2n(grid, m_i(mu + 1))
    n = size(grid, 1)
    b_inds = poly_inds(d, mu)
    B = Array(Float64, n, n)
    for ind = 1:n
        B[:, ind] = reduce(.*, {slice(Ts, :, k, b_inds[ind, k]) for k=1:d})
    end
    return B
end


## ----------------- ##
#- Type: SmolyakGrid -#
## ----------------- ##


### type SmolyakGrid
type SmolyakGrid
    d::Int  # number of dimensions
    mu::Union(Int, Vector{Int})  # density. Int or d element vector of Int
    grid::Matrix{Float64}  # Smolyak grid
    inds::Array{Any, 1}  # Smolyak indices
    B::Matrix{Float64}  # matrix representing interpoland
    B_L::Matrix{Float64}  # L from LU decomposition of B
    B_U::Matrix{Float64}  # U from LU decomposition of B

    # Isotropic constructor
    function SmolyakGrid(d::Int, mu::Int)
        if d < 2
            error("You passed d = $d. d must be greater than 1")
        end
        if mu < 1
            error("You passed mu = $mu. mu must be greater than 1")
        end

        grid = build_grid(d, mu)
        inds = smol_inds(d, mu)
        B = build_B(d, mu, grid)
        lu_B = lu(B)
        B_L = lu_B[1]
        B_U = lu_B[2]

        new(d, mu, grid, inds, B, B_L, B_U)
    end

    # Anisotropic constructor
    function SmolyakGrid(d::Int, mu::Array{Int, 1})
        if d < 2
            error("You passed d = $d. d must be greater than 1")
        end
        if any(mu .< 1)
            error("You passed mu = $mu. each element must be greater than 1")
        end

        if length(mu) != d
            error("mu must have d elements. It has $(length(mu))")
        end

        grid = build_grid(d, mu)
        inds = smol_inds(d, mu)
        B = build_B(d, maximum(mu), grid)
        lu_B = lu(B)
        B_L = lu_B[1]
        B_U = lu_B[2]

        new(d, mu, grid, inds, B, B_L, B_U)
    end

    # Default constructor, just in case someone gets brave
    function SmolyakGrid(d, mu, grid, inds, B, B_L, B_U)
        new(d, mu, grid, inds, B, B_L, B_U)
    end
end


function show(io::IO, sg::SmolyakGrid)
    "show method for SmolyakGrid type"
    npoints = size(sg.grid, 1)
    non_zero_pts = nnz(sg.B)
    pct_non_zero = non_zero_pts / (npoints ^ 2)
    if isa(sg.mu, Array)
        mu_str = replace(strip(string(sg.mu)), '\n', ", ")
        msg = "Anisotropic Smolyak Grid:\n"
        msg *= "\td: $(sg.d)\n\tmu: $(mu_str)\n\tnpoints: $npoints"
        msg *= @sprintf "\n\tB: %.2f%% non-zero" pct_non_zero
    else
        msg = "Smolyak Grid:\n\td: $(sg.d)\n\tmu: $(sg.mu)\n\tnpoints: $npoints"
        msg *= @sprintf "\n\tB: %.2f%% non-zero" pct_non_zero
    end
    print(io, msg)
end


## ---------------- ##
#- Utilty funcitons -#
## ---------------- ##

function plot(sg::SmolyakGrid)
    g = sg.grid
    if size(g, 2) == 2
        PyPlot.plot(g[:, 1], g[:, 2], "b.")
        PyPlot.xlim((-1.5, 1.5))
        PyPlot.ylim((-1.5, 1.5))
        PyPlot.title("Smolyak grid: \$d=$(sg.d), \\; \\mu=$(sg.mu)\$")
    elseif size(g, 2) == 3
        PyPlot.scatter3D(g[:, 1], g[:, 2], g[:, 3], "b.")
        PyPlot.title("Smolyak grid: \$d=$(sg.d), \\; \\mu=$(sg.mu)\$")
    else
        error("ERROR: can only plot 2d or 3d grids")
    end
end

#----------
# Test case from dolo
#----------

# fun(x, y) = exp(- x .^ 2 - y .^ 2)
