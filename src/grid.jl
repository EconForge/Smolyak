"""
This file contains a class that builds a Smolyak Grid.  The hope is that
it will eventually contain the interpolation routines necessary so that
the given some data, this class can build a grid and use the Chebychev
polynomials to interpolate and approximate the data.

Method based on Judd, Maliar, Maliar, Valero 2013 (W.P)

Authors
=======

- Chase Coleman (ccoleman@stern.nyu.edu)
- Spencer Lyon (slyon@stern.nyu.edu)

References
==========
Judd, Kenneth L, Lilia Maliar, Serguei Maliar, and Rafael Valero. 2013.
    "Smolyak Method for Solving Dynamic Economic Models: Lagrange
    Interpolation, Anisotropic Grid and Adaptive Domain".

Krueger, Dirk, and Felix Kubler. 2004. "Computing Equilibrium in OLG
    Models with Stochastic Production." Journal of Economic Dynamics and
    Control 28 (7) (April): 1411-1436.

"""
## --------------- ##
#- Building Blocks -#
## --------------- ##

# union type to handle isotropic/anisotropic mu parameters
IntSorV = Union(Int, Vector{Int})

function num_grid_pts(d::Int, mu::Int)
    """
    Checks the number of grid points for a given d, mu combination.

    Parameters
    ----------
    d, mu : int
        The parameters d and mu that specify the grid

    Returns
    -------
    num : int
        The number of points that would be in a grid with params d, mu

    Notes
    -----
    This function is only defined for mu = 1, 2, or 3

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
    Compute one plus the "total degree of the interpolating
    polynoimals" (Kruger & Kubler, 2004). This shows up many times in
    Smolyak's algorithm. It is defined as:

    .. math::

        m_i = \\begin\{cases\}
        1 \\quad & \\text\{if \} i = 1 \\\\
        2^\{i-1\} + 1 \\quad & \\text\{if \} i \\geq 2
        \\end\{cases\}

    Parameters
    ----------
    i : int
        The integer i which the total degree should be evaluated

    Returns
    -------
    num : int
        Return the value given by the expression above

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
    Computes the first :math:`n+1` Chebychev polynomials of the first
    kind evaluated at each point in :math:`x` .

    Parameters
    ----------
    x : float or array(float)
        A single point (float) or an array of points where each
        polynomial should be evaluated

    n : int
        The integer specifying which Chebychev polynomial is the last
        to be computed

    kind : float, optional(default=1.0)
        The "kind" of Chebychev polynomial to compute. Only accepts
        values 1 for first kind or 2 for second kind

    Returns
    -------
    results : array (float, ndim=x.ndim+1)
        The results of computation. This will be an :math:`(n+1 \\times
        dim \\dots)` where :math:`(dim \\dots)` is the shape of x. Each
        slice along the first dimension represents a new Chebychev
        polynomial. This dimension has length :math:`n+1` because it
        includes :math:`\\phi_0` which is equal to 1 :math:`\\forall x`

    """

    dim = length(x)::Int64
    results = Array(T, dim, n + 1)::Array{Float64, 2}
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
    Finds the set :math:`S_n` , which is the :math:`n` th Smolyak set of
    Chebychev extrema

    Parameters
    ----------
    n : int
        The index :math:`n` specifying which Smolyak set to compute

    Returns
    -------
    s_n : array (float, ndim=1)
        An array containing all the Chebychev extrema in the set
        :math:`S_n`

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
    Finds all of the unidimensional disjoint sets of Chebychev extrema
    that are used to construct the grid.  It improves on past algorithms
    by noting  that :math:`A_{n} = S_{n}` [evens] except for :math:`A_1
    = \{0\}`  and :math:`A_2 = \{-1, 1\}` . Additionally, :math:`A_{n} =
    A_{n+1}` [odds] This prevents the calculation of these nodes
    repeatedly. Thus we only need to calculate biggest of the S_n's to
    build the sequence of :math:`A_n` 's

    Parameters
    ----------
    n : int
      This is the number of disjoint sets from Sn that this should make

    Returns
    -------
    A_chain : dict (int -> list)
      This is a dictionary of the disjoint sets that are made.  They are
      indexed by the integer corresponding

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
    """
    For each number in 1 to n, compute the Smolyak indices for the
    corresponding basis functions. This is the :math:`n` in
    :math:`\\phi_n`

    Parameters
    ----------
    n : int
        The last Smolyak index :math:`n` for which the basis polynomial
        indices should be found

    Returns
    -------
    aphi_chain : dict (int -> list)
        A dictionary whose keys are the Smolyak index :math:`n` and
        values are lists containing all basis polynomial subscripts for
        that Smolyak index

    """
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
    Finds all of the indices that satisfy the requirement that
    :math:`d \leq \sum_{i=1}^d \leq d + \mu`.

    Parameters
    ----------
    d : int
        The number of dimensions in the grid

    mu : int or array (int, ndim=1)
        The parameter mu defining the density of the grid. If an array,
        there must be d elements and an anisotropic grid is formed

    Returns
    -------
    true_inds : array
        A 1-d Any array containing all d element arrays satisfying the
        constraint

    Notes
    -----
    This function is used directly by build_grid and poly_inds

    """

    p_vals = [1:mu + 1]

    # PERF: size_hint here if it is slow
    poss_inds = {}

    for el in Task(()-> comb_with_replacement(p_vals, d))
        if d < sum(el) <= d + mu
            push!(poss_inds, el)
        end
    end

    # PERF: size_hint here if it is slow
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

    # PERF: size_hint here if it is slow
    poss_inds = {}

    for el in Task(()-> comb_with_replacement(p_vals, d))
        if d < sum(el) <= d + mu_max
            push!(poss_inds, el)
        end
    end

    # PERF: size_hint here if it is slow
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


function poly_inds(d::Int, mu::IntSorV, inds::Array{Any}={})
    """
    Build indices specifying all the Cartesian products of Chebychev
    polynomials needed to build Smolyak polynomial

    Parameters
    ----------
    d : int
        The number of dimensions in grid / polynomial

    mu : int
        The parameter mu defining the density of the grid

    inds : list (list (int)), optional (default=None)
        The Smolyak indices for parameters d and mu. Should be computed
        by calling `smol_inds(d, mu)`. If None is given, the indices
        are computed using this function call

    Returns
    -------
    phi_inds : array : (int, ndim=2)
        A two dimensional array of integers where each row specifies a
        new set of indices needed to define a Smolyak basis polynomial

    Notes
    -----
    This function uses phi_chain. The output of this function is used by
    build_B to construct the B matrix

    """
    if length(inds) == 0
        inds = smol_inds(d, mu)
    end
    phi_n = phi_chain(maximum(mu) + 1)
    return vcat([ cartprod([phi_n[i] for i in el]) for el in inds]...)
end


# Build grid
function build_grid(d::Int, mu::IntSorV, inds::Array{Any}={})
    """
    Use disjoint Smolyak sets to construct Smolyak grid of degree d and
    density parameter :math:`mu`

    The return value is an :math:`n \\times d` Array, where :math:`n`
    is the number of points in the grid

    Parameters
    ----------
    d : int
        The number of dimensions in the grid

    mu : int
        The density parameter for the grid

    inds : list (list (int)), optional (default=None)
        The Smolyak indices for parameters d and mu. Should be computed
        by calling `smol_inds(d, mu)`. If None is given, the indices
        are computed using this function call

    Returns
    -------
    grid : array (float, ndim=2)
        The Smolyak grid for the given d, :math:`mu`

    """
    if length(inds) == 0
        inds = smol_inds(d, mu)
    end
    An = a_chain(maximum(mu) + 1)  # use maximum in case mu is Vector
    return vcat([ cartprod([An[i] for i in el]) for el in inds]...)
end


# Build B-matrix
function build_B(d::Int, mu::IntSorV, pts::Array{Float64, 2},
                 b_inds::Array{Int64, 2})
    """
    Given polynomial indices, construct the matrix B from equation
    22 in JMMV 2013

    Optimized translation of
    `dolo.numeric.interpolation.smolyak.SmolyakBasic`

    Parameters
    ----------
    d : int
        The number of dimensions on the grid

    mu : int or array (int, ndim=1, legnth=d)
        The mu parameter used to define grid

    pts : array (float, dims=2)
        Arbitrary d-dimensional points. Each row is assumed to be a new
        point. Often this is the smolyak grid returned by calling
        `build_grid(d, mu)`

    b_inds : array (int, ndim=2)
        The polynomial indices for parameters a given grid. These should
        be  computed by calling `poly_inds(d, mu)`.

    Returns
    -------
    B : array (float, 2)
        If `pts` is a Smolyak grid, then this is the matrix B
        representing  the Smolyak polynomial for the grid. If `pts` is
        not a Smolyak grid, then this is each Smolyak polynomial
        evaluated at each point in pts.

    """
    Ts = cheby2n(pts, m_i(maximum(mu) + 1))::Array{Float64, 3}
    npolys = size(b_inds, 1)::Int64
    npts = size(pts, 1)::Int64
    B = ones(npts, npolys)::Array{Float64, 2}
    for ind = 1:npolys
        for k in 1:d
            b = b_inds[ind,k]::Int64
            for i in 1:npts
                @inbounds B[i,ind] *= Ts[i, k, b]
            end
        end
    end
    return B
end


function dom2cube(pts::Matrix{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
    """
    Takes a point(s) and transforms it(them) into the [-1, 1]^d domain
    """
    dim = size(pts, 2)::Int64

    centers = lb + (ub - lb)./2
    radii = (ub - lb)./2

    cube_points = (pts .- centers')./radii'

    return cube_points
end


function cube2dom(pts::Matrix{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
    """
    Takes a point(s) and transforms it(them) from domain [-1, 1]^d
    back into the desired domain
    """
    dim = size(pts, 2)::Int64

    centers = lb + (ub - lb)./2
    radii = (ub - lb)./2

    dom_points = centers' .+ pts.*radii'

    return dom_points
end


## ----------------- ##
#- Type: SmolyakGrid -#
## ----------------- ##


### type SmolyakGrid
type SmolyakGrid
    d::Int  # number of dimensions
    mu::IntSorV  # density. Int or d element vector of Int
    lb::Array{Float64, 1} # This is the lower bound for the grid (the domain)
    ub::Array{Float64, 1} # This is the upper bound for the grid (the domain)
    cube_grid::Matrix{Float64}  # Smolyak grid on [-1, 1]^d
    grid::Matrix{Float64} # Smolyak grid on original domain
    inds::Array{Any, 1}  # Smolyak indices
    pinds::Matrix{Int64}  # Polynomial indices
    B::Matrix{Float64}  # matrix representing interpoland
    B_fact::LU{Float64}  # LU(copy(B)) -- LU factorization of B

    # Isotropic constructor
    function SmolyakGrid(d::Int, mu::Int, lb::Array{Float64, 1}=ones(d).*-1,
                         ub::Array{Float64, 1}=ones(d))
        if d < 2
            error("You passed d = $d. d must be greater than 1")
        end
        if mu < 1
            error("You passed mu = $mu. mu must be greater than 1")
        end

        inds = smol_inds(d, mu)
        pinds = poly_inds(d, mu, inds)
        cube_grid = build_grid(d, mu, inds)
        grid = cube2dom(cube_grid, lb, ub)
        B = build_B(d, mu, cube_grid, pinds)
        B_fact = lufact(B)

        new(d, mu, lb, ub, cube_grid, grid, inds, pinds, B, B_fact)
    end

    # Anisotropic constructor
    function SmolyakGrid(d::Int, mu::Array{Int, 1},
                         lb::Array{Float64, 1}=ones(d).*-1,
                         ub::Array{Float64, 1}=ones(d))
        if d < 2
            error("You passed d = $d. d must be greater than 1")
        end
        if any(mu .< 1)
            error("You passed mu = $mu. each element must be greater than 1")
        end

        if length(mu) != d
            error("mu must have d elements. It has $(length(mu))")
        end

        inds = smol_inds(d, mu)
        pinds = poly_inds(d, mu, inds)
        cube_grid = build_grid(d, mu, inds)
        grid = cube2dom(cube_grid, lb, ub)
        B = build_B(d, mu, cube_grid, pinds)
        B_fact = lufact(B)

        new(d, mu, lb, ub, cube_grid, grid, inds, pinds, B, B_fact)
    end

    # Default constructor, just in case someone gets brave
    function SmolyakGrid(d, mu, lb, ub, cube_grid, grid, inds, pinds, B, B_fact)
        new(d, mu, lb, ub, cube_grid, grid, inds, pinds, B, B_fact)
    end
end


function SmolyakGrid(d::Int, mu::IntSorV, lb::Real, ub::Real)
    lb = ones(d) * lb
    ub = ones(d) * ub
    return SmolyakGrid(d, mu, lb, ub)
end


function SmolyakGrid(d::Int, mu::IntSorV, lub::Real)
    lub = abs(lub)
    lb = - ones(d) * lub
    ub = ones(d) * lub
    return SmolyakGrid(d, mu, lb, ub)
end


function show(io::IO, sg::SmolyakGrid)
    "show method for SmolyakGrid type"
    npoints = size(sg.cube_grid, 1)
    non_zero_pts = nnz(sg.B)
    pct_non_zero = (non_zero_pts / (npoints ^ 2)) * 100
    if isa(sg.mu, Array)
        mu_str = replace(strip(string(sg.mu)), '\n', " x ")
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
    mu  = sg.mu
    d = sg.d
    mu_str = isa(mu, Array) ? replace(strip(string(mu)),'\n'," \\times ") : mu
    if size(g, 2) == 2
        PyPlot.scatter(g[:, 1], g[:, 2])
        PyPlot.title("Smolyak grid: \$d=$(d), \\; \\mu=$(mu_str)\$")
    elseif size(g, 2) == 3
        PyPlot.scatter3D(g[:, 1], g[:, 2], g[:, 3], "b.")
        PyPlot.title("Smolyak grid: \$d=$(d), \\; \\mu=$(mu_str)\$")
    else
        error("ERROR: can only plot 2d or 3d grids")
    end
end

