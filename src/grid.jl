#=

Chase Coleman (ccoleman@stern.nyu.edu)
Spencer Lyon (spencer.lyon@stern.nyu.edu)

References
----------

Judd, Kenneth L, Lilia Maliar, Serguei Maliar, and Rafael Valero. 2013.
    "Smolyak Method for Solving Dynamic Economic Models: Lagrange
    Interpolation, Anisotropic Grid and Adaptive Domain".

Krueger, Dirk, and Felix Kubler. 2004. "Computing Equilibrium in OLG
    Models with Stochastic Production." Journal of Economic Dynamics and
    Control 28 (7) (April): 1411-1436.

=#

# union type to handle isotropic/anisotropic mu parameters
typealias IntSorV Union{Int, Vector{Int}}

num_grid_points(d::Int, mu::Int) =
    mu == 1 ? 2d - 1 :
    mu == 2 ? Int(1 + 4d + 4d*(d-1)/2 ):
    mu == 3 ? Int(1 + 8d + 12d*(d-1)/2 + 8d*(d-1)*(d-2)/6) :
    error("We only know the number of grid points for mu ∈ [1, 2, 3]")

m_i(i::Int) = i < 0 ?  error("DomainError: i must be positive") :
              i == 0 ? 0 :
              i == 1 ? 1 : 2^(i - 1) + 1

function cheby2n{T<:Number}(x::Array{T}, n::Int, kind::Int=1)
    out = Array(T, size(x)..., n+1)
    cheby2n!(out, x, n, kind)
end

@noinline function cheby2n!{T<:Number,N}(out::Array{T}, x::Array{T,N}, n::Int,
                                         kind::Int)
    if size(out) != tuple(size(x)..., n+1)
        error("out must have dimensions $(tuple(size(x)..., n+1))")
    end

    R = CartesianRange(size(x))
    # fill first element with ones
    @inbounds @simd for I in R
        out[I, 1] = one(T)
        out[I, 2] = kind*x[I]
    end

    @inbounds for i in 3:n+1
        @simd for I in R
            out[I, i] = 2x[I] * out[I, i - 1] - out[I, i - 2]
        end
    end
    out
end

function s_n(n::Int)
    n < 1 && error("DomainError: n must be positive")

    if n == 1
        return [0.0]
    end

    m = m_i(n)
    j = 1:m
    pts = cos(pi * (j - 1.0) / (m - 1.0))
    @inbounds @simd for i in eachindex(pts)
        pts[i] = abs(pts[i]) < 1e-12 ? 0.0 : pts[i]
    end

    pts
end


function a_n(n::Int)
    n < 1 && error("DomainError: n must be positive")

    if n == 1
        return [0.0]
    elseif n == 2
        return [-1.0, 1.0]
    end

    sn = s_n(n)

    sn[2:2:end]
end

function a_chain(n::Int)
    sn = s_n(n)
    a = Dict{Int,Vector{Float64}}()
    sizehint!(a, n)

    # These are constant and don't follow the pattern.
    a[1] = [0.0]
    a[2] = [-1.0, 1.0]

    for i=n:-1:3
        # push!(a, sn[[2:2:size(sn, 1)]])
        a[i] = sn[2:2:end]
        sn = sn[1:2:end]
    end

    a
end

function phi_chain(n::Int)
    max_ind = m_i(n)
    phi = Dict{Int, UnitRange{Int64}}()
    phi[1] = 1:1
    phi[2] = 2:3
    low_ind = 4  # current lower index

    for i = 3:n
        high_ind = m_i(i)
        phi[i] = low_ind:high_ind
        low_ind = high_ind + 1
    end

    phi
end

## ---------------------- ##
#- Construction Utilities -#
## ---------------------- ##

# Isotropic inds
function smol_inds(d::Int, mu::Int)

    p_vals = 1:(mu+1)

    # PERF: size_hint here if it is slow
    poss_inds = Vector{Int}[]

    for el in with_replacement_combinations(p_vals, d)
        if d < sum(el) <= d + mu
            push!(poss_inds, el)
        end
    end

    # PERF: size_hint here if it is slow
    true_inds = Vector{Int}[ones(Int64, d)]  # we will always have (1, 1, ...,  1)
    for val in poss_inds
        for el in @task pmute(val)
            push!(true_inds, el)
        end
    end

    return true_inds
end

# Ansotropic inds
function smol_inds(d::Int, mu::Vector{Int})
    # Compute indices needed for anisotropic smolyak grid given number of
    # dimensions d and a vector of mu parameters mu

    length(mu) != d &&  error("ValueError: mu must have d elements.")

    mu_max = maximum(mu)
    mup1 = mu + 1

    p_vals = 1:(mu_max+1)

    poss_inds = Vector{Int}[]

    for el in with_replacement_combinations(p_vals, d)
        if d < sum(el) <= d + mu_max
            push!(poss_inds, el)
        end
    end

    true_inds = Vector{Int}[ones(Int64, d)]
    for val in poss_inds
        for el in @task pmute(val)
            if all(el .<= mup1)
                push!(true_inds, el)
            end
        end
    end

    return true_inds
end

function poly_inds(d::Int, mu::IntSorV, inds::Vector{Vector{Int}}=smol_inds(d, mu))
    phi_n = phi_chain(maximum(mu) + 1)
    vcat([cartprod([phi_n[i] for i in el]) for el in inds]...)::Matrix{Int}
end

# Build grid
function build_grid(d::Int, mu::IntSorV, inds::Vector{Vector{Int}}=smol_inds(d, mu))
    An = a_chain(maximum(mu) + 1)  # use maximum in case mu is Vector
    vcat([cartprod([An[i] for i in el]) for el in inds]...)::Matrix{Float64}
end

# Build B-matrix
function build_B!(out::AbstractMatrix, d::Int, mu::IntSorV,
                  pts::Matrix{Float64}, b_inds::Matrix{Int64})

    npolys = size(b_inds, 1)
    npts = size(pts, 1)
    size(out) == (npts, npolys) || error("Out should be size $((npts, npolys))")
    Ts = cheby2n(pts, m_i(maximum(mu) + 1))

    @inbounds for ind in 1:npolys, k in 1:d
        b = b_inds[ind, k]
        for i in 1:npts
            out[i, ind] *= Ts[i, k, b]
        end
    end

    return out
end

build_B(d::Int, mu::IntSorV, pts::Matrix{Float64}, b_inds::Matrix{Int64}) =
    build_B!(ones(size(pts, 1), size(b_inds, 1)), d, mu, pts, b_inds)

function dom2cube(pts::AbstractMatrix, lb::AbstractVector, ub::AbstractVector)
    centers = lb + (ub - lb)./2
    radii = (ub - lb)./2
    cube_points = (pts .- centers')./radii'

    cube_points
end

function cube2dom(pts::AbstractMatrix, lb::AbstractVector, ub::AbstractVector)
    centers = lb + (ub - lb)./2
    radii = (ub - lb)./2
    dom_points = centers' .+ pts.*radii'

    dom_points
end

## ----------------- ##
#- Type: SmolyakGrid -#
## ----------------- ##

### type SmolyakGrid
# TODO: decide if we should keep B or sparse(B)
type SmolyakGrid{Tmu<:IntSorV,Tlb<:Number,Tub<:Number}
    d::Int  # number of dimensions
    mu::Tmu  # density. Int or d element vector of Int
    lb::Vector{Tlb} # This is the lower bound for the grid (the domain)
    ub::Vector{Tub} # This is the upper bound for the grid (the domain)
    cube_grid::Matrix{Float64}  # Smolyak grid on [-1, 1]^d
    grid::Matrix{Float64} # Smolyak grid on original domain
    inds::Vector{Vector{Int}}  # Smolyak indices
    pinds::Matrix{Int64}  # Polynomial indices
    B::Matrix{Float64}  # matrix representing interpoland
    B_fact::Base.LinAlg.LU{Float64,Matrix{Float64}}  # LU factorization of B

    function SmolyakGrid(d::Int, mu::IntSorV, lb::Vector=-ones(d),
                         ub::Vector=ones(d))
        d < 2 && error("You passed d = $d. d must be greater than 1")
        mu < 1 && error("You passed mu = $mu. mu must be greater than 1")
        if length(mu) > 1
            # working on building an anisotropic grid
            if length(mu) != d
                error("mu must have d elements. It has $(length(mu))")
            end
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

# outer constructor to infer types
function SmolyakGrid{Tmu<:IntSorV,Tlb<:Number,Tub<:Number}(d::Int,
                                                           mu::Tmu,
                                                           lb::Vector{Tlb}=-ones(d),
                                                           ub::Vector{Tub}=ones(d))
    SmolyakGrid{Tmu,Tlb,Tub}(d, mu, lb, ub)
end

function SmolyakGrid(d::Int, mu::IntSorV, lb::Number, ub::Number)
    lb = fill(lb, d)
    ub = fill(ub, d)
    SmolyakGrid(d, mu, lb, ub)
end

function SmolyakGrid(d::Int, mu::IntSorV, lub::Number)
    lub = abs(lub)
    lb = fill(-lub, d)
    ub = fill(lub, d)
    SmolyakGrid(d, mu, lb, ub)
end

function Base.show(io::IO, sg::SmolyakGrid)
    npoints = size(sg.cube_grid, 1)
    non_zero_pts = countnz(sg.B)
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
