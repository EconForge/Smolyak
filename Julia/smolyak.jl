using Iterators
using Cartesian
# import PyPlot
reload("util")
import Base.show

## ------------- ##
#- Generic Tools -#
## ------------- ##

counting_coef(d, mu, i) = (-1) ^ (d + mu - i) * choose(d - 1, d + mu - i)


function num_grid_pts(d::Int, mu::Int)

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
    # Evaluates the first n+1 Chebyshev polynomials of the 'kind' kind at x
    # NOTE: This will only work when kind = 1 or kind = 2
    # NOTE: We evaluate the first n+1, with the first dimension being all 1's

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


function S_n(n::Int)
    # Compute the set S_n. All points are extrema of Chebyshev polynomials

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


function A_n(n::Int)
    if n < 1
        error("DomainError: n must be positive")
    end

    if n == 1
        return [0.]
    elseif n == 2
        return [-1., 1.]
    end

    sn = S_n(n)
    return sn[[2:2:size(sn, 1)]]
end


function A_chain(n::Int)
    # Return Dict of all A_n from n to 1. Keys are n, values are arrays

    sn = S_n(n)

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


function cheb_dict(mu::Int)
    pts = S_n(mu + 1)
    max_poly = length(pts)

    # TODO: Work on figuring out the types here
    # cd = Dict{Int, Dict{Float64, Float64}}()
    cd = Dict()

    for phi_n=1:max_poly
        cd[phi_n] = {pt => t_pt for (pt, t_pt) in
                     zip(pts, cheby_eval(pts, phi_n))}
    end

    return cd
end


function cheby_eval(x::Array{Float64}, n::Int)
    pv = ones(size(x)...)
    cv = x

    if n == 1
        return pv
    elseif n == 2
        return cv
    else
        for i=3:n
            tmp = 2 .* x .* cv - pv
            pv = cv
            cv = tmp
        end
        return cv
    end
end


## --------------- ##
#- Isotropic Grids -#
## --------------- ##


function smol_inds(d::Int, mu::Int)
    p_vals = [1:mu + 1]

    poss_inds = cell(1)  # see below

    for el in Task(()-> comb_with_replacement(p_vals, d))
        if d < sum(el) <= d + mu
            push!(poss_inds, el)
        end
    end

    # cell(1) gives undefined reference in first slot. We remove it here
    poss_inds = poss_inds[2:]

    true_inds = Any[ones(Int64, d)]
    for val in poss_inds
        for el in Task(()->pmute(val))
            push!(true_inds, el)
        end
    end

    return true_inds
end


function base_inds(n::Int)
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


function sparse_grid(d::Int, mu::Int)
    # Use disjoint Smolyak sets to construct Smolyak grid

    true_inds = smol_inds(d, mu)
    An = A_chain(mu + 1)
    return vcat([ cartprod([An[i] for i in el]) for el in true_inds]...)
end


function grid_base(d::Int, mu::Int)
    # Build B matrix for interpolating polynomial

    true_inds = smol_inds(d, mu)
    An = base_inds(mu + 1)
    return vcat([ cartprod([An[i] for i in el]) for el in true_inds]...)
end

### type SGrid

type SGrid
    d::Int
    mu::Int
    grid::Array{Float64, 2}
    inds::Array{Any, 1}

    function SGrid(d::Int, mu::Int, grid::Array{Float64, 2}, inds::Array{Any, 1})
        if d < 2
            error("You passed d = $d. d must be greater than 1")
        end
        if mu < 1
            error("You passed mu = $mu. mu must be greater than 1")
        end
        new(d, mu, grid, inds)
    end
end


# Add convenience constructor to just pass d, mu
SGrid(d::Int, mu::Int) = SGrid(d, mu, sparse_grid(d, mu), smol_inds(d, mu))


function show(io::IO, sg::SGrid)
    npoints = size(sg.grid, 1)
    msg = "Smolyak Grid:\n\td: $(sg.d)\n\tmu: $(sg.mu)\n\tnpoints: $npoints"
    print(io, msg)
end


### Build B-matrix for grid object


function Bmat(sg::SGrid)
    # Compute the matrix B from equation 22 in JMMV 2013
    # Naive translation of dolo.numeric.interpolation.smolyak.SmolyakBasic

    # TODO: check to see if I should swap the k and : in reduce. It might be
    #       more scalable if I do. I would need pass sg.grid' to cheby2n

    Ts = cheby2n(sg.grid, m_i(sg.mu + 1))
    n = size(sg.grid, 1)
    d = sg.d
    b_inds = grid_base(sg.d, sg.mu)
    B = Array(Float64, n, n)
    for ind = 1:n
        B[:, ind] = reduce(.*, {slice(Ts, :, k, b_inds[ind, k]) for k=1:d})
    end
    return B
end


## ----------------- ##
#- Anisotropic Grids -#
## ----------------- ##


function smol_inds(d::Int, mu::Array{Int, 1})
    # Compute indices needed for anisotropic smolyak grid given number of
    # dimensions d and a vector of mu parameters mu

    if length(mu) != d
        error("ValueError: mu must have d elements.")
    end

    mu_max = maximum(mu)
    mup1 = mu + 1

    p_vals = [1:mu_max + 1]

    poss_inds = cell(1)  # see below

    for el in Task(()-> comb_with_replacement(p_vals, d))
        if d < sum(el) <= d + mu_max
            push!(poss_inds, el)
        end
    end

    # cell(1) gives undefined reference in first slot. We remove it here
    poss_inds = poss_inds[2:]

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



# function plot(sg::SGrid)
#     g = sg.grid
#     if size(g, 2) == 2
#         PyPlot.plot(g[:, 1], g[:, 2], "b.")
#         PyPlot.xlim((-1.5, 1.5))
#         PyPlot.ylim((-1.5, 1.5))
#         PyPlot.title("Smolyak grid: \$d=$(sg.d), \\; \\mu=$(sg.mu)\$")
#     elseif size(g, 2) == 3
#         PyPlot.scatter3D(g[:, 1], g[:, 2], g[:, 3], "b.")
#         PyPlot.title("Smolyak grid: \$d=$(sg.d), \\; \\mu=$(sg.mu)\$")
#     else
#         error("ERROR: can only plot 2d or 3d grids")
#     end
# end

#----------
# Test case from dolo
#----------

# fun(x, y) = exp(- x .^ 2 - y .^ 2)
