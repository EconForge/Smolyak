using Iterators
using Cartesian
import PyPlot

choose(n, k) = factorial(n) / (factorial(k) * factorial(n - k))

count_coef(d, mu, i) = (-1) ^ (d + mu - i) * choose(d - 1, d + mu - i)

# takes the ith slice of A in the dth dimension and squeezes along d
slice_sqz(A, d, i) = squeeze(slicedim(A, d, i), d)


function m_i(i::Int)
    if i < 1
        error("DomainError: i must be positive")

    elseif i == 1
        return 1
    else
        return 2^(i - 1) + 1
    end
end


function cheby2n(x, n::Int; kind=1)
    # Evaluates the first n+1 Chebyshev polynomials of the 'kind' kind at x
    # NOTE: This will only work when kind = 1 or kind = 2
    # NOTE: We evaluate the first n+1, with the first dimension being all 1's
    dim = size(x)
    results = zeros(n + 1, dim...)
    results[1, :] = 1.
    results[2, :] = kind * x
    for i=3:n+1
        results[i, :] = (2 * x .* slice_sqz(results, 1, i-1) -
                         slice_sqz(results, 1, i-2))
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


function dense_grid(d::Int, mu::Int)
    # Use nested Smolyak sets to construct Smolyak grid

    p_vals = [1:mu+1]

    # TODO: This can probably be optimized to not be of type Any
    points = Any[]

    # Smolyak.find_points(self)
    for el in product([p_vals for i in [1:d]]...)
        if sum(el) == d + mu
            temp_points = map(S_n, el)
            append!(points, [collect(product(temp_points...))...])
        end
    end

    # Build Smolyak grid
    p_set = Set()

    for elem in points
        push!(p_set, [elem...])
    end

    unique_pts = collect(p_set)
    n_pts = size(unique_pts, 1)
    grid = Array(Float64, n_pts, d)
    for i=1:n_pts
        grid[i, :] = unique_pts[i]
    end

    return grid
end


function num_grid_pts(d::Int, mu::Int)

    if mu == 1
        return int(2 * s - 1)
    elseif mu == 2
        return int(1 + 4*d + 4*d*(d-1)/2)
    else
        return int(1 + 8*d + 12*d*(d-1)/2 + 8*d*(d-1)*(d-2)/6)
    end

end

function sparse_grid(d::Int, mu::Int)
    # Use disjoint Smolyak sets to construct Smolyak grid

    p_vals = ones(Int64, d)*(mu+1)
    An = A_chain(mu + d)

    # TODO: This can probably be optimized to not be of type Any
    points = Any[]

    n_pts = num_grid_pts(d, mu)::Int
    grid = Array(Float64, n_pts, d)
    n = 1

    @forcartesian el p_vals begin
        if d <= sum(el) <= d + mu
            temp_points = [An[i] for i in el]
            for pt in product(temp_points...)
                grid[n, :] = [pt...]
                n += 1
            end
            # append!(points, [collect(product(temp_points...))...])
        end
    end

    # n_pts = size(points, 1)::Int
    # grid = Array(Float64, n_pts, d)
    # for i=1:n_pts
    #     grid[i, :] = [points[i]...]
    # end

    return grid
end


function plot_grid(g)
    if size(g, 2) == 2
        PyPlot.plot(g[:, 1], g[:, 2], "b.")
        PyPlot.xlim((-1.5, 1.5))
        PyPlot.ylim((-1.5, 1.5))
    elseif size(g, 2) == 3
        PyPlot.scatter3D(g[:, 1], g[:, 2], g[:, 3], "b.")
    else
        throw("ERROR: can only plot 2d or 3d grids")
    end
end
