using Iterators
import PyPlot

function S_n(n::Int)
    # Compute the set S_n. All points are extrema of Chebyshev polynomials
    if n < 1
        error("DomainError: n must be positive")
    end

    if n == 1
        return [0]
    end

    m_i = 2^(n - 1) + 1
    j = [1:m_i]
    pts = cos(pi * (j - 1.) / (m_i - 1))::Array{Float64, 1}
    for i = 1:size(pts, 1)
        pts[i] = abs(pts[i]) < 1e-12 ? 0.0 : pts[i]
    end
    return pts

end

function A_n(n::Int)
    # make a_n
    if n < 1
        error("DomainError: n must be positive")
    elseif n == 1
        return [0]
    end

    return setdiff(S_n(n), S_n(n-1))
end




function dense_grid(d::Int, mu::Int)

    p_vals = [1:mu+1]

    # TODO: This can probably be optimized
    points = Any[]

    # Smolyak.find_points(self)
    for el in product([p_vals for i in [1:d]]...)
        if sum(el) == d + mu
            temp_points = map(S_n, el)
            append!(points, [collect(product(temp_points...))...])
        end
    end

    p_set = Set()

    # Build Smolyak grid

    for elem in points
        push!(p_set, [elem...])
    end

    p_set = unique(p_set)
    n_pts = size(p_set, 1)
    grid = Array(Float64, n_pts, d)
    for i=1:n_pts
        grid[i, :] = p_set[i]
    end
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

