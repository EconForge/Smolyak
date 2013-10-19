using Iterators
import PyPlot


function S_n(n::Int)
    # Compute the set S_n. All points are extrema of Chebyshev polynomials
    if n < 1
        error("DomainError: n must be positive")
    end

    if n == 1
        return [0.]
    end

    m_i = 2^(n - 1) + 1
    j = [1:m_i]
    pts = cos(pi * (j - 1.) / (m_i - 1))::Array{Float64, 1}
    for i = 1:size(pts, 1)
        pts[i] = abs(pts[i]) < 1e-12 ? 0.0 : pts[i]
    end
    return pts
end


# function A_n(n::Int)
#     # make a_n
#     if n < 1
#         error("DomainError: n must be positive")
#     elseif n == 1
#         return [0]
#     end

#     return setdiff(S_n(n), S_n(n-1))
# end


function A_n_prime(n::Int)
    sn = S_n(n)
    return sn[[2:2:size(sn, 1)]]
end


function A_chain(n::Int)
    # Return an Any array of all A_n from n to 1
    sn = S_n(n)
    a = Any[]
    for i=n-1:-1:2
        #pass
        push!(a, sn[[2:2:size(sn, 1)]])
        sn = sn[[1:2:size(sn, 1)]]
    end
    # add A_2 and A_1
    push!(a, [-1, 1])
    push!(a, [0])
    return a
end


function dense_grid(d::Int, mu::Int)

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


function sparse_grid(d::Int, mu::Int)

    # We only need to check values up to mu + 1 because each of the other
    # (n-1) dimensions is at least 1 and they need to be bounded above by
    # d + mu
    p_vals = [1:mu + 1]
    An = A_chain(mu + 1)

    # TODO: This can probably be optimized to not be of type Any
    points = Any[]

    # TODO: Re-write this function. The key innovation behind the disjoint
    #       method is that we will no longer obtain duplicate points when
    #       evaluating the tensor product. That means we just need to stack
    #       all the A_i together.

    # NOTE: There is definitely a way to use A_chain to compute all the A_n.
    #       Then we can "stack" the points in the A_n that satisfy
    #       d <= |i| <= d + u. We will not have to do any set operations!!
    #       (either in constructing A_n *or* in finding unique points)

    # Check all the combinations of values that come from possible values
    for el in product([p_vals for i in [1:d]]...)
        if d <= sum(el) <= d + mu
            # We want to grab all the corresponding sets in An
            # We pull from the end because the sets are made in reverse
            # order.
            temp_points = [An[end-el[i]] for i in el]
            append!(points, [collect(product(temp_points...))...])
        end
    end

    return points

    # Build Smolyak grid
    # p_set = Set()

    # for elem in points
    #     push!(p_set, [elem...])
    # end

    # unique_pts = collect(p_set)
    # n_pts = size(unique_pts, 1)
    # grid = Array(Float64, n_pts, d)
    # for i=1:n_pts
    #     grid[i, :] = unique_pts[i]
    # end

    # return grid
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
