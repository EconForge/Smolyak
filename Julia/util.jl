# module util

import Cartesian

# export choose, slice_sqz, pmute, cartprod, cartesian, comb_with_replacement

choose(n, k) = factorial(n) / (factorial(k) * factorial(n - k))

# takes the ith slice of A in the dth dimension and squeezes along d
slice_sqz(A, d, i) = squeeze(slicedim(A, d, i), d)


function pmute(a::Union(Array{Float64, 1}, Array{Int64, 1}))
    # Return all unique permutations of the vector a, which must be a 1d
    # numerical array.

    sort!(a)

    # produce(a)  # we do this below, just before returning

    # just initializing here so these are available in all `while` scopes
    i = 0
    first = 1
    alen = length(a)

    while true
        i = alen
        while true
            i -= 1

            if a[i] < a[i + 1]
                j = alen
                while a[i] >= a[j]
                    j -= 1 # j--
                end
                a[i], a[j] = a[j], a[i] # swap(a[j], a[i])
                t = a[(i + 1):end]
                reverse!(t)
                a[(i + 1):end] = t

                # Output current.
                produce(collect(tuple(a...)))

                break # next.
            end

            if i == first
                reverse!(a)
                produce(collect(tuple(a...)))
                return
            end
        end
    end
end


function my_repeat(a, n)
    # mimics numpy.repeat
    m = size(a, 1)
    out = Array(eltype(a), n * m)
    out[1:n] = a[1]
    for i=2:m
        out[(i-1)*n+1:i*n] = a[i]
    end
    return out
end


function cartprod(arrs,
                  out=Array(eltype(arrs[1]),
                            prod([length(a) for a in arrs]),
                            length(arrs)))
    sz = Int[length(a) for a in arrs]
    narrs = length(arrs)
    Cartesian.@forcartesian I sz begin
        k = sub2ind(sz, I)
        for i = 1:narrs
            out[k,i] = arrs[i][I[i]]
        end
    end
    out
end


function cartesian(arrs; out=None)
    called=1
    dtype = eltype(arrs[1])

    n = prod([size(i, 1) for i in arrs])::Int

    if is(out, None)
        out = Array(dtype, n, length(arrs))
    end

    m = int(n / size(arrs[1], 1))
    out[:, 1] = my_repeat(arrs[1], m)

    if length(arrs[2:]) > 0
        out_end = size(out, 2)
        cartesian(arrs[2:], out=sub(out, 1:m, 2:out_end))
        for j = 1:size(arrs[1], 1)-1
            out[(j*m + 1):(j+1)*m, 2:] = out[1:m, 2:]
        end
    end

    return out
end

function comb_with_replacement(itr, r::Int)
    # Combinations with replacement. From algorithm in python docs.

    pool = tuple(itr...)
    n = length(pool)::Int
    indicies = ones(Int64, r)
    produce([pool[i] for i in indicies])

    # NOTE: need to define done and i here so they are available in all parts of
    #       the while loop
    done = false
    i = 0

    while true
        for i=r:-1:1
            if indicies[i] != n
                done = false
                break
            else
                done = true
            end
        end

        if is(done, true)
            return
        end
        indicies[i:] = indicies[i] + 1
        produce([pool[i] for i in indicies])
    end
end


# end # module
