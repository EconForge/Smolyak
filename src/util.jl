function pmute(a::Union(Array{Float64, 1}, Array{Int64, 1}))
    # Return all unique permutations of the vector a, which must be a 1d
    # numerical array.

    sort!(a)  # Sort so we can deal with repeated elements

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
                    j -= 1  # j--
                end

                a[i], a[j] = a[j], a[i]  # swap(a[j], a[i])
                t = a[(i + 1):end]
                reverse!(t)
                a[(i + 1):end] = t

                # Output current.
                produce(copy(a))

                break  # next.
            end

            if i == first
                reverse!(a)
                produce(copy(a))
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


function cartprod(arrs, out=Array(eltype(arrs[1]),
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


function tensordot{T, S, N}(a::Array{T, N}, b::Array{S, N}, axes::Array{Int, 1})
    if length(axes) != 2
        error("Haven't written this one yet")
    end
    a, b = axes
end

## ------------- ##
#- Testing Tools -#
## ------------- ##

function all_close(x::Array, y::Array, rtol::Float64=1.e-5, atol::Float64=1.e-8)
    """
    Returns True if two arrays are element-wise equal within a tolerance.

    The tolerance values are positive, typically very small numbers.  The
    relative difference (`rtol` * abs(`b`)) and the absolute difference
    `atol` are added together to compare against the absolute difference
    between `a` and `b`.

    If either array contains one or more NaNs, False is returned.
    Infs are treated as equal if they are in the same place and of the same
    sign in both arrays.

    Parameters
    ----------
    x, y : array_like
        Input arrays to compare.
    rtol : float
        The relative tolerance parameter (see Notes).
    atol : float
        The absolute tolerance parameter (see Notes).

    Returns
    -------
    allclose : bool
        Returns True if the two arrays are equal within the given
        tolerance; False otherwise.

    Notes
    -----
    If the following equation is element-wise True, then allclose returns
    True.

     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    The above equation is not symmetric in `a` and `b`, so that
    `allclose(a, b)` might be different from `allclose(b, a)` in
    some rare cases.

    Examples
    --------
    julia> allclose([1e10,1e-7], [1.00001e10,1e-8])
    false
    julia> allclose([1e10,1e-8], [1.00001e10,1e-9])
    true
    julia> allclose([1e10,1e-8], [1.0001e10,1e-9])
    false

    References
    ==========
    This function, and docstring, are taken directly from numpy.allclose

    """
    xinf = isinf(x)
    yinf = isinf(y)
    if any(xinf) || any(yinf)
        if !all(xinf .== yinf)
            return false
        end

        if !all(x[xinf] .== y[yinf])
            return false
        end

        x = x[~xinf]
        y = y[~xinf]
    end

    return all(.<=(abs(x - y), atol + rtol * abs(y)))
end

