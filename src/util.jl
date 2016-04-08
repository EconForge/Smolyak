#=

Utility functions used throughout the library.

=#

function pmute(a::AbstractVector)
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
    k = 1
    for v in product(arrs...)
        i = 1
        for el in v
            out[k, i] = el
            i += 1
        end
        k += 1
    end

    return out
end

immutable WithReplacementCombinations{T}
    a::T
    t::Int
end

Base.eltype{T}(::Type{WithReplacementCombinations{T}}) = Vector{eltype(T)}

Base.length(c::WithReplacementCombinations) = binomial(length(c.a)+c.t-1, c.t)

"generate all combinations with replacement of size t from an array a."
with_replacement_combinations(a, t::Integer) = WithReplacementCombinations(a, t)

Base.start(c::WithReplacementCombinations) = [1 for i in 1:c.t]
function Base.next(c::WithReplacementCombinations, s)
    n = length(c.a)
    t = c.t
    comb = [c.a[si] for si in s]
    if t > 0
        s = copy(s)
        changed = false
        for i in t:-1:1
            if s[i] < n
                s[i] += 1
                for j in (i+1):t
                    s[j] = s[i]
                end
                changed = true
                break
            end
        end
        !changed && (s[1] = n+1)
    else
        s = [n+1]
    end
    (comb, s)
end
Base.done(c::WithReplacementCombinations, s) = !isempty(s) && s[1] > length(c.a) || c.t < 0

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
