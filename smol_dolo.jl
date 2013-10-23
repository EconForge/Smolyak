# http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

using Cartesian

function cheb_extrema(n::Int)
    j = [1:n]
    return cos(pi * (j -1) / (n - 1))
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


function cartprod(arrs,
                  out=Array(eltype(arrs[1]),
                            prod([length(a) for a in arrs]),
                            length(arrs)))
    sz = Int[length(a) for a in arrs]
    narrs = length(arrs)
    @forcartesian I sz begin
        k = sub2ind(sz, I)
        for i = 1:narrs
            out[k,i] = arrs[i][I[i]]
        end
    end
    out
end


# function enum(d::Int, l::Int)
#     r = [1:l]
