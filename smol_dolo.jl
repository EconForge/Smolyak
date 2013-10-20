# http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

function cheb_extrema(n::Int)
    j = [1:n]
    return cos(pi * (j -1) / (n - 1))
end

function my_repeat(a, n)
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
        println("HERE!")
        # out = Array(Float64, n, length(arrs))
        out = zeros(n, length(arrs))
    end

    m = int(n / size(arrs[1], 1))
    out[:, 1] = my_repeat(arrs[1], m)
    println("OUTSIDE LOOP: out =\n$out\n")

    if length(arrs[2:]) > 0
        cartesian(arrs[2:], out=out[1:m, 2:])
        called = called + 1
        println("called = $called")
        println("Above LOOP: out =\n$out\n")
        for j = 1:size(arrs[1], 1)-1
            println("j = $j, m = $m")
            println("(j*m + 1):(j+1)*m = \n$([(j*m + 1):(j+1)*m])")
            println("out = \n$out")
            out[(j*m + 1):(j+1)*m, 2:] = out[1:m, 2:]
        end
    end

    return out
end


# function enum(d::Int, l::Int)
#     r = [1:l]
