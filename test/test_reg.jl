function all_close(x::Array, y::Array, rtol::Float64=1.e-5, atol::Float64=1.e-8)
    """
    Returns True if two arrays are element-wise equal within a tolerance.

    References
    ==========
    This function, and docstring, are taken directly from numpy.allclose

    See that documentation for more info

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

reload("../stable_reg.jl")
using MAT

# load in matlab results
ml = matread("num_stab_approxMATLAB.mat")

# Grab data
X = ml["X"]
Y = ml["Y"]

# compare to matlab answers

# OLS
@assert all_close(OLS(X, Y), ml["b_ols_yes"])
@assert all_close(OLS(X, Y, false), ml["b_ols_no"])

# SVD
@assert all_close(LS_SVD(X, Y), ml["b_ols_yes"])
@assert all_close(LS_SVD(X, Y, false), ml["b_ols_no"])

# LAD-PP
@assert all_close(LAD_PP(X, Y), ml["b_ladpp_yes"])
@assert all_close(LAD_PP(X, Y, false), ml["b_ladpp_no"])

# LAD-DP
@assert all_close(LAD_DP(X, Y), ml["b_laddp_yes"])
@assert all_close(LAD_DP(X, Y, false), ml["b_laddp_no"])

# Regularized least square Tikhonov
@assert all_close(RLS_T(X, Y), ml["b_ols_yes"])

# Regularized SVD Tikhonov
@assert all_close(RLS_TSVD(X, Y), ml["b_ols_yes"])
