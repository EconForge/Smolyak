"""
Created November 2, 2013

Author: Spencer Lyon

This python script implements routines described in 'Numerically Stable and
Accurate Stochastic Simulation Approaches for Solving Dynamic Economic Models'
by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011),
Quantitative Economics 2/2, 173-210 (henceforth, JMM, 2011).'

As such, it is adapted from the file 'Num_Stab_Approx.m'.

"""

function normalize_data(X, Y, intercept=true)
    T, n = size(X)
    X1 = (X[:, 2:n] .- ones(T) * mean(X[:, 2:n], 1)) ./ (ones(T) * std(X[:, 2:n], 1))
    Y1 = (Y .- ones(T) * mean(Y)) ./ (ones(T) * std(Y))
    return X1, Y1
end


function de_normalize(X, Y, beta)
    # De-normalize the coefficients
    T, n = size(X)
    B = Array(Float64, size(beta, 1) + 1, size(beta, 2))
    B[2:end, :] = (1. / std(X[:,2:n], 1)') * std(Y) .* beta
    B[1, :] = mean(Y) - mean(X[:, 2:n], 1) * B[2:end, :]

    return B
end


function OLS(X, Y, normalize=true)
    # Normal OLS.
    # Verified on 11-2-13
    T, n = size(X)
    if normalize
        X1, Y1 = normalize_data(X, Y)
        B1 = inv(X1' * X1) * X1' * Y1
        B = de_normalize(X, Y, B1)
    else
        B = inv(X' * X) * X' * Y
    end

    return B
end


function LS_SVD(X, Y, normalize=true)
    # OLS using singular value decomposition
    # Verified on 11-2-13
    if normalize
        X1, Y1 = normalize_data(X, Y)
        U, S, V = svd(X1, true)
        S_inv = diagm(1. / S)
        B1 = V*S_inv*U'*Y1
        B = de_normalize(X, Y, B1)
    else
        U, S, V = svd(X, true)
        S_inv = diagm(1. / S)
        B = V*S_inv*U'*Y
    end
    return B
end

"""
reload("stable_reg.jl")
using MAT
ml = matread("/Users/sglyon/Desktop/num_stab_approxMATLAB.mat")
X = ml["X"]
Y = ml["Y"]
"""
