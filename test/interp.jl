using Smolyak

@testset "interpolation with derivatives" begin
    d = 2
    N = 2
    mu = 2
    f(x) = vec(sum(x, 2))

    # ub = rand(d) + 6
    # lb = rand(d) - 5
    sg = SmolyakGrid(d, mu, 0, 1)

    values = f(sg.grid)
    si = SmolyakInterp(sg, values)
    pts = rand(N, d)

    @test maxabs(f(pts) - @inferred(si(pts))) < 1e-12

    # now get derivatives
    res, res_s = interp_withderiv(si, pts)
    res_s

    @test size(res) == (N,)
    @test size(res_s) == (N,d)
end
