using SphericalHaarTransform
using Test

function sph_to_cart(x_sph)
    x, θ, φ = x_sph
    rx = x*sin(θ)*cos(φ)
    ry = x*sin(θ)*sin(φ)
    rz = x*cos(θ)
    return [rx, ry, rz]
end

function fj2(nj, qLj)
    if qLj == 0.0
        if nj == 1
            return 1.0
        else
            return 0.0
        end
    end
    qlp = abs(qLj)/π
    t1 = sinc(0.5*(qlp - nj + 1.0)) / (1 + (nj-1)/qlp)
    t2 = sinc(0.5*(qlp - nj - 1.0)) / (1 + (nj+1)/qlp)
    return (t1+t2)^2
end

@testset "SphericalHaarTransform.jl" begin
    qBohr = 511e3 / 137

    nmax = 2^10-1
    ellmax = 36
    umax = 10*qBohr

    a0 = 1 / (qBohr)
    Length = [4.0, 7.0, 10.0] .* a0
    fs2(Lvec, nz, q_xyz) = prod(fj2.([1, 1, nz], Lvec.*q_xyz))
    fs2_model(qSph) = fs2(Length, 2, sph_to_cart(qSph))

    ui, θj, φk = sph_haar_points(nmax, ellmax, umax)
    F = evaluate_f(fs2_model, ui, θj, φk)
    Fp = sph_haar_transform(F)
    Fpp = sph_haar_inverse(Fp)

    C1 = sph_haar_transform(fs2_model, nmax, ellmax, umax; method=:onept)
    C2 = sph_haar_transform(fs2_model, nmax, ellmax, umax; method=:twopt)
    Cn = sph_haar_transform(fs2_model, nmax, ellmax, umax; method=:npt)

    @test isapprox(Fp, C1; rtol=1e-3)
    @test isapprox(C2, C1; rtol=1e-3)
    @test isapprox(C2, Cn; rtol=1e-3)
    @test isapprox(F, Fpp; rtol=1e-3)
end
