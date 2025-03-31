module SphericalHaarTransform

export sph_haar_transform

using QuadGK
using FastSphericalHarmonics
# using Healpix

include("utils.jl")
include("haar.jl")

function sph_haar_transform(f, nmax, ellmax, umax; method=:twopt, 
                            rtol_gquad=1e-6)
    ngen = _n_gens(nmax)
    xs = 0:(2.0^(-ngen)):1
    ths,phs = sph_points(ellmax+1)

    F = zeros(length(xs)-1, length(ths), length(phs))
    _eval_rint!(F, f, xs, umax, ths, phs; method=method, rtol_gquad=rtol_gquad)

    Fp = zeros(size(F))
    for i in 1:size(F)[1]
        @views Fp[i,:,:] = sph_transform(F[i,:,:])
    end

    for k in 1:size(F)[3]
        for j in 1:size(F)[2]
            @views discrete_haar_3D!(Fp[:,j,k])
        end
    end
    return Fp
end

end
