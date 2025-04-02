function sph_haar_transform(F)
    ngen = _n_gens(size(F)[1]-1)
    dx = 2.0^(-ngen)
    xs = 0:dx:1
    xi = _xbari3.(xs[1:end-1], xs[2:end])

    Fp = dx .* F
    for k in 1:size(F)[3]
        for j in 1:size(F)[2]
            @. Fp[:,j,k] *= xi^2
            @views haar3D_transform!(Fp[:,j,k])
        end
    end

    Fpp = zeros(size(Fp))

    for i in 1:size(F)[1]
        Fpp[i,:,:] = sph_transform(Fp[i,:,:])
    end
    return Fpp
end

function sph_haar_transform(f, nmax, ellmax, umax; method=:twopt, 
                            rtol_gquad=1e-6)
    ngen = _n_gens(nmax)
    xs = 0:(2.0^(-ngen)):1
    ths,phs = sph_points(ellmax+1)

    F = zeros(length(xs)-1, length(ths), length(phs))
    _eval_rint!(F, f, xs, umax, ths, phs; method=method, rtol_gquad=rtol_gquad)

    Fp = zeros(size(F))

    for i in 1:size(F)[1]
        Fp[i,:,:] = sph_transform(F[i,:,:])
    end

    for k in 1:size(F)[3]
        for j in 1:size(F)[2]
            @views haar3D_transform!(Fp[:,j,k])
        end
    end

    return Fp
end

function sph_haar_inverse(Fp)
    Fpp = zeros(size(Fp))
    for i in 1:size(Fp)[1]
        Fpp[i,:,:] = sph_evaluate(Fp[i,:,:])
    end
    res = zeros(size(Fp))
    for k in 1:size(Fp)[3]
        for j in 1:size(Fp)[2]
            @views haar3D_inverse!(res[:,j,k], Fpp[:,j,k])
        end
    end
    return res
end
