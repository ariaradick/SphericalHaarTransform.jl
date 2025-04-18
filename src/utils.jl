"""Calculates the pairwise sum of a vector `v` and stores the result in `res`"""
function psum!(res,v)
    for i in eachindex(res)
        res[i] = v[2i-1] + v[2i]
    end
end

const GL1 = -1/sqrt(3)
const GL2 = 1/sqrt(3)

function _x_eval_2pt(a,b)
    xm = 0.5*(b-a)
    xp = 0.5*(b+a)

    x1 = xm*GL1 + xp
    x2 = xm*GL2 + xp

    return x1,x2
end

function _x_eval_2pt!(x,u)
    for i in eachindex(u)[1:end-1]
        a = u[i]
        b = u[i+1]

        xm = 0.5*(b-a)
        xp = 0.5*(b+a)

        x[2i-1] = xm*GL1 + xp
        x[2i] = xm*GL2 + xp
    end
end

function _x_eval_2pt(u)
    res = zeros(2*(length(u)-1))
    _x_eval_2pt!(res,u)
    return res
end

"""2-pt Gauss-Legendre quadrature for fast integrals"""
function gauss2pt(f,a,b)
    xm = 0.5*(b-a)
    xp = 0.5*(b+a)

    x1 = xm*GL1 + xp
    x2 = xm*GL2 + xp

    return xm*(f(x2) + f(x1))
end

"""Evaluation point for the single-point rule from arXiv:2502.17547"""
function _xbari3(a,b)
    0.75 * (b^(4) - a^(4)) / (b^(3) - a^(3))
end

"""Single-point integral from arXiv:2502.17547, accurate for very small
integration ranges and faster than 2pt Gauss quadrature for slow functions."""
function singlept(f,a,b)
    (b-a)*f(_xbari3(a,b))
end

function _eval_rint_1!(y, f, x, umax, θ, φ)
    for i in eachindex(x)[1:end-1]
        y[i] = singlept(z -> z^2*f([z*umax, θ, φ]), x[i], x[i+1])
    end
end

function _eval_rint_2!(y, f, x, umax, θ, φ)
    for i in eachindex(x)[1:end-1]
        y[i] = gauss2pt(z -> z^2*f([z*umax, θ, φ]), x[i], x[i+1])
    end
end

function _eval_rint_n!(y, f, x, umax, θ, φ; rtol=1e-6)
    for i in eachindex(x)[1:end-1]
        y[i] = quadgk(z -> z^2*f([z*umax, θ, φ]), x[i], x[i+1]; rtol=rtol)[1]
    end
end

function _eval_rint!(y, f, x, umax, th, ph; method=:twopt, rtol_gquad=1e-6)
    if method == :twopt
        for j in eachindex(ph)
            for i in eachindex(th)
                @views _eval_rint_2!(y[:,i,j], f, x, umax, th[i], ph[j])
            end
        end
    elseif method == :onept
        for j in eachindex(ph)
            for i in eachindex(th)
                @views _eval_rint_1!(y[:,i,j], f, x, umax, th[i], ph[j])
            end
        end
    elseif method == :npt
        for j in eachindex(ph)
            for i in eachindex(th)
                @views _eval_rint_n!(y[:,i,j], f, x, umax, th[i], ph[j]; 
                              rtol=rtol_gquad)
            end
        end
    end
end

function sph_haar_index(n,ell,m)
    return CartesianIndex(n+1, sph_mode(ell,m)...)
end

function sph_haar_points(nmax, ellmax, umax; method=:onept)
    ngen = _n_gens(nmax)
    xi = 0:(2.0^(-ngen)):1
    if method==:twopt
        edges = umax * xi
        us = _x_eval_2pt(edges)
    elseif method==:onept
        us = @. umax*_xbari3(xi[1:end-1], xi[2:end])
    end
    ths, phs = sph_points(ellmax+1)
    return (us, ths, phs)
end

function evaluate_f(f,u,θ,φ)
    res = zeros(length(u), length(θ), length(φ))
    for k in eachindex(φ)
        for j in eachindex(θ)
            for i in eachindex(u)
                res[i,j,k] = f([u[i],θ[j],φ[k]])
            end
        end
    end
    return res
end