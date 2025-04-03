# SphericalHaarTransform

[![Build Status](https://github.com/ariaradick/SphericalHaarTransform.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ariaradick/SphericalHaarTransform.jl/actions/workflows/CI.yml?query=branch%3Amain)

Transform a function to a spherical Haar wavelet and spherical harmonic basis where each coefficient is given by

$$
f_{n \ell m} = \int d^3 u f(\vec{u}) h_n(u/u_{\textrm{max}}) Y_{\ell m}(\hat{u})
$$

This package uses [`FastSphericalHarmonics.jl`](https://github.com/eschnett/FastSphericalHarmonics.jl) to compute the spherical harmonic transform and inverse.

One can transform a function by calling
```jl
Fp = sph_haar_transform(f, nmax, ellmax, umax; method=:twopt)
```
where `f = f([r,θ,φ])` is the function to transform, `nmax = 2^L-1` where `L` is the number of generations in the Haar transform, `ellmax` is the maximum `ell` for the spherical harmonic transform, and `umax` is the maximum value of `u` to consider (since the Haar transform is on a closed interval). To perform the inverse,
```jl
F = sph_haar_inverse(Fp)
```

One can also precompute the function on a grid,
```jl
us, θs, φs = sph_haar_points(nmax,ellmax,umax)
F = evaluate_f(f, us, θs, φs)
Fp = sph_haar_transform(F)
```

To access a specific (n,ell,m) coefficient use
```jl
Fp[sph_haar_index(n,ell,m)]
```