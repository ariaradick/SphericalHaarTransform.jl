# SphericalHaarTransform

[![Build Status](https://github.com/ariaradick/SphericalHaarTransform.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ariaradick/SphericalHaarTransform.jl/actions/workflows/CI.yml?query=branch%3Amain)

Transform a function to a spherical Haar wavelet and spherical harmonic basis where each coefficient is given by

$$
f_{n \ell m} = \int d^3 r \: f(\textbf{r}) h_n(r) Y_{\ell m}(\hat{r})
$$

This package uses `FastSphericalHarmonics.jl` to compute the spherical harmonic transform and inverse.

One can transform a function by calling
```jl
Fp = sph_haar_transform(f, nmax, ellmax, umax; method=:twopt)
```
where `f = f([r,θ,φ])` is the function to transform, `nmax = 2^L-1` where `L` is the number of generations in the Haar transform, `ellmax` is the maximum `ell` for the spherical harmonic transform, and `umax` is the maximum value of `r` to consider (since the Haar transform is on a closed interval). To perform the inverse,
```jl
F = sph_haar_inverse(Fp)
```

One can also precompute the function on a grid,
```jl
Rs, θs, φs = sph_haar_points(nmax,ellmax)
F = evaluate_f(f, Rs, θs, φs)
Fp = sph_haar_transform(F)
```