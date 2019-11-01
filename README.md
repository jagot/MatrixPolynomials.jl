# MatrixPolynomials.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jagot.github.io/MatrixPolynomials.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jagot.github.io/MatrixPolynomials.jl/dev)
[![Build Status](https://travis-ci.com/jagot/MatrixPolynomials.jl.svg?branch=master)](https://travis-ci.com/jagot/MatrixPolynomials.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jagot/MatrixPolynomials.jl?svg=true)](https://ci.appveyor.com/project/jagot/MatrixPolynomials-jl)
[![Codecov](https://codecov.io/gh/jagot/MatrixPolynomials.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jagot/MatrixPolynomials.jl)

This package aids in the computation of the action of a matrix
polynomial on a vector, i.e. `p(A)v`, where `A` is a (square) matrix
(or a linear operator) that is supplied to the polynomial `p`. The
matrix polynomial `p(A)` is never formed explicitly, instead only its
action on `v` is evaluated. This is commonly used in time-stepping
algorithms for ordinary differential equations (ODEs) and discretized
partial differential equations (PDEs) where `p` is an approximation of
the exponential function (or the related `φ` functions:
`φ₀(z) = exp(z)`, `φₖ₊₁ = [φₖ(z)-φₖ(0)]/z`, `φₖ(0)=1/k!`) on the
field-of-values of the matrix `A`, which for the methods in this
package needs to be known before-hand.

## Alternatives

Other packages with similar goals, but instead based on matrix
polynomials found via Krylov iterations are

- https://github.com/JuliaDiffEq/ExponentialUtilities.jl
- https://github.com/Jutho/KrylovKit.jl

Krylov iterations do not need to know the field-of-values of the
matrix `A` before-hand, instead, an orthogonal basis is built-up
on-the-fly, by repeated action of `A` on test vectors: `Aⁿ*v`. This
process is however very sensitive to the condition number of `A`,
something that can be alleviated by iterating a shifted and inverted
matrix instead: `(A-σI)⁻¹` (rational Krylov). Not all matrices/linear
operators are easily inverted/factorized, however.

Moreover, the Krylov iterations for general matrices (then called
Arnoldi iterations) require long-term recurrences with mutual
orthogonalization along with inner products, all of which can be
costly to compute. Finally, a subspace approximation of the polynomial
`p` of a upper Hessenberg matrix needs to computed. The
real-symmetric/complex-Hermitian case (Lanczos iterations) reduces to
three-term recurrences and a tridiagonal subspace matrix. In contrast,
the polynomial methods of this packages two-term recurrences only, no
orthogonalization (and hence no inner products), and finally no
evaluation of the polynomial on a subspace matrix. This could
potentially mean that the methods are easier to implement on a GPU.
