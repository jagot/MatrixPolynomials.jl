# MatrixPolynomials.jl

The main purpose of this package is to provide polynomial
approximations to ``f(\mat{A})``, i.e. the function of a matrix
``\mat{A}`` for which the field-of-values ``W(\mat{A}) \subset
\Complex`` (or equivalently the distribution of eigenvalues) is known
_a priori_. If this is the case, a polynomial approximation ``p(z)
\approx f(z)`` for ``z \in W(\mat{A})`` can be constructed, and this
can subsequently be used, substituting ``\mat{A}`` for ``z``. This is
in contrast to Krylov-based methods, where the matrix polynomials are
generated on-the-fly, without any prior knowledge of ``W(\mat{A})``
(even though knowledge _can_ be used to speed up the convergence of
the Krylov iterations).

## Index

```@index
```
