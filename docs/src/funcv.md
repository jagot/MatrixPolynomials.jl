# Functions of matrices

```@docs
MatrixPolynomials.FuncV
MatrixPolynomials.FuncV(f::Function, A, m, t=one(eltype(A)); distribution=:leja, leja_multiplier=100, tol=1e-15, kwargs...)
```

## Spectral ranges and shapes

```@docs
MatrixPolynomials.spectral_range
MatrixPolynomials.Shape
MatrixPolynomials.Line
MatrixPolynomials.Rectangle
```
