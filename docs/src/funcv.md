# Functions of matrices

```@docs
MatrixPolynomials.FuncV
MatrixPolynomials.FuncV(f::Function, A, m::Integer, t=one(eltype(A)); distribution=:leja, leja_multiplier=100, tol=1e-15, kwargs...)
LinearAlgebra.mul!(w, funcv::MatrixPolynomials.FuncV, v)
```

## Spectral ranges and shapes

```@docs
MatrixPolynomials.spectral_range
MatrixPolynomials.hermitian_spectral_range
```

### Shapes

The shapes in the complex plane are mainly used to generate suitable
distributions of [Leja points](@ref), which are in turn used to
generate [Newton polynomials](@ref) that efficiently approximate
various functions on the field-of-values of a matrix ``\mat{A}`` which
is contained within the spectral shape.

```@docs
MatrixPolynomials.Shape
```

#### Lines

```@docs
MatrixPolynomials.Line
Base.:(*)(n::Number, l::MatrixPolynomials.Line)
Base.range(l::MatrixPolynomials.Line, n)
Statistics.mean(l::MatrixPolynomials.Line)
Base.union(a::MatrixPolynomials.Line, b::MatrixPolynomials.Line)
```

#### Rectangles

```@docs
MatrixPolynomials.Rectangle
Base.:(*)(n::Number, l::MatrixPolynomials.Rectangle)
Base.range(l::MatrixPolynomials.Rectangle, n)
Statistics.mean(l::MatrixPolynomials.Rectangle)
Base.union(a::MatrixPolynomials.Rectangle, b::MatrixPolynomials.Rectangle)
```
