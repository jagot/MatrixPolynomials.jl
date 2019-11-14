"""
    spectral_range(A[; ctol=√ε, verbosity=0])

Estimate the spectral range of the matrix/linear operator `A` using
[ArnoldiMethod.jl](https://github.com/haampie/ArnoldiMethod.jl), with
`ctol` setting the tolerance for the Krylov iterations. Returns a
spectral [`Shape`](@ref).

# Examples

```julia-repl
julia> A = Diagonal(1.0:6)
6×6 Diagonal{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 1.0   ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅   2.0   ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅   3.0   ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅   4.0   ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅   5.0   ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅   6.0

julia> MatrixPolynomials.spectral_range(A, verbosity=2)
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
[ Info: Imaginary extent of spectral range 0.0 below tolerance 1.4901161193847656e-8, conflating.
[ Info: 0.0 below tolerance 1.4901161193847656e-8, truncating.
[ Info: 0.0 below tolerance 1.4901161193847656e-8, truncating.
MatrixPolynomials.Line{Float64}(0.9999999999999998 + 0.0im, 6.0 + 0.0im)

julia> MatrixPolynomials.spectral_range(exp(im*π/4)*A, verbosity=2)
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
Converged: 1 of 1 eigenvalues in 6 matrix-vector products
MatrixPolynomials.Rectangle{Float64}(0.7071067811865468 + 0.7071067811865478im, 4.242640687119288 + 4.242640687119283im)
```
The second example should also be a [`Line`](@ref), but the algorithm
is not yet clever enough.
"""
function spectral_range(A; ctol=√(eps(real(eltype(A)))), verbosity=0, kwargs...)
    r = map([(SR(),real),(SI(),imag),(LR(),real),(LI(),imag)]) do (which,comp)
        schurQR,history = partialschur(A; which=which, nev = 1, kwargs...)
        verbosity > 0 && println(history)
        comp(first(schurQR.eigenvalues))
    end

    for (label,(i,j)) in [("Real", (1,3)),("Imaginary", (2,4))]
        if abs(r[i]-r[j]) < ctol
            verbosity > 1 && @info "$label extent of spectral range $(abs(r[i]-r[j])) below tolerance $(ctol), conflating."
            r[i] = r[j] = (r[i]+r[j])/2
        end
    end

    for i = 1:4
        if abs(r[i]) < ctol
            verbosity > 1 && @info "$(abs(r[i])) below tolerance $(ctol), truncating."
            r[i] = zero(r[i])
        end
    end
    a,b = (r[1]+im*r[2], r[3]+im*r[4])

    if real(a) == real(b) || imag(a) == imag(b)
        Line(a,b)
    else
        # There could be cases where all the eigenvalues fall on a
        # sloped line in the complex plane, but we don't know how to
        # deduce that yet. The user is free to define such sloped
        # lines manually, though.
        Rectangle(a,b)
    end
end

"""
    spectral_range(t, A)

Finds the spectral range of `t*A`; if `t` is a vector, find the
largest such range.
"""
function spectral_range(t, A; kwargs...)
    λ = spectral_range(A; kwargs...)
    ta,tb = extrema(t)
    ta*λ ∪ tb*λ
end
