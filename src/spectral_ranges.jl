mutable struct Lanczos{T,Op}
    A::Op
    Q::Matrix{T}
    α::Vector{T}
    β::Vector{T}
    k::Int
end

function Lanczos(A, K; kwargs...)
    T = real(eltype(A))
    Q = zeros(T, size(A,1), K)
    α = zeros(T, K)
    β = zeros(T, K)

    reset!(Lanczos(A, Q, α, β, 1); kwargs...)
end

function reset!(l::Lanczos{T}; v=nothing) where T
    if isnothing(v)
        v = rand(T, size(l.A,1))
    end
    l.Q[:,1] = normalize(v)
    l.k = 1
    l
end

function step!(l::Lanczos; verbosity=0)
    @unpack A,Q,α,β,k = l
    k == size(Q,2) && return

    v = view(Q, :, k)
    w = view(Q, :, k+1)

    mul!(w, A, v)
    α[k] = dot(v, w)
    w .-= α[k] .* v
    if k > 1
        w .-= β[k-1] .* view(Q, :, k-1)
    end

    β[k] = norm(w)
    w ./= β[k]

    verbosity > 0 && printfmtln("iter {1:d}, α[{1:d}] {2:e}, β[{1:d}] {3:e}", k, α[k], β[k])
    l.k += 1
end

LinearAlgebra.SymTridiagonal(L::Lanczos) =
    SymTridiagonal(view(L.α, 1:L.k-1), view(L.β, 1:L.k-2))

"""
    hermitian_spectral_range(A;[ K=20])

Estimate the spectral range of a Hermitian operator `A` using Algorithm 1 of

- Zhou, Y., & Li, R. (2011). Bounding the spectrum of large hermitian
  matrices. Linear Algebra and its Applications, 435(3),
  480–493. http://dx.doi.org/10.1016/j.laa.2010.06.034

"""
function hermitian_spectral_range(A; K=min(20,size(A,1)-1), ctol=√(eps(real(eltype(A)))),
                                  verbosity=0, kwargs...)
    verbosity > 0 &&
        @info "Trying to estimate spectral interval for allegedly Hermitian operator"

    l = Lanczos(A, K+1; kwargs...)
    K̃ = min(4,K-1)
    for k = 1:K̃
        step!(l; verbosity=verbosity-1)
    end

    U = real(eltype(A))
    βₖ = zero(U)
    λₘᵢₙₖ = zero(U)
    λₘₐₓₖ = zero(U)
    zₘᵢₙₖ = zero(U)
    zₘₐₓₖ = zero(U)

    for k = K̃+1:K
        step!(l; verbosity=verbosity-1)
        ee = eigen(SymTridiagonal(l))
        Z = ee.vectors

        βₖ = l.β[k]
        zₘᵢₙₖ = abs(Z[k,1])
        zₘₐₓₖ = abs(Z[k,k])

        λₘᵢₙₖ = ee.values[1]
        λₘₐₓₖ = ee.values[k]

        verbosity > 0 && printfmtln("k = {1:3d} β[k] = {2:e} λₘᵢₙ(Tₖ) = {3:+e} λₘₐₓ(Tₖ) = {4:+e} zₘᵢₙ[k]β[k] = {5:e} zₘₐₓ[k]β[k] = {6:e}",
                                    k, βₖ, λₘᵢₙₖ, λₘₐₓₖ, zₘᵢₙₖ*βₖ, zₘₐₓₖ*βₖ)
        if zₘₐₓₖ*βₖ < ctol
            # Zhou et al. (2011), bound (2.8)
            return Line(λₘᵢₙₖ - min(zₘᵢₙₖ, abs(Z[k,2]), abs(Z[k,3]))*βₖ,
                        λₘₐₓₖ + max(zₘₐₓₖ, abs(Z[k,k-1]), abs(Z[k,k-2]))*βₖ)
        end
    end
    # Zhou et al. (2011), mean of bounds (2.5,6)
    Line(λₘᵢₙₖ - (1+zₘᵢₙₖ)*βₖ/2,
         λₘₐₓₖ + (1+zₘₐₓₖ)*βₖ/2)
end

"""
    spectral_range(A[; ctol=√ε, verbosity=0])

Estimate the spectral range of the matrix/linear operator `A` using
[ArnoldiMethod.jl](https://github.com/haampie/ArnoldiMethod.jl). If
the spectral range along the real/imaginary axis is smaller than
`ctol`, it is compressed into a line. Returns a spectral
[`Shape`](@ref).

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
function spectral_range(A; ctol=√(eps(real(eltype(A)))), ishermitian=false, verbosity=0, kwargs...)
    ishermitian && return hermitian_spectral_range(A; ctol=ctol, verbosity=verbosity, kwargs...)
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

function spectral_range(A::SymTridiagonal{<:Real}; kwargs...)
    n = size(A,1)
    a,b = minmax(first(eigvals(A, 1:1)),
                 first(eigvals(A, n:n)))
    Line(a,b)
end

spectral_range(A::Diagonal{<:Real}; kwargs...) =
    Line(minimum(A.diag), maximum(A.diag))

function spectral_range(A::Diagonal{<:Complex}; kwargs...)
    d = A.diag
    rd = real(d)
    id = imag(d)
    a = minimum(rd) + im*minimum(id)
    b = maximum(rd) + im*maximum(id)
    if real(a) == real(b) || imag(a) == imag(b)
        Line(a,b)
    else
        Rectangle(a, b)
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
