"""
    FuncV{f,T}(s⁻¹Amc, c, s)

Structure for applying the action of a polynomial approximation of the
function `f` with a matrix argument acting on a vector, i.e. `w ←
p(A)*v` where `p(z) ≈ f(z)`. Various properties of `f` may be used,
such as shifting and scaling of the argument, to improve convergence
and/or accuracy. `s⁻¹Amc` is the shifted linear operator ``s⁻¹(A-c)``, `c` and `s`
are the shift and scaling, respectively.
"""
struct FuncV{f,Op,P,C,S}
    "Shifted and scaled linear operator that is iterated"
    s⁻¹Amc::Op
    "Polynomial approximation of `f`"
    p::P
    "Numeric shift employed in iterations"
    c::C
    "Scaling employed in iterations"
    s::S
    FuncV(f::Function, s⁻¹Amc::Op, p::P, c::C, s::S) where {Op,P,C,S} =
        new{f,Op,P,C,S}(s⁻¹Amc, p, c, s)
end

# For arbitrary functions, we do not scale or shift, since there are
# no universal scaling and/or shifting laws.
scaling(::Function, λ) = 1
shift(::Function, λ) = 0

scale(A::AbstractMatrix, h) = h*A
shift(A::AbstractMatrix, c) = A - I*c

"""
    FuncV(f, A, m[, t=1; distribution=:leja, kwargs...])

Create a [`FuncV`](@ref) that is used to compute `f(t*A)*v` using a
polynomial approximation of `f` formed using `m` interpolation points.

`kwargs...` are passed on to [`spectral_range`](@ref) which estimates
the range over which `f` has to be interpolated.
"""
function FuncV(f::Function, A, m::Integer, t=one(eltype(A));
               distribution=:leja, leja_multiplier=100,
               tol=1e-15, kwargs...)
    λ = spectral_range(t, A; kwargs...)

    ζ = if distribution == :leja
        points(Leja(range(λ, m*leja_multiplier), m))
    elseif distribution == :fast_leja
        points(FastLeja(λ.a, λ.b, m))
    else
        throw(ArgumentError("Invalid distribution of interpolation points $(distribution); valid choices are :leja and :fast_leja"))
    end

    s = scaling(f, λ)
    c = shift(f, λ)

    s⁻¹Amc = scale(shift(scale(A, t), c), 1/s)

    n = size(A,1)
    p = if n > 1
        d = ⏃(f, ζ, 1, 0, 1)

        np = NewtonPolynomial(ζ, d)
        NewtonMatrixPolynomial(np, n, error_estimator(f, np, n, tol))
    else
        @warn "Scalar case, no interpolation polynomial necessary"
        nothing
    end

    FuncV(f, s⁻¹Amc, p, c, s)
end

# This does not yet consider substepping
matvecs(f::FuncV) = matvecs(f.p)

unshift!(w, funcv::FuncV) = @assert iszero(funcv.c)

single_step!(w, funcv::FuncV, v) =
    mul!(w, funcv.p, funcv.s⁻¹Amc, v)

"""
    mul!(w, funcv::FuncV, v)

Evaluate the action of the matrix polynomial `funcv` on `v` and store
the result in `w`.
"""
function LinearAlgebra.mul!(w, funcv::FuncV, v)
    single_step!(w, funcv, v)
    unshift!(w, funcv)
end

scalar(A::AbstractMatrix) = A[1]
function scalar(A)
    v = ones(eltype(A),1)
    w = similar(v)
    mul!(w, A, v)[1]
end

function single_step!(w, funcv::FuncV{f,<:Any,Nothing}, v) where f
    copyto!(w, v)
    lmul!(f(scalar(funcv.s⁻¹Amc)), w)
end
