# * Scalar Newton polynomial

@doc raw"""
    NewtonPolynomial(ζ, d)

The unique interpolation polynomial of a function in its Newton form,
i.e.
```math
f(z) \approx p(z) = \sum_{j=1}^m \divdiff(\zeta_{1:j})f \prod_{i=1}^{j-1}(z - \zeta_i),
```
where `ζ` are the interpolation points and
`d[j]=⏃(ζ[1:j])f` is the ``j``th divided difference of the
interpolated function `f`.
"""
struct NewtonPolynomial{T,ZT<:AbstractVector{T},DT<:AbstractVector{T}}
    "Interpolation points of the Newton polynomial"
    ζ::ZT
    "Divided differences for the function interpolated by the Newton polynomial"
    d::DT
end

"""
    NewtonPolynomial(f, ζ)

Construct the Newton polynomial interpolating `f` at `ζ`,
automatically deriving the divided differences using [`⏃`](@ref).
"""
NewtonPolynomial(f::Function, ζ::AbstractVector) =
    NewtonPolynomial(ζ, ⏃(f, ζ, 1, 0, 1))

Base.view(np::NewtonPolynomial, args...) =
    NewtonPolynomial(view(np.ζ, args...), view(np.d, args...))

function (np::NewtonPolynomial{T})(z::Number; errors=nothing) where T
    update_error = zero(T)
    p = np.d[1]
    r = z - np.ζ[1]
    m = length(np.ζ)
    for i = 2:m
        p += np.d[i]*r
        update_error = abs(np.d[i])*norm(r)
        r *= z - np.ζ[i]
    end
    if errors isa AbstractVector && length(errors) ≥ 1
        errors[1] = update_error
    end
    p
end

function Base.show(io::IO, np::NewtonPolynomial)
    degree = length(np.ζ)-1
    ar,br = extrema(real(np.ζ))
    ai,bi = extrema(imag(np.ζ))
    compl_str(r,i) = if iszero(i)
        r
    elseif iszero(r) && iszero(i)
        0
    else
        r + im*i
    end

    write(io, "Newton polynomial of degree $(degree) on $(compl_str(ar,ai))..$(compl_str(br,bi))")
end

# * Newton matrix polynomial

"""
    NewtonMatrixPolynomial(np, pv, r, Ar, error)

This structure aids in the computation of the action of a matrix
polynomial on a vector. `np` is the [`NewtonPolynomial`](@ref), `pv`
is the desired result, `r` and `Ar` are recurrence vectors, and
`error` is an optional error estimator algorithm that can be used to
terminate the iterations early.
"""
struct NewtonMatrixPolynomial{T,NP<:NewtonPolynomial{T},Vec,ErrorEstim}
    np::NP
    "Action of the Newton polynomial `np` on a vector `v`"
    pv::Vec
    "Recurrence vector"
    r::Vec
    "Matrix–Recurrence vector product"
    Ar::Vec
    error::ErrorEstim
end

function NewtonMatrixPolynomial(np::NewtonPolynomial{T}, n::Integer, res=nothing) where T
    pv = Vector{T}(undef, n)
    r = Vector{T}(undef, n)
    Ar = Vector{T}(undef, n)
    NewtonMatrixPolynomial(np, pv, r, Ar, res)
end

estimate_converged!(::Nothing, args...) = false

"""
    mul!(w, p::NewtonMatrixPolynomial, A, v)

Compute the action of the [`NewtonMatrixPolynomial`](@ref) `p`
evaluated for the matrix (or linear operator) `A` acting on `v` and
storing the result in `w`, i.e. `w ← p(A)*v`.
"""
function LinearAlgebra.mul!(w, nmp::NewtonMatrixPolynomial, A, v)
    # Equation numbers refer to
    #
    # - Kandolf, P., Ostermann, A., & Rainer, S. (2014). A residual based
    #   error estimate for leja interpolation of matrix functions. Linear
    #   Algebra and its Applications, 456(nil),
    #   157–173. http://dx.doi.org/10.1016/j.laa.2014.04.023

    @unpack pv,r,Ar = nmp
    @unpack d,ζ = nmp.np

    pv .= d[1]*v # Equation (3c)
    r .= v # r is initialized using the normal iteration, below
    m = length(ζ)
    for i = 2:m
        # Equations (3a,b) are applied in reverse order, since at the
        # beginning of each iteration, r is actually lagging one
        # iteration behind, because r is initialized to v, not
        # (A-ζ[1])*v.

        # Equation (3b)
        mul!(Ar, A, r)
        lmul!(-ζ[i-1], r)
        r .+= Ar

        # Equation (3a)
        BLAS.axpy!(d[i], r, pv)

        estimate_converged!(nmp.error, A, pv, v, Ar, i-1) && break
    end

    w .= pv
end

# ** Newton matrix polynomial derivative

struct NewtonMatrixPolynomialDerivative{T,NP<:NewtonPolynomial{T},Vec}
    np::NP
    "Action of the time-derivative of the Newton polynomial `np` on a vector `v`"
    p′v::Vec
    "Recurrence vector"
    r′::Vec
    "Matrix–Recurrence vector product"
    Ar′::Vec
end

function NewtonMatrixPolynomialDerivative(np::NewtonPolynomial{T}, n::Integer) where T
    p′v = Vector{T}(undef, n)
    r′ = Vector{T}(undef, n)
    Ar′ = Vector{T}(undef, n)
    NewtonMatrixPolynomialDerivative(np, p′v, r′, Ar′)
end

function step!(nmpd::NewtonMatrixPolynomialDerivative, A, Ar, i)
    @unpack p′v,r′,Ar′ = nmpd
    @unpack d,ζ = nmpd.np

    if i == 1
        # Equation (18c)
        p′v .= false
        copyto!(r′, Ar)
    end

    # Equation (18a)
    BLAS.axpy!(d[i], r′, p′v)

    # Equation (18b)
    mul!(Ar′, A, r′)
    lmul!(-ζ[i], r′)
    r′ .+= Ar′

    p′v
end

# ** Residual error estimator for φₖ functions

"""
    φₖResidualEstimator{T,k}(nmpd, ρ, vscaled, estimate, tol)

An implementation of the residual error estimate of the φₖ functions,
as presented in

- Kandolf, P., Ostermann, A., & Rainer, S. (2014). A residual based
  error estimate for Leja interpolation of matrix functions. Linear
  Algebra and its Applications, 456(nil), 157–173. [DOI:
  10.1016/j.laa.2014.04.023](http://dx.doi.org/10.1016/j.laa.2014.04.023)

`nmpd` is a [`NewtonMatrixPolynomialDerivative`](@ref) that
successively computes the time-derivative of the
[`NewtonMatrixPolynomial`](@ref) used to interpolate ``\\varphi_k(tA)``
(the time-step ``t`` is subsequently set to unity), `ρ` is the
residual vector, `vscaled` an auxiliary vector for `k>0`, and
`estimate` and `tol` are the estimated error and tolerance,
respectively.
"""
mutable struct φₖResidualEstimator{T,k,NMPD<:NewtonMatrixPolynomialDerivative{T},Vec,R}
    nmpd::NMPD
    "Residual vector"
    ρ::Vec
    "``v/(k-1)!`` cache"
    vscaled::Vec

    estimate::R
    tol::R
end

φₖResidualEstimator(k, nmpd::NMPD, ρ::Vec, vscaled::Vec, estimate::R, tol::R) where {T,NMPD<:NewtonMatrixPolynomialDerivative{T},Vec,R} =
    φₖResidualEstimator{T,k,NMPD,Vec,R}(nmpd, ρ, vscaled, estimate, tol)

function φₖResidualEstimator(k::Integer, np::NewtonPolynomial{T}, n::Integer, tol::R) where {T,R<:AbstractFloat}
    nmpd = NewtonMatrixPolynomialDerivative(np, n)
    ρ = Vector{T}(undef, n)
    vscaled = Vector{T}(undef, k > 0 ? n : 0)
    φₖResidualEstimator(k, nmpd, ρ, vscaled, convert(R, Inf), tol)
end

function estimate_converged!(error::φₖResidualEstimator{T,k}, A, pv, v, Ar, m) where {T,k}
    @unpack ρ, vscaled = error
    if k > 0 && m == 1
        vscaled .= v/Γ(k)
    end

    mul!(ρ, A, pv)
    ρ .-= step!(error.nmpd, A, Ar, m)

    # # TODO: Figure out why this does not work as intended.
    # if k > 0
    #     ρ .+= vscaled
    #     ρ .-= k*pv
    #     @. ρ += vscaled - k*pv
    # end

    error.estimate = norm(ρ)

    if k > 0
        error.estimate /= k
    end

    error.estimate < error.tol
end

error_estimator(::typeof(exp), args...) = φₖResidualEstimator(0, args...)
error_estimator(::typeof(φ₁), args...) = φₖResidualEstimator(1, args...)
error_estimator(fix::Base.Fix1{typeof(φ),<:Integer}, args...) = φₖResidualEstimator(fix.x, args...)

export error_estimator
