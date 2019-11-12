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
    mul!(w, nmp::NewtonMatrixPolynomial, A, v)

Compute the action of the [`NewtonMatrixPolynomial`](@ref) `nmp`
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
