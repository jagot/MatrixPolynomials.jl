@doc raw"""
    NewtonPolynomial(ζ, d)

The unique interpolation polynomial of a function in its Newton form,
i.e.
```math
f(z) \approx p(z) = \sum_{j=1}^m f[\zeta_1,...,\zeta_j] \prod_{i=1}^{j-1}(z - \zeta_i),
```
where `ζ` are the interpolation points and
`d[j]=f[ζ₁,...,ζⱼ]` is the ``j``th divided difference of the
interpolated function `f`.
"""
struct NewtonPolynomial{T,VT<:AbstractVector{T}}
    "Interpolation points of the Newton polynomial"
    ζ::VT
    "Divided differences for the function interpolated by the Newton polynomial"
    d::VT
end

NewtonPolynomial(f::Function, ζ::AbstractVector) =
    NewtonPolynomial(ζ, std_div_diff(f, ζ, 1, 0, 1))

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

struct NewtonMatrixPolynomial{T,NP<:NewtonPolynomial{T},Vec,ErrorEstim}
    np::NP
    "Action of the Newton polynomial `p` on a vector `v`"
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

function LinearAlgebra.mul!(w, nmp::NewtonMatrixPolynomial, A, v)
    @unpack pv,r,Ar = nmp
    @unpack d,ζ = nmp.np

    pv .= d[1]*v
    r .= v
    m = length(ζ)
    for i = 2:m
        mul!(Ar, A, r)
        lmul!(-ζ[i-1], r)
        r .+= Ar

        BLAS.axpy!(d[i], r, pv)
    end

    w .= pv
end
