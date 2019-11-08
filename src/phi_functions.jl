@doc raw"""
    φ₁(z)

Special case of [`φ`](@ref) for `k=1`, taking care to avoid numerical
rounding errors for small ``|z|``.
"""
function φ₁(z::T) where T
    if abs(z) < eps(real(T))
        one(T)
    else
        y = exp(z)
        if abs(z) ≥ one(real(T))
            (y - 1)/z
        else
            (y-1)/log(y)
        end
    end
end

@doc raw"""
    φ(k, z)

Compute the entire function ``\varphi_k(z)``, ``z\in\mathbb{C}``,
which is recursively defined as [Eq. (2.11) of
[Hochbruck2010](http://dx.doi.org/10.1017/s0962492910000048)]
```math
\varphi_{k+1}(z) \equiv \frac{\varphi_k(z)-\varphi_k(0)}{z},
```
with the base cases
```math
\varphi_{0}(z) = \exp(z), \quad
\varphi_{1}(z) = \frac{\exp(z)-1}{z},
```
and the special case
```math
\varphi_k(0) = \frac{1}{k!}.
```

This function, as the base case [`φ₁`](@ref), is implemented to avoid
rounding errors for small ``|z|``.
"""
function φ(k, z::T) where T
    if k == 0
        exp(z)
    elseif k == 1
        φ₁(z)
    else
        abs(z) < eps(real(T)) && return 1/gamma(k+1)

        # Eq. (2.11) of
        #
        # - Hochbruck, M., & Ostermann, A. (2010). Exponential Integrators. Acta
        #   Numerica, 19(nil),
        #   209–286. http://dx.doi.org/10.1017/s0962492910000048
        #
        if abs(z) > k*one(real(T))
            (φ(k-1, z) - φ(k-1, zero(T)))/z
        else
            # Horner's rule applied to the Taylor expansion of φₖ = ∑ zⁱ/(k+i)!
            # Truncating the Taylor expansion at k+1 terms seems to work well.
            n = 10k
            b = one(T)/gamma(k+n+1)
            for i = n-1:-1:0
                b = muladd(b, z, 1/gamma(k+i+1))
            end
            b
        end
    end
end

"""
    φ(k)

Return a function corresponding to `φₖ`.

# Examples

```jldoctest
julia> φ(0)
exp (generic function with 14 methods)

julia> φ(1)
φ₁ (generic function with 1 method)

julia> φ(2)
φ₂ (generic function with 1 method)

julia> φ(15)
φ₁₅ (generic function with 1 method)

julia> φ(15)(5.0 + im)
1.0931836313419128e-12 + 9.301475570434819e-14im
```
"""
function φ(k::Integer)
    if k == 0
        exp
    elseif k == 1
        φ₁
    else
        Base.Fix1(φ, k)
    end
end

Base.string(f::Base.Fix1{typeof(φ),<:Integer}) = "φ$(to_subscript(f.x))"

function Base.show(io::IO, f::Base.Fix1{typeof(φ),<:Integer})
    write(io, "φ")
    write(io, to_subscript(f.x))
    n = length(methods(f))
    write(io, " (generic function with $n method$(n > 1 ? "s" : ""))")
end

Base.show(io::IO, ::MIME"text/plain", f::Base.Fix1{typeof(φ),<:Integer}) =
    show(io, f)
