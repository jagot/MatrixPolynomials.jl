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
