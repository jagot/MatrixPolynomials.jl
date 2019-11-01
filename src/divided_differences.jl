function std_div_diff(f, ξ::AbstractVector{T}, h, c, γ) where T
    m = length(ξ)
    d = Vector{T}(undef, m)
    for i = 1:m
        d[i] = f(h*(c + γ*ξ[i]))
        for j = 2:i
            d[i] = (d[i]-d[j-1])/(ξ[i]-ξ[j-1])
        end
    end
    d
end

# * Special cases for φₖ(z)

function φ₁(x::T) where T
    if abs(x) < eps(real(T))
        one(T)
    else
        y = exp(x)
        if abs(x) ≥ one(real(T))
            (y - 1)/x
        else
            (y-1)/log(y)
        end
    end
end

function φ(k, z)
    if k == 0
        exp(z)
    elseif k == 1
        φ₁(z)
    else
        iszero(z) && return 1/factorial(k)
        # Eq. (2.11) of
        #
        # - Hochbruck, M., & Ostermann, A. (2010). Exponential Integrators. Acta
        #   Numerica, 19(nil),
        #   209–286. http://dx.doi.org/10.1017/s0962492910000048
        #
        # This is probably not numerically stable for abs(z) ≈ 0
        (φ(k-1, z) - φ(k-1, 0))/z
    end
end

# Algorithm in Table 2 of
#
# - Caliari, M. (2007). Accurate evaluation of divided differences for
#   polynomial interpolation of exponential propagators. Computing,
#   80(2), 189–201. http://dx.doi.org/10.1007/s00607-007-0227-1
function φₖ_ts_div_diff(k, ξ::AbstractVector{T}, h, c, γ, τ) where T
    m = length(ξ)
    F = Matrix{T}(undef, m, m)
    for i = 1:m
        for j = 1:i
            F[i,j] = F[j,i] = (τ*h*γ)^i/gamma(i+k-j+1)
        end
    end
    for l = 2:17
        for j = 1:m-1
            F[j,j] *= τ*h*(c + γ*ξ[j])/(l+k-1)
            for i = j+1:m
                F[j,i] = τ*h*((c + γ*ξ[i])*F[j,i] + γ*F[j,i-1])/(l + i - j + k - 1)
                F[i,j] += F[j,i]
            end
        end
    end
    for j = 1:m
        F[j,j] = φ(k,τ*h*(c + γ*ξ[j]))
    end
    F[:,1]
end
