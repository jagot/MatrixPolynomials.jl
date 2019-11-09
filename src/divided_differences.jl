"""
    std_div_diff(f, ξ, h, c, γ)

Compute the divided differences of `f` at `h*(c .+ γ*ξ)`, where `ξ` is
a vector of (possibly complex) interpolation points, using the
standard recursion formula.
"""
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

"""
    ⏃(f, ξ, args...)

Compute the divided differences of `f` at `ξ`, using a method that is
optimized for the function `f`, if one is available, otherwise
fallback to [`MatrixPolynomials.std_div_diff`](@ref).
"""
⏃(args...) = std_div_diff(args...)

# * Special cases for φₖ(z)

"""
    φₖ_ts_div_diff(k, ξ, h, c, γ[, τ=1])

Compute the divided differences of `φₖ` at `τ*h*(c .+ γ*ξ)`, where `ξ`
is a vector of (possibly complex) interpolation points, using the
algorithm in Table 2 of

- Caliari, M. (2007). Accurate evaluation of divided differences for
  polynomial interpolation of exponential propagators. Computing,
  80(2), 189–201. [DOI: 10.1007/s00607-007-0227-1](http://dx.doi.org/10.1007/s00607-007-0227-1)
"""
function φₖ_ts_div_diff(k, ξ::AbstractVector{T}, h, c, γ, τ=one(T)) where T
    τxmax = maximum(ξᵢ -> abs(τ*h*(c + γ*ξᵢ)), ξ)
    τxmax > 1.59 && @warn("Taylor series only valid for |τ*max(x)| < 1.59, got $(τxmax)")
    m = length(ξ)
    F = Matrix{T}(undef, m, m)
    for i = 1:m
        for j = 1:i
            F[i,j] = F[j,i] = (τ*h*γ)^i/Γ(i+k-j+1)
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

⏃(::typeof(exp), args...) = φₖ_ts_div_diff(0, args...)
⏃(::typeof(φ₁), args...) = φₖ_ts_div_diff(1, args...)
⏃(fix::Base.Fix1{typeof(φ),<:Integer}, args...) = φₖ_ts_div_diff(fix.x, args...)

# These are linear fits that are always above the values of Table 3.1 of
#
# - Al-Mohy, A. H., & Higham, N. J. (2011). Computing the action of the
#   matrix exponential, with an application to exponential
#   integrators. SIAM Journal on Scientific Computing, 33(2),
#   488–511. http://dx.doi.org/10.1137/100788860

"""
    min_degree(::typeof(exp), θ)

Minimum degree of Taylor polynomial to represent `exp` to machine
precision, within a circle of radius `θ`.
"""
min_degree(::typeof(exp), θ::Float64) =
    ceil(Int, 4.1666θ + 15.0)

min_degree(::typeof(exp), θ::Float32) =
    ceil(Int, 3.7037θ + 6.8519)

"""
    taylor_series(::Type{T}, ::typeof(exp), n; s=1, θ=3.5) where T

Compute the Taylor series of `exp(z/s)`, with `n` terms, or as many
terms as required to achieve convergence within a circle of radius
`θ`, whichever is largest.
"""
function taylor_series(::Type{T}, ::typeof(exp), n; s=1, θ=3.5) where T
    N = max(n, min_degree(exp, θ))
    vcat(one(T), one(T) ./ [Γ(s*k+1) for k = 1:N])
end

"""
    div_diff_table(f, ζ[; kwargs...])

Construct the table of divided differences of `f` at the interpolation
points `ζ`, based on the algorithm on page 26 of

- Zivcovich, F. (2019). Fast and accurate computation of divided
  differences for analytic functions, with an application to the
  exponential function. Dolomites Research Notes on Approximation,
  12(1), 28–42.
"""
function div_diff_table(f, ζ::AbstractVector{T}; s=1, kwargs...) where T
    n = length(ζ)-1
    ts = taylor_series(T, f, n+1; s=s, kwargs...)
    N = length(ts)-1

    F = zeros(T, n+1, n+1)
    for i = 1:n
        F[i+1:n+1,i] .= ζ[i] .- ζ[i+1:n+1]
    end

    for j = n:-1:0
        for k = N:-1:(n-j+1)
            ts[k] += ζ[j+1]*ts[k+1]
        end
        for k = (n-j):-1:1
            ts[k] += F[k+j+1,j+1]*ts[k+1]
        end
        F[j+1,j+1:n+1] .= ts[1:n-j+1]
    end
    F[1:n+2:(n+1)^2] .= f.(ζ/s)

    UpperTriangular(F)
end

"""
    φₖ_div_diff(k, ξ[; θ=3.5, s=1])

Specialized interface to [`div_diff_table`](@ref) for the `φₖ`
functions. `θ` is the desired radius of convergence of the Taylor
series of `φₖ`, and `s` is the scaling-and-squaring parameter, which
if set to zero, will be calculated to fulfill `θ`.
"""
function φₖ_div_diff(k, ξ::AbstractVector{T}; θ=real(T(3.5)), s=1) where T
    μ = mean(ξ)
    z = vcat(zeros(k), ξ) .- μ
    n = length(z) - 1

    # Scaling
    if s == 0
        Δz = maximum(a -> maximum(b -> abs(a-b), z), z)
        s = max(1, ceil(Int, Δz/θ))
    end

    # The Taylor series of φₖ is just a shifted version of exp.
    F = div_diff_table(exp, z; s=s, θ=θ)

    dd = F[1,:]

    # Squaring
    for j = 1:s-1
        lmul!(F', dd)
    end

    exp(μ)*dd[k+1:end]
end
