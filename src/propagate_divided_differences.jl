num_steps(f, xmax) = 1
num_steps(::Union{typeof(exp),typeof(φ₁)# ,Base.Fix1{typeof(φ),<:Integer}
                  }, xmax) =
    max(1, ceil(Int, xmax/0.3)) # This number is more conservative
                                # that actually is necessary, however,
                                # since ts_div_diff_table currently is
                                # implemented via repeated matrix
                                # powers, truncation of very small
                                # numbers occur. With a proper routine
                                # for powers of Bidiagonal matrices
                                # (c.f. McCurdy 1984), this value can
                                # be increased.

function num_steps(::Union{typeof(sin),typeof(cos)}, xmax)
    J = 1
    while xmax/J > 1.59
        J = nextpow(2, J+1)
    end
    J
end

"""
    propagate_div_diff(::typeof(exp), expτH, J, args...)

Find the divided differences of `exp` by utilizing that
``\\exp(a+b)=\\exp(a)\\exp(b)``.
"""
function propagate_div_diff(::typeof(exp), expτH, J, args...)
    d = expτH[:,1]
    for j = 1:J-1
        lmul!(expτH, d)
    end
    d
end

@doc raw"""
    propagate_div_diff(::typeof(φ₁), φ₁H, J, H, τ)

Find the divided differences of `φ₁` by solving the ODE
```math
\dot{\vec{y}}(t) = \mat{H} \vec{y}(t) + \vec{e}_1, \quad \vec{y}(0) = 0,
```
by iterating
```math
\vec{y}_{j+1} = \vec{y}_j + \tau\varphi_1(\tau\mat{H})(\mat{H}\vec{y}_j + \vec{e}_1),
\quad j=0,...,J-1.
```
"""
function propagate_div_diff(::typeof(φ₁), φ₁H, J, H, τ)
    d = φ₁H[:,1]
    tmp = similar(d)
    lmul!(τ, d)
    for j = 1:J-1
        mul!(tmp, H, d)
        tmp[1] += 1
        d .+= lmul!(τ, lmul!(φ₁H, tmp))
    end
    d
end

@doc raw"""
    propagate_div_diff_sin_cos(sinH, cosH, J)

Find the divided differences tables of `sin` and `cos` simultaneously,
by utilizing the double-angle formulas
```math
\sin2\theta = 2\sin\theta\cos\theta, \quad
\cos2\theta = 1 - \sin^2\theta,
```
recursively, doubling the angle at each iteration until the desired
angle is achieved.
"""
function propagate_div_diff_sin_cos(sinH, cosH, J)
    S = 2sinH*cosH
    C = I - 2sinH^2

    while J > 2
        tmp = 2S*C
        C = I - 2S^2
        S = tmp
        J >>= 1
    end

    S,C
end

"""
    propagate_div_diff(::typeof(sin), sinH, J, H, τ)

Find the divided differences of `sin`; see
[`propagate_div_diff_sin_cos`](@ref).
"""
function propagate_div_diff(::typeof(sin), sinH, J, H, τ)
    @assert ispow2(J)
    propagate_div_diff_sin_cos(sinH, taylor_series(cos)(τ*H), J)[1][:,1]
end

"""
    propagate_div_diff(::typeof(cos), cosH, J, H, τ)

Find the divided differences of `cos`; see
[`propagate_div_diff_sin_cos`](@ref).
"""
function propagate_div_diff(::typeof(cos), cosH, J, H, τ)
    @assert ispow2(J)
    propagate_div_diff_sin_cos(taylor_series(sin)(τ*H), cosH, J)[2][:,1]
end
