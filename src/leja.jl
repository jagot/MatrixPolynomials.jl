# This is an implementation of the algorithm described in
#
# - Reichel, L. (1990). Newton Interpolation At Leja Points. BIT, 30(2),
#   332–346. http://dx.doi.org/10.1007/bf02017352

struct Leja{T}
    S::Vector{T}
    ζ::Vector{T}
    ∏ζ::Vector{T}
end

# https://github.com/JuliaLang/julia/issues/27639#issuecomment-399518134
function Base.argmax(f, itr)
    r = iterate(itr)
    r === nothing && error("empty collection")
    m, state = r
    f_m = f(m)
    while true
        r = iterate(itr, state)
        r === nothing && break
        x, state = r
        f_x = f(x)
        isless(f_m, f_x) || continue
        m, f_m = x, f_x
    end
    return m
end

function leja!(l::Leja, n)
    @unpack S,ζ,∏ζ = l
    curn = length(ζ)
    n-curn > length(S) &&
        throw(DimensionMismatch("Cannot generate more Leja points than underlying discretization contains"))
    if curn < n
        resize!(ζ, n)
        resize!(∏ζ, n)
    end

    if curn == 0 && n > 0
        maxi = argmax(eachindex(S))
        ζ[1] = S[maxi]
        deleteat!(S, maxi)
        ∏ζ[1] = 0
        curn += 1
    end

    for i = curn+1:n
        maxi = argmax(eachindex(S)) do j
            ζs = S[j]
            prod(ζₖ -> abs(ζs-ζₖ), view(ζ, 1:i-1))
        end
        ζ[i] = S[maxi]
        deleteat!(S, maxi)
        ∏ζ[i] = prod(j -> abs(ζ[j]-ζ[i]), 1:i-1)
    end

    l
end

function Leja(S::AbstractVector{T}, n::Integer) where T
    l = Leja(collect(S), Vector{T}(), Vector{T}())
    leja!(l, n)
end

points(l::Leja) = l.ζ
