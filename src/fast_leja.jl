# This is a Julia port of the Matlab algorithm published in
#
# - Baglama, J., Calvetti, D., & Reichel, L. (1998). Fast Leja
#   points. Electron. Trans. Numer. Anal, 7(124-140), 119–120.

struct FastLeja{T}
    ζ::Vector{T}
    ∏ζ::Vector{T}
    ζs::Vector{T}
    ia::Vector{Int}
    ib::Vector{Int}
end

meanζ(ζ) = (i,j) -> (ζ[i]+ζ[j])/2

function fast_leja!(fl::FastLeja, n)
    @unpack ζ, ∏ζ, ζs, ia, ib = fl
    mζ = meanζ(ζ)

    curn = length(ζ)
    if curn < n
        resize!(ζ, n)
        resize!(∏ζ, n)
        resize!(ζs, n)
        resize!(ia, n)
        resize!(ib, n)
    end

    for i = curn+1:n
        maxi = argmax(abs.(view(∏ζ, 1:i-2)))
        ζ[i] = ζs[maxi]

        ia[i-1] = i
        ib[i-1] = ib[maxi]
        ib[maxi] = i

        ζs[maxi] = mζ(ia[maxi], ib[maxi])
        ζs[i-1] = mζ(ia[i-1], ib[i-1])

        sel = 1:i-1
        ∏ζ[maxi] = prod(ζs[maxi] .- ζ[sel])
        ∏ζ[i-1] = prod(ζs[i-1] .- ζ[sel])
        ∏ζ[sel] .*= ζs[sel] .- ζ[i]
    end

    fl
end

function FastLeja(a, b, n)
    T = float(promote_type(typeof(a),typeof(b)))
    ζ = zeros(T, 3)
    ∏ζ = zeros(T, 3)
    ζs = zeros(T, 3)
    ia = zeros(Int, 3)
    ib = zeros(Int, 3)

    ζ[1:2] = abs(a) > abs(b) ? [a,b] : [b,a]
    ζ[3] = (a+b)/2

    mζ = meanζ(ζ)

    ζs[1] = mζ(2,3)
    ζs[2] = mζ(3,1)

    ∏ζ[1] = prod(ζs[1] .- ζ[1:3])
    ∏ζ[2] = prod(ζs[2] .- ζ[1:3])

    ia[1] = 2
    ib[1] = 3
    ia[2] = 3
    ib[2] = 1

    fl = FastLeja(ζ, ∏ζ, ζs, ia, ib)
    fast_leja!(fl, n)
end

points(fl::FastLeja) = fl.ζ
