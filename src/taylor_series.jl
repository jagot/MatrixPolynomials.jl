@doc raw"""
    TaylorSeries(d, c)

Represents the Taylor series of a function as

```math
f(x) = \sum_{k=0}^\infty c_k x^{d_k},
```

where `dₖ = d(k)` and `cₖ = c(k)`.
"""
struct TaylorSeries{D,C}
    d::D
    c::C
end

function Base.show(io::IO, ts::TaylorSeries)
    for k = 0:3
        d = ts.d(k)
        c = ts.c(k)
        k > 0 && write(io, " ")
        write(io, c < 0 ? "- " : (k > 0 ? "+ " : ""))
        !isone(abs(c)) && write(io, "$(abs(c))")
        if d == 0
            isone(abs(c)) && write(io, "1")
        elseif d == 1
            write(io, "x")
        else
            write(io, "x^$(d)")
        end
    end
    write(io, " + ...")
end

function bisect_find_last(f, r::UnitRange{<:Integer})
    f(r[1]) || return nothing
    f(r[end]) && return r[end]
    while length(r) > 2
        i = div(length(r),2)
        if isodd(i)
            i = length(r) - i
        end
        fhalf = f(r[i])
        r = fhalf ? (r[i]:r[end]) : (r[1]:r[i])
    end
    r[1]
end

"""
    (ts)(x; max_degree=17)

Evaluate the Taylor series represented by `ts` up to a maximum degree
in `x` (default 17).
"""
function (ts::TaylorSeries)(x; max_degree=17)
    kmax = bisect_find_last(k -> ts.d(k) ≤ max_degree, 0:max_degree)
    v = zero(closure(x), size(x)...) + ts.c(kmax)*I
    d_prev = ts.d(kmax)

    # This is a special version of Horner's rule that takes into
    # account that some powers may be absent from the Taylor
    # polynomial.
    for k = kmax-1:-1:0
        d = ts.d(k)
        for j = 1:(d_prev-d)
            v *= x
        end
        
        d_prev = d
        v += ts.c(k)*I

        # Handle odd cases, e.g. sine
        if k == 0 && d == 1
            v *= x
        end
    end

    v
end

macro taylor_series(f, d, c)
    docstring = """
    taylor_series(::typeof($f))

Generates the [`TaylorSeries`](@ref) of `$f(x) = ∑ₖ x^($d) $c`.

# Example

```julia-repl
julia> taylor_series($f)

"""
    quote
        @doc $docstring*string(TaylorSeries(k -> $d, k -> $c))*"\n```"
        taylor_series(::typeof($f)) = TaylorSeries(k -> $d, k -> $c)
    end |> esc
end

@taylor_series exp k 1/Γ(k+1)
@taylor_series sin 2k+1 (-1)^k/Γ(2k+2)
@taylor_series cos 2k (-1)^k/Γ(2k+1)

@taylor_series sinh 2k+1 1/Γ(2k+2)
@taylor_series cosh 2k 1/Γ(2k+1)

@taylor_series φ₁ k 1/Γ(k+2)

taylor_series(fix::Base.Fix1{typeof(φ),<:Integer}) = TaylorSeries(k -> k, k -> 1/Γ(k+fix.x+1))
