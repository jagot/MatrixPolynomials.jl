"""
    Shape

Abstract base type for different shapes in the complex plane
encircling the spectra of linear operators.
"""
abstract type Shape{T} end

# * Line

"""
    Line(a, b)

For spectra falling on a line in the complex plane from `a` to `b`.
"""
struct Line{T} <: Shape{T}
    a::Complex{T}
    b::Complex{T}
    function Line(a::A, b::B) where {A,B}
        T = real(promote_type(A,B))
        new{T}(Complex{T}(a), Complex{T}(b))
    end
end

"""
    n * l::Line

Scale the [`Line`](@ref) `l` by `n`.
"""
Base.:(*)(n::Number, l::Line) = Line(n*l.a, n*l.b)
Base.iszero(l::Line) = iszero(l.a) && iszero(l.b)

"""
    range(l::Line, n)

Generate a uniform distribution of `n` values along the [`Line`](@ref)
`l`. If the line is on the real axis, the resulting values will be
real as well.

# Examples

```julia-repl
julia> range(MatrixPolynomials.Line(0, 1.0im), 5)
5-element LinRange{Complex{Float64}}:
 0.0+0.0im,0.0+0.25im,0.0+0.5im,0.0+0.75im,0.0+1.0im

julia> range(MatrixPolynomials.Line(0, 1.0), 5)
0.0:0.25:1.0
```
"""
function Base.range(l::Line, n)
    a,b = if isreal(l.a) && isreal(l.b)
        real(l.a), real(l.b)
    else
        l.a, l.b
    end
    range(a,stop=b,length=n)
end

"""
    mean(l::Line)

Return the mean value along the [`Line`](@ref) `l`.
"""
function Statistics.mean(l::Line)
    μ = (l.a + l.b)/2
    isreal(μ) ? real(μ) : μ
end

function LinearAlgebra.normalize(z::Number)
    N = norm(z)
    iszero(N) ? z : z/N
end

"""
    a::Line ∪ b::Line

Form the union of the [`Line`](@ref)s `a` and `b`, which need to be
collinear.
"""
function Base.union(a::Line, b::Line)
    a == b && return a
    da = normalize(a.b - a.a)
    db = normalize(b.b - b.a)

    # Every line is considered "parallel" with the origin
    if !iszero(a) && !iszero(b)
        cosθ = da'db
        cosθ ≈ 1 || cosθ ≈ -1 ||
            throw(ArgumentError("Lines $a and $b not parallel"))
    end

    # If the lengths of both lines are zero, then they are "collinear"
    # by definition. Otherwise, the points of one line have to lie on
    # the extension of the other line.
    if !iszero(da) || !iszero(db)
        t = (b.a - a.a)/(iszero(da) ? db : da)
        isapprox(imag(t), zero(t), atol=eps(real(t))) ||
            throw(ArgumentError("Lines $a and $b not collinear"))
    end

    ra,rb = extrema([real(a.a),real(a.b),real(b.a),real(b.b)])
    ia,ib = extrema([imag(a.a),imag(a.b),imag(b.a),imag(b.b)])
    Line(ra+im*ia, rb+im*ib)
end

# * Rectangle

"""
    Rectangle(a,b)

For spectra falling within a rectangle in the complex plane with
corners `a` and `b`.
"""
struct Rectangle{T} <: Shape{T}
    a::Complex{T}
    b::Complex{T}
    function Rectangle(a::A, b::B) where {A,B}
        T = real(promote_type(A,B))
        new{T}(Complex{T}(a), Complex{T}(b))
    end
end

"""
    n * r::Rectangle

Scale the [`Rectangle`](@ref) `r` by `n`.
"""
Base.:(*)(n::Number, r::Rectangle) = Rectangle(n*r.a, n*r.b)


"""
    range(r::Rectangle, n)

Generate a uniform distribution of `n` values along the diagonal of
the [`Rectangle`](@ref) `r`.

This assumes that the eigenvalues lie on the diagonal of the
rectangle, i.e. that the spread is negligible. It would be more
correct to instead generate samples along the sides of the rectangle,
however, [`spectral_range`](@ref) needs to be modified to correctly
identify spectral ranges falling on a line that is not lying on the
real or imaginary axis.
"""
Base.range(r::Rectangle, n) = range(r.a,stop=r.b,length=n)

"""
    mean(r::Rectangle)

Return the middle value of the [`Rectangle`](@ref) `r`.
"""
function Statistics.mean(r::Rectangle)
    μ = (r.a + r.b)/2
    isreal(μ) ? real(μ) : μ
end

"""
    a::Rectangle ∪ b::Rectangle

Find the smalling [`Rectangle`](@ref) encompassing `a` and `b`.
"""
function Base.union(a::Rectangle, b::Rectangle)
    ra,rb = extrema([real(a.a),real(a.b),real(b.a),real(b.b)])
    ia,ib = extrema([imag(a.a),imag(a.b),imag(b.a),imag(b.b)])
    Rectangle(ra+im*ia, rb+im*ib)
end
