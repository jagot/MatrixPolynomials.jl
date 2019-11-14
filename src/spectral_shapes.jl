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

Base.:(*)(n::Number, l::Line) = Line(n*l.a, n*l.b)
Base.iszero(l::Line) = iszero(l.a) && iszero(l.b)

function Base.range(l::Line, n)
    a,b = if isreal(l.a) && isreal(l.b)
        real(l.a), real(l.b)
    else
        l.a, l.b
    end
    range(a,stop=b,length=n)
end

function LinearAlgebra.normalize(z::Number)
    N = norm(z)
    iszero(N) ? z : z/N
end

function Base.union(a::Line, b::Line)
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

Base.:(*)(n::Number, r::Rectangle) = Rectangle(n*r.a, n*r.b)

# This assumes that the eigenvalues lie on the diagonal of the
# rectangle, i.e. that the spread is negligible. It may be more useful
# in the future to instead generate samples along the sides of the
# rectangle.
Base.range(r::Rectangle, n) = range(r.a,stop=r.b,length=n)

function Base.union(a::Rectangle, b::Rectangle)
    ra,rb = extrema([real(a.a),real(a.b),real(b.a),real(b.b)])
    ia,ib = extrema([imag(a.a),imag(a.b),imag(b.a),imag(b.b)])
    Rectangle(ra+im*ia, rb+im*ib)
end
