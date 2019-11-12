"""
    closure(x::Number)

Generates the closure type of `xⁿ` as `n → ∞`, i.e. a scalar.
"""
closure(::T) where {T<:Number} = T

for Mat = [:Matrix, :Diagonal, :LowerTriangular, :UpperTriangular]
    docstring = """
    closure(x::$Mat)

Generates the closure type of `xⁿ` as `n → ∞`, i.e. a `$Mat`.
"""
    @eval begin
        @doc $docstring
        closure(::M) where {M<:$Mat} = M
    end
end

for Mat = [:Tridiagonal, :SymTridiagonal]
    docstring = """
    closure(x::$Mat)

Generates the closure type of `xⁿ` as `n → ∞`, i.e. a `Matrix`.
"""
    @eval closure(::$Mat{T}) where T = Matrix{T}
end

"""
    closure(x::Bidiagonal)

Generates the closure type of `xⁿ` as `n → ∞`, i.e. a
`UpperTriangular` or `LowerTriangular`, depending on `x.uplo`.
"""
function closure(B::Bidiagonal{T}) where T
    if B.uplo == 'L'
        LowerTriangular{T,Matrix{T}}
    else
        UpperTriangular{T,Matrix{T}}
    end
end

function Base.zero(::Type{Mat}, m, n) where {T,Mat<:Diagonal{T}}
    @assert m == n
    Mat(zeros(m))
end

function Base.zero(::Type{Mat}, m, n) where {T,Mat<:Tridiagonal{T}}
    @assert m == n
    Mat(zeros(T,m-1),zeros(T,m),zeros(T,m-1))
end

function Base.zero(::Type{Mat}, m, n) where {T,Mat<:SymTridiagonal{T}}
    @assert m == n
    Mat(zeros(T,m),zeros(T,m-1))
end

Base.zero(::Type{Mat}, m, n) where {T,Mat<:AbstractMatrix{T}} = Mat(zeros(T, m, n))
