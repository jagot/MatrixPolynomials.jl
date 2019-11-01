module MatrixPolynomials

using Parameters
using LinearAlgebra
using ArnoldiMethod

include("spectral_shapes.jl")
include("spectral_ranges.jl")

include("leja.jl")
include("fast_leja.jl")
end # module
