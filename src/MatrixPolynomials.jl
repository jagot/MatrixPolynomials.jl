module MatrixPolynomials

using Parameters
using LinearAlgebra
using ArnoldiMethod
using SpecialFunctions
const Γ = gamma
using Statistics
using UnicodeFun

include("spectral_shapes.jl")
include("spectral_ranges.jl")

include("leja.jl")
include("fast_leja.jl")

include("phi_functions.jl")
include("divided_differences.jl")
include("newton.jl")

end # module
