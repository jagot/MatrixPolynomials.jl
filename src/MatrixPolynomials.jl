module MatrixPolynomials

using Parameters
using LinearAlgebra
using ArnoldiMethod
using SpecialFunctions
const Γ = gamma
const lnΓ = loggamma
using Statistics
using UnicodeFun
using Formatting
using Compat

include("spectral_shapes.jl")
include("spectral_ranges.jl")

include("leja.jl")
include("fast_leja.jl")

include("phi_functions.jl")
include("matrix_closures.jl")
include("taylor_series.jl")

include("propagate_divided_differences.jl")
include("divided_differences.jl")
include("newton.jl")

include("funcv.jl")

end # module
