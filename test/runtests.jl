using MatrixPolynomials
using Test

using LinearAlgebra
import MatrixPolynomials: Line, Rectangle, spectral_range,
    Leja, leja!, FastLeja, fast_leja!, points,
    φ₁, φ, Γ,
    φₖ_ts_div_diff, φₖ_div_diff, ⏃,
    NewtonPolynomial, NewtonMatrixPolynomial

include("spectral.jl")
include("leja.jl")
include("divided_differences.jl")
include("newton.jl")
