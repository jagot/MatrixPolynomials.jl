using MatrixPolynomials
using Test

using LinearAlgebra
import MatrixPolynomials: Line, Rectangle, spectral_range,
    Leja, leja!, FastLeja, fast_leja!, points,
    φ₁, φ, Γ,
    std_div_diff, ts_div_diff_table, φₖ_div_diff_basis_change, ⏃,
    NewtonPolynomial, NewtonMatrixPolynomial,
    FuncV

include("spectral.jl")
include("leja.jl")
include("divided_differences.jl")
include("newton.jl")

include("funcv.jl")
