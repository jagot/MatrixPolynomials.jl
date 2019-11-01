using MatrixPolynomials
using Test

using LinearAlgebra
import MatrixPolynomials: Line, Rectangle, spectral_range,
    Leja, leja!, FastLeja, fast_leja!, points,
    NewtonPolynomial, NewtonMatrixPolynomial, φ₁

include("spectral.jl")
include("leja.jl")
include("newton.jl")
