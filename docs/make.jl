using Documenter, MatrixPolynomials, LinearAlgebra, Statistics

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(;
    modules=[MatrixPolynomials],
    format = Documenter.HTML(assets = ["assets/latex.js"],
                             mathengine = Documenter.MathJax()),
    pages=[
        "Home" => "index.md",
        "Functions of matrices" => "funcv.md",
        "Leja points" => "leja.md",
        "Divided differences" => "divided_differences.md",
        "Newton polynomials" => "newton_polynomials.md",
        "φₖ functions" => "phi_functions.md",
    ],
    repo="https://github.com/jagot/MatrixPolynomials.jl/blob/{commit}{path}#L{line}",
    sitename="MatrixPolynomials.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
)

deploydocs(;
    repo="github.com/jagot/MatrixPolynomials.jl",
)
