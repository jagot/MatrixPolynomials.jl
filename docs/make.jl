using Documenter, MatrixPolynomials

makedocs(;
    modules=[MatrixPolynomials],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jagot/MatrixPolynomials.jl/blob/{commit}{path}#L{line}",
    sitename="MatrixPolynomials.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
)

deploydocs(;
    repo="github.com/jagot/MatrixPolynomials.jl",
)
