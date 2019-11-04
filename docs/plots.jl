using PyPlot
using Jagot.plotting
plot_style("ggplot")

import MatrixPolynomials: φ₁, φ
using SpecialFunctions

function φ₁_accuracy()
    φnaïve(x) = (exp(x) - 1)/x
    x = 10 .^ range(-18, stop=0, length=1000)

    cfigure("φ₁") do
        csubplot(211,nox=true) do
            semilogx(x, φnaïve.(x))
            semilogx(x, φ₁.(x), "--")
        end
        csubplot(212) do
            loglog(x, abs.(φnaïve.(x) - φ₁.(x))./abs.(φ₁.(x)))
            xlabel(L"x")
            ylabel("Relative error")
        end
    end

    savefig("docs/src/figures/phi_1_accuracy.svg")
end

function φₖ_accuracy()
    x = vcat(0,10 .^ range(-18,stop=2.5,length=1000))

    cfigure("φ") do
        for k = 100:-1:0
            loglog(x, φ.(k,x))
        end
        xlabel(L"x")
    end

    savefig("docs/src/figures/phi_k_accuracy.svg")

    function φnaïve(k,z)
        if k == 0
            exp(z)
        elseif k == 1
            (exp(z)-1)/z
        else
            (φnaïve(k-1,z) - 1/gamma(k))/z
        end
    end

    x = vcat(0,10 .^ range(-18,stop=0,length=1000))

    cfigure("φ naïve") do
        for k = 4:-1:0
            loglog(x, φnaïve.(k,x))
        end
        ylim(1e-3,10)
        xlabel(L"x")
    end

    savefig("docs/src/figures/phi_k_naive_accuracy.svg")
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
mkpath("docs/src/figures")
@echo φ₁_accuracy()
@echo φₖ_accuracy()
