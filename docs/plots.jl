using PythonPlot
using Jagot.plotting
plot_style("ggplot")

import MatrixPolynomials: φ₁, φ, std_div_diff, ⏃,
    Leja, FastLeja, points, NewtonPolynomial
using SpecialFunctions

function leja()
    m = 10
    a,b = -2,2
    l = Leja(range(a, stop=b, length=1000), m)
    fl = FastLeja(a, b, m)
    ζ = points(l)
    fζ = points(fl)

    cfigure("leja") do
        csubplot(211,nox=true) do
            plot(1:m, ζ, ".-", label="Leja points")
            plot(1:m, fζ, ".--", label="Fast Leja points")
            ylabel(L"\zeta_m")
            legend()
        end
        csubplot(212) do
            plot(1:m, abs.(l.∏ζ).^(1 ./ (1:m)), ".-", label="Leja points")
            plot(1:m, abs.(fl.∏ζ).^(1 ./ (1:m)), ".--", label="Fast Leja points")
            xlabel(L"m")
            ylabel(L"C(\{\zeta_{1:m}\})")
        end
    end

    savefig("docs/src/figures/leja_points.svg")
end

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

function div_differences_cancellation()
    x = range(-2, stop=2, length=100)
    ξ = collect(x)
    f = exp
    d_std = @time std_div_diff(f, ξ, 1, 0, 1)
    d_std_big = @time std_div_diff(f, big.(ξ), 1, 0, 1)
    d_auto = @time ⏃(f, ξ, 1, 0, 1)

    cfigure("div differences cancellation") do
        loglog(d_std, label="Recursive")
        loglog(Float64.(d_std_big), label="Recursive, BigFloat")
        loglog(d_auto, "--", color="black", label="Taylor series")
        xlabel(L"j")
        ylabel(L"\Delta\!\!\!|\,(\zeta_{1:j})\exp")
    end
    legend(framealpha=0.75)

    savefig("docs/src/figures/div_differences_cancellation.svg")
end

function div_differences_sine()
    μ = 10.0 # Extent of interval
    m = 40 # Number of Leja points
    ζ = points(Leja(μ*range(-1,stop=1,length=1000),m))
    d = ⏃(sin, ζ, 1, 0, 1)
    np = NewtonPolynomial(sin, ζ)

    x = range(-μ, stop=μ, length=1000)
    f_np = np.(x)
    f_exact = sin.(x)

    cfigure("div differences sine") do
        csubplot(211, nox=true) do
            plot(x, f_np, "-", label=L"$\sin(x)$ approximation")
            plot(x, f_exact, "--", label=L"\sin(x)")
            legend()
        end
        csubplot(212) do
            semilogy(x, abs.(f_np - f_exact), label="Absolute error")
            semilogy(x, abs.(f_np - f_exact)./abs.(f_exact), label="Relative error")
            xlabel(L"x")
            legend()
        end
    end

    savefig("docs/src/figures/div_differences_sine.svg")
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
mkpath("docs/src/figures")
@echo leja()
@echo φ₁_accuracy()
@echo φₖ_accuracy()
@echo div_differences_cancellation()
@echo div_differences_sine()
