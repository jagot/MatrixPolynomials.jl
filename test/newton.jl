function test_scalar_newton_leja(f, x, m, x̃, h, h̃)
    ξ = points(Leja(x, m))
    @info "$m $(eltype(ξ)) Leja points"
    np = NewtonPolynomial(f, ξ)
    y = h.(x̃,Ref(np))
    ỹ = h̃.(x̃)

    ms = 1:m
    errors = zeros(m)
    error_estimates = zeros(m,1)
    for m = ms
        np′ = view(np, 1:m)
        p = x -> begin
            val,err = np′(x, true)
            error_estimates[m,1] = max(error_estimates[m,1], err)
            val
        end
        y′ = similar(y)
        for i = eachindex(x̃)
            y′[i] = h(x̃[i],p)
        end
        errors[m] = norm(y′ - h̃.(x̃))
    end
    norm(y-ỹ),errors,error_estimates
end

function test_mat_newton_leja(f, x, m, x̃, h!, h̃)
    ξ = points(Leja(x, m))
    @info "$m $(eltype(ξ)) Leja points"
    np = NewtonPolynomial(f, ξ)
    @info "Reference solution"
    ỹ = @time reduce(hcat, h̃.(x̃))
    n = size(ỹ, 1)

    nmp = NewtonMatrixPolynomial(np, n)
    y = similar(ỹ)
    hh! = h!(nmp)
    @info "Leja/Newton solution"
    @time for i = eachindex(x̃)
        hh!(view(y,:,i), x̃[i])
    end

    ms = 1:m
    errors = zeros(m)
    for m = ms
        np′ = view(np, 1:m)
        nmp′ = NewtonMatrixPolynomial(np′, n)
        hh! = h!(nmp′)

        y′ = similar(y)
        for i = eachindex(x̃)
            hh!(view(y′,:,i), x̃[i])
        end
        errors[m] = norm(y′ - ỹ)
    end
    norm(y-ỹ),errors
end


@testset "Newton polynomials" begin
    @testset "Scalar polynomials" begin
        dx = 2
        x = dx*range(-1,stop=1,length=1000)
        @testset "Quadratic function" begin
            ξ = points(Leja(x, 10))
            f = x -> x^2
            d = std_div_diff(f, ξ, 1, 0, 1)
            np = NewtonPolynomial(ξ, d)
            @test np.d[1:3] ≈ [4,0,1]
            @test all(d -> isapprox(d, 0, atol=√(eps(d))), np.d[4:end])
        end
        @testset "Exponential function" begin
            Δy,errors,error_estimates = test_scalar_newton_leja(exp, x, 20, x, (t, p) -> p(t), exp)
            @test Δy < 7e-14
            @test all(errors[end-2:end] .< 1e-13)
        end
        @testset "Inhomogeneous ODE, $kind" for (kind,m,tol) in [(:real,43,1e-13), (:complex,60,5e-13)]
            y₀ = 1.0
            g = -3.0
            tmax = 10.0

            b,tmin = if kind == :real
                -2, 0
            else
                -2im, -tmax
            end

            t = range(tmin, stop=tmax, length=1000)

            h = (t, p) -> y₀ + t*p(t*b)*(b*y₀ + g)
            h̃ = t -> exp(t*b)*(y₀ + g/b) - g/b

            Δy,errors,error_estimates = test_scalar_newton_leja(φ₁, b*t, m, t, h, h̃)
            @test Δy < tol
            @test all(errors[end-12:end] .< 3tol)
            # Should also look at error estimates
        end
    end

    @testset "Matrix polynomials" begin
        @testset "Inhomogeneous coupled ODEs, $kind" for (kind,m,tol) in [(:real,43,1e-12), (:complex,60,1e-12)]
            n = 10 # Number of ODEs
            Y₀ = 1.0*ones(kind == :real ? Float64 : ComplexF64, n)
            G = -3*ones(n) # Inhomogeneous terms

            n_discr = 1000 # Number of points spanning eigenspectrum interval
            m = 60 # Number of Leja points

            tmax = 10.0

            b,c,tmin = if kind == :real
                -2, 0.2, 0
            else
                -2im, 0.2im, -tmax
            end

            Bdiag = Diagonal(b./(1:n))
            o = ones(n)
            B = Bdiag + Tridiagonal(c*o[2:end], 0c*o, c*o[2:end])

            t = range(tmin, stop=tmax, length=1000)

            @show λ = spectral_range(t, B, verbosity=2)

            H! = function(p)
                (w,t) -> BLAS.axpy!(1, Y₀, lmul!(t, mul!(w, p, t*B, B*Y₀ + G)))
            end
            H̃ = t -> exp(t*Matrix(B))*(Y₀ + B\G) - B\G

            Δy,errors = test_mat_newton_leja(φ₁, range(λ, n_discr), m, t, H!, H̃)
            @test errors[end] < 5e-12
        end
    end
end
