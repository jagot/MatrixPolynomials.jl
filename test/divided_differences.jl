@testset "Divided differences" begin
    @testset "Infinitesimal divided differences" begin
        # This will lead to catastrophic cancellation for the standard
        # recursive formulation of divided differences; the goal is to
        # ascertain the accuracy of the optimized methods, which
        # should converge to the Taylor expansion of φₖ(z) for
        # infinitesimal |z|.
        @testset "$label" for (label,μ) in [("Real", 1.0),
                                            ("Imaginary", 1.0im),
                                            ("Complex", exp(im*π/4))]
            m = 100
            ξ = μ*eps(Float64)*range(-1,stop=1,length=m)
            x = μ*eps(Float64)*range(-1,stop=1,length=1000)
            @testset "k = $k" for k = 0:6
                f = φ(k)
                f_exact = f.(x)
                taylor_expansion = vcat(1.0 ./ [Γ(k+n+1) for n = 0:m-1])
                @testset "$method" for (method, func) in [
                    ("Taylor series", (k,ξ) -> φₖ_ts_div_diff(k, ξ, 1, 0, 1)),
                    ("Basis change", (k,ξ) -> φₖ_div_diff(k, ξ)),
                    ("Auto", (k,ξ) -> ⏃(f, ξ, 1, 0, 1))
                ]
                    d = func(k, ξ)
                    @test d ≈ taylor_expansion atol=1e-15
                    # The Taylor expansion coefficients for φₖ(z)
                    # should all be real, to machine precision, even
                    # for complex z.
                    @test norm(imag(d)) ≈ 0 atol=1e-15

                    f_d = NewtonPolynomial(ξ, d).(x)
                    @test f_d ≈ f_exact atol=1e-14
                end
            end
        end
    end

    @testset "Finite divided differences" begin
        @testset "$label" for (label,μ) in [("Real", 1.0),
                                            ("Imaginary", 1.0im),
                                            ("Complex", exp(im*π/4))]
            m = 100
            dx = 1.0
            ξ = μ*dx*range(-1,stop=1,length=m)
            x = μ*dx*range(-1,stop=1,length=1000)
            @testset "k = $k" for k = 0:6
                f = φ(k)
                f_exact = f.(x)
                @testset "$method" for (method, func) in [
                    ("Taylor series", (k,ξ) -> φₖ_ts_div_diff(k, ξ, 1, 0, 1)),
                    ("Basis change", (k,ξ) -> φₖ_div_diff(k, ξ)),
                    ("Auto", (k,ξ) -> ⏃(f, ξ, 1, 0, 1))
                ]
                    d = func(k, ξ)
                    f_d = NewtonPolynomial(ξ, d).(x)
                    @test f_d ≈ f_exact atol=1e-12
                end
            end
        end
    end
end
