@testset "Leja points" begin
    @testset "$llabel Leja points" for (llabel,LejaType,extra_func!) = [("Discretized", (a,b,n)->Leja(range(a,stop=b,length=1001),n), leja!),
                                                                        ("Fast", (a,b,n)->FastLeja(a,b,n), fast_leja!)]
        @testset "$label Leja points" for (label,comp,factor) = [("Real", real, 1.0),
                                                                 ("Imaginary", imag, 1.0im)]
            a = -2*factor
            b = 2*factor
            n = 100
            l = LejaType(a, b, n)

            ζ = points(l)
            @test abs(ζ[1]) == 2
            @test ζ[2] == -ζ[1]
            @test ζ[3] == 0
            @test length(ζ) == n
            @test allunique(ζ)
            @test all(comp(a) .≤ comp(ζ) .≤ comp(b))

            extra_func!(l, 300)
            @test length(ζ) == 300
            @test allunique(ζ)
            @test all(comp(a) .≤ comp(ζ) .≤ comp(b))
        end
    end
end
