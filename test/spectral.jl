@testset "Spectral shapes" begin
    @testset "Lines" begin
        l = Line(0.1, 1.0)

        @test 2l == Line(0.2, 2.0)
        @test range(l, 11) == range(0.1, stop=1.0, length=11)
        @test range(im*l, 11) == range(0.1im, stop=1.0im, length=11)

        @test l ∪ l == l
        @test l ∪ Line(0.0, 0.9) == Line(0.0, 1.0)
        @test l ∪ Line(0.0, 0.0) == Line(0.0, 1.0)
        @test Line(0.5(1+im),1+im) ∪ Line(0.0, 0.0) == Line(0.0, 1+im)
        @test Line(0.1+im, 0.1+im) ∪ Line(0.0,0.0) == Line(0,0.1+im)

        @test_throws ArgumentError l ∪ Line(0.0, 1+im)
        @test_throws ArgumentError l ∪ Line(0.1+im, 1+im)
        @test_throws ArgumentError Line(0.1+im, 1+im) ∪ Line(0.0,0.0)
    end

    @testset "Rectangles" begin
        r = Rectangle(0.0, 1.0+im)
        @test 0.5r == Rectangle(0.0, 0.5+0.5im)
        @test range(r, 11) == range(0.0, stop=1.0+im, length=11)

        @test r ∪ Rectangle(0.5*(1+im), 1.5*(1+im)) == Rectangle(0, 1.5*(1+im))
    end
end

@testset "Spectral ranges" begin
    A = Diagonal([1.0, 2])
    λ = spectral_range(A)
    @test λ isa Line
    @test λ.a ≈ 1.0
    @test λ.b ≈ 2.0

    λim = spectral_range(-im*A)
    @test λim isa Line
    @test λim.a ≈ -2.0im
    @test λim.b ≈ -1.0im

    λcomp = spectral_range(exp(-im*π/4)*A)
    @test λcomp isa Rectangle
    @test λcomp.a ≈ √2*(0.5 - im)
    @test λcomp.b ≈ √2*(1 - 0.5im)

    λt = spectral_range(-1:0.1:1, A)
    @test λt isa Line
    @test λt.a ≈ -2.0
    @test λt.b ≈ 2.0
end
