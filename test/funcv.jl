function test_stepping(f, x, m, μ, x̃, h̃; kwargs...)
    @info "Reference solution"
    ỹ = @time reduce(hcat, h̃.(x̃))
    n = size(ỹ,1)

    t = μ*step(x̃)
    
    f̂ = FuncV(f, x, m, t; kwargs...)

    y = similar(ỹ)
    y[:,1] = ỹ[:,1]
    
    @info "Leja/Newton solution"
    @time for i = 2:length(x̃)
        mul!(view(y,:,i), f̂, view(y,:,i-1))
    end
    
    global_error = abs.(y-ỹ)
    # This is only a rough estimate of the local error
    local_error = vcat(0, abs.(diff(global_error, dims=2)))
    
    @info "Maximum global error: $(maximum(global_error))"
    @info "Maximum local error: $(maximum(local_error))"
    
    global_error, local_error
end

function tdse(N, ρ, ℓ)
    j = 1:N
    r = (j .- 1/2)*ρ

    j² = j.^2
    α = (j²./(j² .- 1/4))[1:end-1]
    β = (j² - j .+ 1/2)./(j² - j .+ 1/4)

    T = Tridiagonal(α, -2β, α)/(-2ρ^2)

    V = Diagonal(-1 ./ r + ℓ*(ℓ + 1) ./ 2r.^2)
    
    ψ₀ = exp.(-r.^2)
    lmul!(1/√ρ, normalize!(ψ₀))

    T,V,ψ₀
end

@testset "FuncV" begin
    @testset "TDSE" begin
        N = 7
        ρ = 0.1
        L = 1
        tmax = 1.0
        t = range(0, stop=tmax, length=1000)
        
        T,V,ψ₀ = tdse(N, ρ, L)
        H = T+V
        B = -im*H

        # Exact solution
        F̃ = t -> exp(t*Matrix(B))*ψ₀
        
        m = 40 # Number of Leja points
        
        @testset "Tolerance = $tol" for (tol,exp_error) in [(3e-14,7e-12),
                                                            (1e-12,1e-9)]
            global_error, local_error = test_stepping(exp, H, m, -im, t, F̃, tol=tol)
            @test all(global_error .≤ exp_error)
            @test all(local_error .≤ 10tol)
            # Should test number of Leja points used
        end
    end
end
