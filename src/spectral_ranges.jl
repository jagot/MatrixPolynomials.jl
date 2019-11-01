function spectral_range(A; ctol=√(eps(real(eltype(A)))), verbosity=0, kwargs...)
    r = map([(SR(),real),(SI(),imag),(LR(),real),(LI(),imag)]) do (which,comp)
        schurQR,history = partialschur(A; which=which, nev = 1, kwargs...)
        verbosity > 0 && println(history)
        comp(first(schurQR.eigenvalues))
    end

    for (label,(i,j)) in [("Real", (1,3)),("Imaginary", (2,4))]
        if abs(r[i]-r[j]) < ctol
            verbosity > 1 && @info "$label extent of spectral range $(abs(r[i]-r[j])) below tolerance $(ctol), conflating."
            r[i] = r[j] = (r[i]+r[j])/2
        end
    end

    for i = 1:4
        if abs(r[i]) < ctol
            verbosity > 1 && @info "$(abs(r[i])) below tolerance $(ctol), truncating."
            r[i] = zero(r[i])
        end
    end
    a,b = (r[1]+im*r[2], r[3]+im*r[4])

    if real(a) == real(b) || imag(a) == imag(b)
        Line(a,b)
    else
        # There could be cases where all the eigenvalues fall on a
        # sloped line in the complex plane, but we don't know how to
        # deduce that yet. The user is free to define such sloped
        # lines manually, though.
        Rectangle(a,b)
    end
end

function spectral_range(t, A; kwargs...)
    λ = spectral_range(A; kwargs...)
    ta,tb = extrema(t)
    ta*λ ∪ tb*λ
end
