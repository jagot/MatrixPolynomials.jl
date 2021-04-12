# Leja points

A common problem in polynomial interpolation of functions, is that
when the number of interpolation points is increased, the
interpolation polynomial becomes ill-conditioned
([overfitting](https://en.wikipedia.org/wiki/Overfitting)). It can be
shown that interpolation at the roots of the [Chebyshev polynomials
](https://en.wikipedia.org/wiki/Chebyshev_polynomials) yields the best
approximation, however, it is difficult to generate successively
better approximations, since the roots of the Chebyshev polynomial of
degree ``m`` are not related to those of the polynomial of degree
``m-1``.

The Leja points [^Leja] ``\{\zeta_i\}`` are generated from a set ``E
\subset \Complex`` such that the next point in the sequence is
maximally distant from all previously generated points:

```math
\begin{equation}
w(\zeta_j)
\prod_{k=0}^{j-1} \abs{\zeta_j-\zeta_k} =
\max_{\zeta\in E}
w(\zeta)
\prod_{k=0}^{j-1}
\abs{\zeta - \zeta_k},
\end{equation}
```
with ``w(\zeta)`` being an optional weight function (unity
hereinafter). Interpolating a function on the Leja points largely
avoids the overfitting problems and performs similarly to Chebyshev
interpolation [^Reichel], while still allowing
for iteratively improved approximation by the addition of more
interpolation points.

MatrixPolynomials.jl provides two methods for generating the Leja
points, [`MatrixPolynomials.Leja`](@ref) and
[`MatrixPolynomials.FastLeja`](@ref). The figure below illustrates the
distribution of Leja points using both methods, on the line
``[-2,2]``, for the [`MatrixPolynomials.Leja`](@ref), an underlying
discretization of 1000 points was employed, and 10 Leja points were
generated. The lower part of the plot shows the estimation of the
[capacity](https://en.wikipedia.org/wiki/Capacity_of_a_set),
calculated as
```math
C(\{\zeta_{1:m}\}) \approx
\left|\left(\prod_{i=1}^{m-1} |\zeta_m-\zeta_i|\right)\right|^{1/m}.
```
For the set ``[-2,2]``, the capacity is unity, which is approached for
increasing values of ``m``.

```julia-repl
julia> import MatrixPolynomials: Leja, FastLeja

julia> m = 10
10

julia> a,b = -2,2
(-2, 2)

julia> l = Leja(range(a, stop=b, length=1000), m)
Leja{Float64}([-1.995995995995996, -1.991991991991992, -1.987987987987988, -1.983983983983984, -1.97997997997998, -1.975975975975976, -1.971971971971972, -1.967967967967968, -1.9639639639639639, -1.95995995995996  …  1.95995995995996, 1.9639639639639639, 1.967967967967968, 1.971971971971972, 1.975975975975976, 1.97997997997998, 1.983983983983984, 1.987987987987988, 1.991991991991992, 1.995995995995996], [2.0, -2.0, -0.002002002002002002, 1.155155155155155, -1.3193193193193193, 1.6796796796796796, -1.7397397397397398, -0.6106106106106106, 0.6426426426426426, 1.887887887887888], [0.0, 4.0, 3.9999959919879844, 3.084537289340691, 7.36488275292736, 3.118030920568761, 7.038861956228758, 7.143962613999413, 7.199339458696, 4.549146401863414])

julia> fl = FastLeja(a, b, m)
FastLeja{Float64}([2.0, -2.0, 0.0, -1.0, 1.0, -1.5, 1.5, 0.5, -1.75, 1.75], [-3.111827946268022, -1.5140533447265625, 7.91015625, -1.3255691528320312, -3.0929946899414062, 0.6718902150169015, 1.1896133422851562, 1.2691259616985917, -1.8015846004709601, 2.6076411906e-314], [-1.875, 0.25, -0.5, 1.25, -1.25, 1.625, 0.75, -1.625, 1.875, 1.5e-323], [2, 3, 4, 5, 6, 7, 8, 9, 10, 2], [9, 8, 3, 7, 4, 10, 5, 6, 1, 4570435120])
```

![Leja points](figures/leja_points.svg)

## Reference

```@docs
MatrixPolynomials.Leja
MatrixPolynomials.Leja(S::AbstractVector{T}, n::Integer) where T
MatrixPolynomials.leja!
MatrixPolynomials.FastLeja
MatrixPolynomials.fast_leja!
MatrixPolynomials.points
```

## Bibliography

[^Leja]: Leja, F. (1957). Sur certaines suites liées aux ensembles
    plans et leur application à la représentation conforme. Annales
    Polonici Mathematici, 4(1), 8–13. [DOI:
    10.4064/ap-4-1-8-13](http://dx.doi.org/10.4064/ap-4-1-8-13)

[^Reichel]: Reichel, L. (1990). Newton Interpolation At Leja
    Points. BIT, 30(2), 332–346. [DOI:
    10.1007/bf02017352](http://dx.doi.org/10.1007/bf02017352)

[^Baglama]: Baglama, J., Calvetti, D., & Reichel, L. (1998). Fast Leja
    points. Electron. Trans. Numer. Anal, 7(124-140), 119–120. [URL:
    https://elibm.org/article/10006464](https://elibm.org/article/10006464)
