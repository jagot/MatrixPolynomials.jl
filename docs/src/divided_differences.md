# Divided differences

The divided differences of a function ``f`` with respect to a set of
interpolation points ``\{\zeta_i\}`` is defined as
([McCurdy 1984](#Bibliography-1))
```math
\begin{equation}
\label{eqn:div-diff-def}
\divdiff(\zeta_{i:j})f \defd
\frac{1}{2\pi\im}
\oint
\diff{z}
\frac{f(z)}{(z-\zeta_i)(z-\zeta_{i+1})...(z-\zeta_j)},
\end{equation}
```
where the integral is taken along a simple contour encircling the
poles once. A common approach to evaluate the divided differences of
``f``, and an alternative definition, is the recursive scheme
```math
\begin{equation}
\label{eqn:div-diff-recursive}
\tag{\ref{eqn:div-diff-def}*}
\divdiff(\zeta_{i:j},z)f \defd
\frac{\divdiff(\zeta_{i:j-1},z)f-\divdiff(\zeta_{i:j})f}{z - \zeta_j}, \quad
\divdiff(\zeta_i,z)f \defd
\frac{\divdiff(z)f-\divdiff(\zeta_i)f}{z - \zeta_i}, \quad
\divdiff(z)f \defd f(z),
\end{equation}
```
which, however, is prone to catastrophic cancellation for very small
``\abs{\zeta_i-\zeta_j}``. This can be partially alleviated by
employing `BigFloat`s, but that will only postpone the breakdown,
albeit with ~40 orders of magnitude, which might be enough for
practical purposes (but much slower).

For specific choices of ``f``, efficient and accurate algorithms can
be devised. [`MatrixPolynomials.φₖ_ts_div_diff`](@ref) is based upon
the fact the divided differences in a third way can be computed as
([McCurdy 1984](#Bibliography-1), [Opitz 1964](#Bibliography-1))
```math
\begin{equation}
\label{eqn:div-diff-mat-fun}
\tag{\ref{eqn:div-diff-def}†}
\divdiff(\zeta_{i:j})f \defd
\vec{e}_1^\top
f(\mat{Z}_{i:j}),
\end{equation}
```
i.e. the first row of the function ``f`` applied to the matrix
```math
\begin{equation}
\mat{Z}_{i:j}\defd
\bmat{
\zeta_i&1&\\
&\zeta_{i+1}&1\\
&&\ddots&\ddots\\
&&&\ddots&1\\
&&&&\zeta_j}.
\end{equation}
```
The right-eigenvectors are given by ([Opitz 1964](#Bibliography-1))
```math
\begin{equation}
\label{eqn:div-diff-mat-right-eigen}
\mat{Q}_\zeta = \{q_{ik}\}, \quad
q_{ik} =
\begin{cases}
\prod_{j=i}^{k-1} (\zeta_k - \zeta_j)^{-1}, & i < k,\\
1, & i = k,\\
0, & \textrm{else},
\end{cases}
\end{equation}
```
and similarly, the left-eigenvectors are given by
```math
\begin{equation}
\label{eqn:div-diff-mat-left-eigen}
\tag{\ref{eqn:div-diff-mat-right-eigen}*}
\mat{Q}_\zeta^{-1} = \{\conj{q}_{ik}\}, \quad
\conj{q}_{ik} =
\begin{cases}
\prod_{j=i+1}^k (\zeta_i - \zeta_j)^{-1}, & i < k,\\
1, & i = k,\\
0, & \textrm{else},
\end{cases}
\end{equation}
```
such that
```math
\begin{equation}
\divdiff(\zeta_{i:j})f=
\mat{Q}_\zeta\mat{F}_\zeta\mat{Q}_\zeta^{-1},\quad
\mat{F}_\zeta \defd \bmat{f(\zeta_i)\\&f(\zeta_{i+1})\\&&\ddots\\&&&f(\zeta_j)}.
\end{equation}
```
However, straight evaluation of
``(\ref{eqn:div-diff-mat-right-eigen},\ref{eqn:div-diff-mat-left-eigen})``
is prone to the same kind of catastrophic cancellation as is
``\eqref{eqn:div-diff-recursive}``, so to evaluate
``\eqref{eqn:div-diff-mat-fun}``, one instead turns to Taylor or
[Padé](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant) expansions
of ``f(\mat{Z}_{i:j})`` ([McCurdy 1984](#Bibliography-1), [Caliari
2004](#Bibliography-1)), or interpolation polynomial basis changes
([Zivcovich 2019](#Bibliography-1)).

As an illustration, we show the divided differences of `exp` over 100
points uniformly spread over ``[-2,2]``, calculated using
``\eqref{eqn:div-diff-recursive}``, in `Float64` and `BigFloat`
precision, along with a Taylor expansion of
``\eqref{eqn:div-diff-mat-fun}``:

![Illustration of divided differences accuracy](figures/div_differences_cancellation.svg)

It can clearly be seen that the Taylor expansion is not susceptible to
the catastrophic cancellation; it is however, not valid outside the
interval ``[-1.59,1.59]``, so the domain of interest has to be
rescaled prior to its usage.

```@docs
MatrixPolynomials.⏃
MatrixPolynomials.std_div_diff
MatrixPolynomials.φₖ_ts_div_diff
MatrixPolynomials.φₖ_div_diff
MatrixPolynomials.div_diff_table
MatrixPolynomials.min_degree
MatrixPolynomials.taylor_series
```

## Bibliography

- Caliari, M. (2007). Accurate evaluation of divided differences for
  polynomial interpolation of exponential propagators. Computing,
  80(2), 189–201. [DOI:
  10.1007/s00607-007-0227-1](http://dx.doi.org/10.1007/s00607-007-0227-1)


- McCurdy, A. C., Ng, K. C., & Parlett, B. N. (1984). Accurate
  computation of divided differences of the exponential
  function. Mathematics of Computation, 43(168), 501–501. [DOI:
  10.1090/s0025-5718-1984-0758198-0](http://dx.doi.org/10.1090/s0025-5718-1984-0758198-0)

- Opitz, G. (1964). Steigungsmatrizen. ZAMM - Journal of Applied
  Mathematics and Mechanics / Zeitschrift für Angewandte Mathematik
  und Mechanik, 44(S1), [DOI:
  10.1002/zamm.19640441321](http://dx.doi.org/10.1002/zamm.19640441321)

- Zivcovich, F. (2019). Fast and accurate computation of divided
  differences for analytic functions, with an application to the
  exponential function. Dolomites Research Notes on Approximation,
  12(1), 28–42. [PDF:
  Zivcovich_2019_FAC.pdf](https://drna.padovauniversitypress.it/system/files/papers/Zivcovich_2019_FAC.pdf)
