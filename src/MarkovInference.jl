# MarkovInference — §4.4: conditional increment moments and their standard functionals
# (Kramers–Moyal coefficients, drift, diffusion, the Pawula/Q diagnostic, the Markov test). Defined
# and estimated WITHOUT assuming the data is a diffusion — model-neutrality is the invariant here;
# the SDE reading stays consumer-side (`docs/ARCHITECTURE.md` §4.4).

module MarkovInference

using ..Histograms: Histogram2D
using ..Regression: _column_sums, _local_linear

export increment_moment, kramers_moyal, drift, diffusion, pawula_quotient, markov_test


# Conditioning on X_t = x makes the increment ΔX = Y − x with Y = X_{t+τ} (the response stored in the
# joint summary). From the order-2 power-sums only M_1 and M_2 are recoverable:
#   M_1(x) = E[Y|x] − x,   M_2(x) = E[(Y−x)²|x] = E[Y²|x] − 2x E[Y|x] + x².
# E[Y|x] and E[Y²|x] are the local-linear fits already used by `Regression` (E[Y²|x] is the fit of the
# Y² response). One per-column reduction feeds both. Returns (M_1, M_2) at a single x.
@inline function _increment_moments(cx, N, sumY, sumY2, x, invh)
  m = _local_linear(cx, N, sumY, x, invh)       # Ê[Y | X = x]
  ey2 = _local_linear(cx, N, sumY2, x, invh)    # Ê[Y² | X = x]
  return (m - x), (ey2 - 2 * x * m + x^2)
end

@inline _select_moment(M1, M2, n) = n == 1 ? M1 : M2


# n ≥ 3 needs per-cell Y power-sums of order n, beyond the order-2 Histogram2D — a recorded coverage
# gap (invariant 5), thrown rather than silently approximated.
function _check_order(n::Integer)
  1 ≤ n ≤ 2 || throw(ArgumentError(
    "increment-moment order n = $n is outside {1, 2}: the order-2 Histogram2D resolves only M_1 " *
    "(drift) and M_2 (diffusion). Higher n (the higher KM coefficients, and the Pawula Q's M_4) needs the joint summary raised to " *
    "per-cell Y power-sums of order n — a recorded coverage gap (docs/ARCHITECTURE.md §4.4)."))
  return nothing
end

_check_bandwidth(h::Real) =
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))


"""
    increment_moment(summary::Histogram2D, x, n; h)   # ≡ M_n(x, τ)

**Estimand.** the conditional increment moment
``M_n(x, \\tau) = \\mathbb{E}\\!\\left[(X_{t+\\tau} - X_t)^n \\mid X_t = x\\right]`` — the ``n``-th moment
of the lag-``\\tau`` increment ``\\Delta X = X_{t+\\tau} - X_t`` conditioned on the present value ``x``.
A pure conditional expectation; no model is assumed.

**Estimator.** writing ``Y = X_{t+\\tau}`` for the response stored in the joint
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) at lag ``\\tau`` and conditioning on
``X_t = x`` (so ``\\Delta X = Y - x``), the first two moments are conditional moments of ``Y``:

```math
M_1(x, \\tau) = \\hat m_h(x) - x, \\qquad
M_2(x, \\tau) = \\widehat{\\mathbb{E}}[Y^2 \\mid x] - 2x\\,\\hat m_h(x) + x^2 ,
```

with ``\\hat m_h(x) = \\widehat{\\mathbb{E}}[Y \\mid x]`` the local-linear
[`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean) and
``\\widehat{\\mathbb{E}}[Y^2 \\mid x]`` the same fit applied to the ``Y^2`` response — both exact from
the per-cell power-sums. Hence ``n \\in \\{1, 2\\}`` (drift and diffusion); ``n \\ge 3`` (the higher
Kramers–Moyal coefficients, and the Pawula ``Q``, which needs ``M_4``) needs the summary raised to
order-``n`` per-cell ``Y`` power-sums — a recorded
coverage gap, so it throws rather than silently approximate. `h` is the supplied smoothing bandwidth;
a vector `x` returns a vector.
"""
function increment_moment(summary::Histogram2D, x::Real, n::Integer; h::Real)
  _check_order(n)
  _check_bandwidth(h)
  cx, N, sumY, sumY2 = _column_sums(summary)
  M1, M2 = _increment_moments(cx, N, sumY, sumY2, x, inv(float(h)))
  return _select_moment(M1, M2, n)
end

function increment_moment(summary::Histogram2D, xs::AbstractVector, n::Integer; h::Real)
  _check_order(n)
  _check_bandwidth(h)
  cx, N, sumY, sumY2 = _column_sums(summary)
  invh = inv(float(h))
  return [_select_moment(_increment_moments(cx, N, sumY, sumY2, x, invh)..., n) for x in xs]
end


"""
    kramers_moyal(summary::Histogram2D, x, n; h, τ = summary.lag)   # ≡ D⁽ⁿ⁾(x)

**Estimand.** the ``n``-th Kramers–Moyal coefficient
``D^{(n)}(x) = \\dfrac{1}{n!}\\,\\lim_{\\tau \\to 0} \\dfrac{M_n(x, \\tau)}{\\tau}``, with ``M_n`` the
conditional [`increment_moment`](@ref GBGStreamingStatistics.MarkovInference.increment_moment).

**Estimator.** the finite-``\\tau`` ratio
``\\widehat{D}^{(n)}(x) = \\widehat{M}_n(x, \\tau)\\,/\\,(n!\\,\\tau)`` at the summary's lag ``\\tau``.
The keyword `τ` is the **physical lag time**; it defaults to the integer `lag` label of the summary
(exact when the sampling interval is the time unit). For small ``\\tau`` this approximates ``D^{(n)}``
to ``O(\\tau)``; the all-orders finite-time (Bell-polynomial) correction that inverts the
``\\tau``-expansion rather than extrapolating (Rydin Gorjão et al. 2021) is forthcoming. Supports
``n \\in \\{1, 2\\}`` — see [`drift`](@ref GBGStreamingStatistics.MarkovInference.drift) and
[`diffusion`](@ref GBGStreamingStatistics.MarkovInference.diffusion).

> **Model-neutral.** ``D^{(n)}`` is defined and estimated without assuming the data *is* a diffusion.
> That ``D^{(1)}`` "is" drift, that the ``D^{(n)}`` are coefficients of a Fokker–Planck operator —
> these are the *consumer's* interpretation; the package supplies the estimator and the tests of the
> hypothesis (see [`pawula_quotient`](@ref GBGStreamingStatistics.MarkovInference.pawula_quotient),
> [`markov_test`](@ref GBGStreamingStatistics.MarkovInference.markov_test)).

Reference: A. Rydin Gorjão et al. (2021), [*Arbitrary-Order Finite-Time Corrections for the Kramers–Moyal Operator*](https://doi.org/10.3390/e23050517), Entropy 23(5).
"""
function kramers_moyal(summary::Histogram2D, x, n::Integer; h::Real, τ::Real = summary.lag)
  τ > 0 || throw(ArgumentError(
    "kramers_moyal needs a positive lag time τ (got τ = $τ; the summary has lag = $(summary.lag)). " *
    "Pass τ = (physical lag time), or accumulate the Histogram2D at a lag > 0."))
  return increment_moment(summary, x, n; h = h) ./ (factorial(n) * τ)
end


"""
    drift(summary::Histogram2D, x; h, τ = summary.lag)   # ≡ D⁽¹⁾(x)

**Estimand.** the first Kramers–Moyal coefficient
``D^{(1)}(x) = \\lim_{\\tau \\to 0} M_1(x, \\tau)/\\tau`` — the ``n = 1``
[`kramers_moyal`](@ref GBGStreamingStatistics.MarkovInference.kramers_moyal) coefficient. For an Itô SDE
``\\mathrm{d}X = D^{(1)}(X)\\,\\mathrm{d}t + \\sqrt{2 D^{(2)}(X)}\\,\\mathrm{d}W`` a consumer reads it
as the drift; the package asserts no such model.

**Estimator.** ``\\widehat{D}^{(1)}(x) = \\bigl(\\hat m_h(x) - x\\bigr)/\\tau`` — the local-linear
increment mean over the lag time ``\\tau``. `h` is the supplied bandwidth; a vector `x` returns a
vector.
"""
drift(summary::Histogram2D, x; h::Real, τ::Real = summary.lag) =
  kramers_moyal(summary, x, 1; h = h, τ = τ)


"""
    diffusion(summary::Histogram2D, x; h, τ = summary.lag)   # ≡ D⁽²⁾(x)

**Estimand.** the second Kramers–Moyal coefficient
``D^{(2)}(x) = \\tfrac{1}{2}\\lim_{\\tau \\to 0} M_2(x, \\tau)/\\tau`` — the ``n = 2``
[`kramers_moyal`](@ref GBGStreamingStatistics.MarkovInference.kramers_moyal) coefficient (the
``\\tfrac{1}{2} = 1/2!`` is the Kramers–Moyal normalization). For the SDE
``\\mathrm{d}X = D^{(1)}\\,\\mathrm{d}t + \\sqrt{2 D^{(2)}}\\,\\mathrm{d}W`` a consumer reads
``\\sqrt{2 D^{(2)}}`` as the noise amplitude (so ``D^{(2)} = \\sigma^2/2`` for additive
``\\sigma\\,\\mathrm{d}W``); the package asserts no such model.

**Estimator.** ``\\widehat{D}^{(2)}(x) = M_2(x, \\tau)/(2\\tau)`` with ``M_2`` the conditional
[`increment_moment`](@ref GBGStreamingStatistics.MarkovInference.increment_moment). `h` is the supplied
bandwidth; a vector `x` returns a vector.
"""
diffusion(summary::Histogram2D, x; h::Real, τ::Real = summary.lag) =
  kramers_moyal(summary, x, 2; h = h, τ = τ)


"""
    pawula_quotient(summary, x; ...)   # ≡ Q(x)

**Estimand.** the Pawula / ``Q`` diagnostic — the standardized fourth conditional-increment moment
``Q(x) = \\mu_4(x,\\tau) / \\bigl(3\\,\\mu_2(x,\\tau)^2\\bigr)``, with ``\\mu_n`` the *central* conditional
increment moments (about the conditional mean increment ``M_1``, so the drift does not contaminate the
ratio; recovered from the [`increment_moment`](@ref GBGStreamingStatistics.MarkovInference.increment_moment)
``M_n`` by the shift formula). A diffusion has Gaussian conditional increments, for which
``\\mu_4 = 3\\mu_2^2``, so ``Q(x) \\to 1`` — a *test* of whether the higher Kramers–Moyal coefficients
vanish (Pawula's theorem), not a premise: ``Q \\approx 1`` supports a diffusion reading, while the excess
``Q - 1`` (the increment kurtosis) flags non-Gaussian, jump-like dynamics.

!!! note "Forthcoming — needs the order raise"
    Declared; not yet implemented. ``\\mu_4`` requires the joint summary raised to per-cell ``Y``
    power-sums of order 4, beyond the current order-2
    [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) — a recorded coverage gap.
"""
function pawula_quotient end


"""
    markov_test(summaries...; ...)

A ``\\chi^2`` **Markov / Chapman–Kolmogorov test**: a hypothesis test of whether the data is
Markovian, comparing the two-step transition density against the Chapman–Kolmogorov composition of
one-step densities (see
[`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density)).
Intrinsically finite-sample — a test, not an estimand.

!!! note "Forthcoming"
    Declared; not yet implemented — it composes the now-available §4.3
    [`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density).
"""
function markov_test end

end # module MarkovInference
