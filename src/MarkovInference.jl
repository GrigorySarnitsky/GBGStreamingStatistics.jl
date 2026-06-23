# MarkovInference — §4.4: conditional increment moments and their standard functionals
# (Kramers–Moyal coefficients, drift, diffusion, the Pawula/Q diagnostic, the Markov test). Defined
# and estimated WITHOUT assuming the data is a diffusion — model-neutrality is the invariant here;
# the SDE reading stays consumer-side (`docs/ARCHITECTURE.md` §4.4).

module MarkovInference

using ..Histograms: Histogram2D

export increment_moment, kramers_moyal, drift, diffusion, pawula_quotient, markov_test


"""
    increment_moment(summary::Histogram2D, x, n; ...)   # ≡ M_n(x,τ)

**Estimand.** the conditional increment moment
``M_n(x, \\tau) = \\mathbb{E}\\!\\left[(X_{t+\\tau} - X_t)^n \\mid X_t = x\\right]`` — a pure
statistical object (a conditional expectation of a power of an increment).

**Estimator.** [`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean) applied to the response ``(\\Delta X)^n`` over the per-lag
joint histogram `summary`.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function increment_moment end


"""
    kramers_moyal(summary, x, n; ...)   # ≡ D⁽ⁿ⁾(x)

**Estimand.** the Kramers–Moyal coefficient
``D^{(n)}(x) = \\lim_{\\tau \\to 0} \\frac{1}{n!\\,\\tau}\\, M_n(x, \\tau)``, with ``M_n`` the
conditional [`increment_moment`](@ref GBGStreamingStatistics.MarkovInference.increment_moment).

**Estimator.** the finite-``\\tau`` increment moments with the **all-orders Bell-polynomial
correction** (Rydin Gorjão et al. 2021), which inverts the ``\\tau``-expansion rather than
extrapolating ``\\tau \\to 0``.

> **Model-neutral.** ``D^{(n)}`` is defined and estimated without assuming the data *is* a diffusion.
> That ``D^{(1)}`` "is" drift, that the ``D^{(n)}`` are coefficients of a Fokker–Planck operator —
> these are the *consumer's* interpretation; the package supplies the estimator and the tests of the
> hypothesis (see [`pawula_quotient`](@ref GBGStreamingStatistics.MarkovInference.pawula_quotient), [`markov_test`](@ref GBGStreamingStatistics.MarkovInference.markov_test)).

Reference: A. Rydin Gorjão et al. (2021), *Kramers–Moyal coefficients with finite-time corrections*.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function kramers_moyal end


"""
    drift(summary, x; ...)   # ≡ D⁽¹⁾(x)

**Estimand.** the first Kramers–Moyal coefficient ``D^{(1)}(x)`` — the ``n = 1`` case of
[`kramers_moyal`](@ref GBGStreamingStatistics.MarkovInference.kramers_moyal). A consumer reads it as the SDE drift; the package does not assert the model.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function drift end


"""
    diffusion(summary, x; ...)   # ≡ D⁽²⁾(x)

**Estimand.** the second Kramers–Moyal coefficient ``D^{(2)}(x)`` — the ``n = 2`` case of
[`kramers_moyal`](@ref GBGStreamingStatistics.MarkovInference.kramers_moyal). A consumer reads it as the SDE diffusion; the package does not assert the
model.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function diffusion end


"""
    pawula_quotient(summary, x; ...)   # ≡ Q(x)

**Estimand.** the Pawula / ``Q`` diagnostic
``Q(x) = M_6(x, \\tau) / \\bigl(5\\, M_4(x, \\tau)^2\\bigr)`` (with ``M_n`` the conditional
[`increment_moment`](@ref GBGStreamingStatistics.MarkovInference.increment_moment)) — a *test* of whether the higher Kramers–Moyal coefficients vanish
(Pawula's theorem), not a premise. A small ``Q`` supports a diffusion reading; a large one refutes
it.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function pawula_quotient end


"""
    markov_test(summaries...; ...)

A ``\\chi^2`` **Markov / Chapman–Kolmogorov test**: a hypothesis test of whether the data is
Markovian, comparing the two-step transition density against the Chapman–Kolmogorov composition of
one-step densities (see [`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density)). Intrinsically finite-sample — a test, not an
estimand.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function markov_test end

end # module MarkovInference
