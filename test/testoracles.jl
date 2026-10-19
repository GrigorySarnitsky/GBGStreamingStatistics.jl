# Shared analytic oracles for the test suite: closed-form densities, the Ornstein–Uhlenbeck transition
# laws, and small numeric helpers (trapezoid integral, ISE). Pure functions — no `Test` or package
# coupling — so every test file can `include("testoracles.jl"); using .TestOracles` and stay
# standalone-runnable. The single source of truth for the analytic targets the estimators are checked
# against (so a density's pdf or the OU drift is written once, not re-derived per file).
module TestOracles

using SpecialFunctions: loggamma, logbeta

export gaussian_pdf, exponential_pdf, lognormal_pdf, gamma_pdf, beta_pdf, student_t_pdf, uniform_pdf,
       inverse_gaussian_pdf, claw_pdf,
       ou_stationary_std, ou_drift, ou_diffusion, ou_transition,
       trapz, ise


# --- analytic densities (zero outside the support) ---

"Gaussian density ``\\mathcal{N}(\\mu, \\sigma^2)``."
gaussian_pdf(x; mu = 0.0, sigma = 1.0) = exp(-((x - mu) / sigma)^2 / 2) / (sigma * sqrt(2π))

"Exponential density, rate ``\\lambda`` (support ``[0, \\infty)``)."
exponential_pdf(x; rate = 1.0) = x ≥ 0 ? rate * exp(-rate * x) : 0.0

"Log-normal density with log-mean ``\\mu`` and log-sd ``\\sigma`` (support ``(0, \\infty)``)."
lognormal_pdf(x; mu = 0.0, sigma = 1.0) =
  x > 0 ? exp(-((log(x) - mu)^2) / (2 * sigma^2)) / (x * sigma * sqrt(2π)) : 0.0

"Gamma density, shape ``k`` and scale ``\\theta`` (support ``(0, \\infty)``)."
gamma_pdf(x; shape, scale = 1.0) =
  x > 0 ? exp((shape - 1) * log(x) - x / scale - shape * log(scale) - loggamma(shape)) : 0.0

"Beta density with shape parameters ``a, b`` (support ``(0, 1)``)."
beta_pdf(x; a, b) =
  (0 < x < 1) ? exp((a - 1) * log(x) + (b - 1) * log1p(-x) - logbeta(a, b)) : 0.0

"Inverse-Gaussian density, mean ``\\mu`` and shape ``\\lambda`` (support ``(0,\\infty)``)."
inverse_gaussian_pdf(x; mu, lambda) =
  x > 0 ? sqrt(lambda / (2π * x^3)) * exp(-lambda * (x - mu)^2 / (2 * mu^2 * x)) : 0.0

"Student-``t`` density with ``\\nu`` degrees of freedom."
student_t_pdf(x; nu) =
  exp(loggamma((nu + 1) / 2) - loggamma(nu / 2)) / sqrt(nu * π) * (1 + x^2 / nu)^(-(nu + 1) / 2)

"Uniform density on ``[\\texttt{lower}, \\texttt{upper})``."
uniform_pdf(x; lower = 0.0, upper = 1.0) = (lower ≤ x < upper) ? 1.0 / (upper - lower) : 0.0

"Marron–Wand (1992) **claw** density: ``\\tfrac12 \\mathcal{N}(0,1) + \\sum_{l=0}^{4} \\tfrac1{10} \\mathcal{N}(l/2-1,\\,(1/10)^2)`` — a broad Gaussian with five sharp spikes, the standard hard-multimodal benchmark for bandwidth selection."
claw_pdf(x) = 0.5 * gaussian_pdf(x) + sum(0.1 * gaussian_pdf(x; mu = l / 2 - 1, sigma = 0.1) for l in 0:4)


# --- Ornstein–Uhlenbeck  dX = −κ X dt + σ dW : stationary 𝒩(0, σ_X²) with σ_X² = σ²/(2κ); the lag-τ
# transition is X_{t+τ} | X_t = x ∼ 𝒩(x e^{−κτ}, σ_X²(1 − e^{−2κτ})), so the finite-τ Kramers–Moyal
# coefficients are closed-form. The canonical drift/diffusion oracle (test/oracles.jl). ---

"Stationary standard deviation ``\\sigma_X = \\sqrt{\\sigma^2/(2\\kappa)}`` of the OU process."
ou_stationary_std(; kappa, sigma) = sqrt(sigma^2 / (2 * kappa))

"OU finite-``\\tau`` drift ``D^{(1)}(x,\\tau) = x(e^{-\\kappa\\tau}-1)/\\tau``."
ou_drift(x; kappa, tau) = x * (exp(-kappa * tau) - 1) / tau

"OU finite-``\\tau`` diffusion ``D^{(2)}(x,\\tau) = [\\sigma_X^2(1-e^{-2\\kappa\\tau}) + x^2(e^{-\\kappa\\tau}-1)^2]/(2\\tau)``."
function ou_diffusion(x; kappa, sigma, tau)
  variance_X = sigma^2 / (2 * kappa)
  return (variance_X * (1 - exp(-2 * kappa * tau)) + x^2 * (exp(-kappa * tau) - 1)^2) / (2 * tau)
end

"OU lag-``\\tau`` transition density ``p_\\tau(y\\mid x) = \\mathcal{N}(x e^{-\\kappa\\tau},\\, \\sigma_X^2(1-e^{-2\\kappa\\tau}))``."
function ou_transition(y, x; kappa, sigma, tau)
  v = (sigma^2 / (2 * kappa)) * (1 - exp(-2 * kappa * tau))
  return exp(-(y - x * exp(-kappa * tau))^2 / (2 * v)) / sqrt(2π * v)
end


# --- numeric helpers (uniform grid of spacing Δ) ---

"Trapezoidal integral of `values` sampled on a uniform grid of spacing `Δ`."
trapz(values, Δ) = (sum(values) - (first(values) + last(values)) / 2) * Δ

"Integrated squared error ``\\int (\\hat f - f)^2`` of `estimate` against `truth` on a grid of spacing `Δ`."
ise(estimate, truth, Δ) = sum(abs2, estimate .- truth) * Δ

end # module TestOracles
