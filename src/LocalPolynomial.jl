# LocalPolynomial — the shared local-polynomial jet-smoother: the kernel-weighted local-polynomial
# design solve over a binned 1-D measure, returning the local Taylor coefficients (the jet). This is
# the single *linear* primitive on which the local smoothers are built — `ConditionalExpectations`'s local-linear
# fit is one `local_fit` call, and `MarginalDensity`'s local-likelihood density is `local_fit`'s
# iteratively reweighted (IRLS) inner solve: a local GLM is this linear solve iterated
# (`docs/ARCHITECTURE.md` §2, §11). Pure host-side numerics over recovered arrays. Kernel-agnostic via the
# shared `Kernels._kernel_weight` atom — any symmetric callable at any degree (its normalization cancels);
# an asymmetric boundary kernel has no jet, so it is degree-0 only.

module LocalPolynomial

using LinearAlgebra

using ..Kernels: _kernel_weight, AsymmetricKernel


# `local_fit` (below) is the package-internal shared primitive — not exported and not a manual entry
# (like `ConditionalExpectations._column_sums`); consumers reach it by explicit `using ..LocalPolynomial: local_fit`.


# Kernel-weighted design moments of a binned measure about the query point `x`:
#   moments[p+1] = Σ_j _kernel_weight(c̃_j, x, kernel, h) · weights[j] · (c̃_j − x)^p,   p = 0 … maxorder,
# with c̃_j = `locations[j]`.
# `weights` is the per-bin weight — the bin counts for the design Gram, a response total for the
# right-hand side, or the IRLS working weight (n̂_j, or the residual count n_j − n̂_j) for local likelihood.
# The weight comes from the shared `_kernel_weight` atom (the 1/h it folds in cancels in the solve);
# powers are built iteratively in one pass over the bins.
function _design_moments(locations, weights, x, kernel, maxorder, h)
  moments = zeros(Float64, maxorder + 1)
  @inbounds for j in eachindex(locations)
    d = locations[j] - x
    w = _kernel_weight(locations[j], x, kernel, h) * weights[j]
    p = w
    moments[1] += p
    for k in 1:maxorder
      p *= d
      moments[k + 1] += p
    end
  end
  return moments
end


# Solve the local-polynomial normal equations: the (d+1)×(d+1) symmetric design matrix
# U = [U_{r+s}]_{r,s=0}^d (the Hankel matrix of U_0 … U_{2d}; in code M[i, j] = moments[i+j-1]) against
# `rhs` (= V_0 … V_d), returning the jet β̂ = (β̂_0, …, β̂_d). The design matrix is positive
# semidefinite; a degenerate local design (singular — too little local data near x) yields an all-NaN
# jet, the caller's signal that x is unsupported.
function _hankel_solve(moments, rhs)
  d1 = length(rhs)
  M = Matrix{Float64}(undef, d1, d1)
  @inbounds for j in 1:d1, i in 1:d1
    M[i, j] = moments[i + j - 1]
  end
  beta = try
    M \ rhs
  catch
    return fill(NaN, d1)
  end
  return all(isfinite, beta) ? beta : fill(NaN, d1)
end


"""
    local_fit(locations, mass, response, x, kernel, degree, h)

The *local-polynomial jet* of a binned 1-D measure at the query point `x` is the kernel-weighted
least-squares fit of a degree-``d`` polynomial in the bin-location offsets ``\\tilde c_j - x`` to the bin data, returning the
coefficient vector (jet) ``\\hat\\beta = (\\hat\\beta_0, \\dots, \\hat\\beta_d)`` (with ``d`` the `degree` argument).
From bin locations ``\\tilde c_j``, design weights ``n_j`` (`mass`) and response totals ``n_j \\bar y_j``
(`response`, count times mean response), and `kernel` scaled by the bandwidth `h`,
``K_h(v) \\propto K(v/h)``, form the design and response
moments

```math
\\begin{aligned}
U_p &= \\sum_j K_h(\\tilde c_j - x)  n_j  (\\tilde c_j - x)^p, \\\\
V_p &= \\sum_j K_h(\\tilde c_j - x)  n_j \\bar y_j  (\\tilde c_j - x)^p,
\\end{aligned}
```

and solve the normal equations ``[ U_{r+s} ]_{r,s=0}^{d} \\hat\\beta = (V_0, \\dots, V_d)``. The
jet is the local Taylor expansion: ``\\hat\\beta_0`` is the fitted value, ``\\hat\\beta_1`` the slope, and
``2\\hat\\beta_2`` the curvature — so a regression fit gives the regression-function estimate
``\\hat m(x) = \\hat\\beta_0`` and ``\\hat m''(x) = 2\\hat\\beta_2``. The kernel weight's normalization cancels in
the solve, so any symmetric callable `kernel` works at any `degree`; an asymmetric boundary
kernel has no jet and is admitted only at `degree = 0` (the Nadaraya–Watson ratio). A degenerate local design (too little data near `x`) returns an all-NaN jet.

This is the shared linear primitive of the local smoothers:
``\\hat m(x)`` is ``\\hat\\beta_0`` of one `local_fit`, and the local-likelihood density reuses the same
kernel-weighted-moment solve as the iteratively reweighted (IRLS) inner step of its local-Poisson fit.

# Arguments

- `locations` — the bin locations ``\\tilde c_j`` (bin centers, or centroids from the stored sum of
  within-bin offsets).
- `mass` — the per-bin design weight ``n_j`` (the bin counts).
- `response` — the per-bin response total ``n_j \\bar y_j``, with ``\\bar y_j`` the bin's mean response; passing
  `mass` as `response` smooths the counts (as the local-likelihood density does) to a
  count level.
- `x` — the query point at which the jet is evaluated.
- `kernel` — the smoothing kernel: a symmetric callable ``K`` (admitted at any `degree`), or an
  [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel) (`degree = 0` only).
- `degree` — the local-polynomial order ``d``: `0` (local-constant / Nadaraya–Watson), `1`
  (local-linear), or `2` (local-quadratic); the returned jet then has `degree + 1` entries.
- `h` — the smoothing bandwidth (``h > 0``); it scales the kernel argument.
"""
function local_fit(locations, mass, response, x, kernel, degree, h)
  (kernel isa AsymmetricKernel && degree ≥ 1) && throw(ArgumentError(
    "an asymmetric kernel is a degree-0 boundary kernel — it has no jet; `degree` must be 0 (got $degree)"))
  T = _design_moments(locations, mass, x, kernel, 2 * degree, float(h))
  P = _design_moments(locations, response, x, kernel, degree, float(h))
  return _hankel_solve(T, P)
end


# The leverage ℓ(x) of the local-polynomial value β̂_0 at x: Σ_j N_j ω_j², the sum of squared effective
# regression weights (ω_j the weight `local_fit`'s β̂_0 places on column j's response, N_j = `weights[j]`).
# This is the variance factor of β̂_0 — for a response whose per-column conditional variance is σ²(x),
# Var(β̂_0) ≈ σ²(x)·ℓ(x) — so a standard error is √(ℓ(x) · noise-scale). Writing
# ω_j = e_0ᵀ U⁻¹ (K_h(c_j−x)(c_j−x)^p)_{p=0..d} with U = [U_{r+s}]_{r,s=0}^d the design Hankel,
# U_p = Σ_j N_j K_h(c_j−x)(c_j−x)^p, gives ℓ(x) = Σ_j N_j ω_j² = e_0ᵀ U⁻¹ W U⁻¹ e_0, where W = [W_{r+s}]
# is the Hankel of the squared-kernel moments W_p = Σ_j N_j K_h(c_j−x)² (c_j−x)^p
# (the 1/h in K_h cancels between U and W). Degenerate design (too little local data) ⇒ NaN.
function _leverage(locations, weights, x, kernel, degree, h)
  hf = float(h)
  U_p = zeros(Float64, 2 * degree + 1)         # plain-kernel design moments U_p = Σ_j N_j K_h(c_j−x)(c_j−x)^p
  W_p = zeros(Float64, 2 * degree + 1)         # squared-kernel moments      W_p = Σ_j N_j K_h(c_j−x)²(c_j−x)^p
  @inbounds for j in eachindex(locations)
    wj = weights[j]
    iszero(wj) && continue
    kw = _kernel_weight(locations[j], x, kernel, hf)
    iszero(kw) && continue
    dloc = locations[j] - x
    p = 1.0
    for l in 0:2 * degree
      U_p[l + 1] += wj * kw * p
      W_p[l + 1] += wj * kw^2 * p
      p *= dloc
    end
  end
  d1 = degree + 1
  U = Matrix{Float64}(undef, d1, d1)
  W = Matrix{Float64}(undef, d1, d1)
  @inbounds for j in 1:d1, i in 1:d1
    U[i, j] = U_p[i + j - 1]
    W[i, j] = W_p[i + j - 1]
  end
  e_0 = zeros(Float64, d1)
  e_0[1] = 1.0
  u = try
    U \ e_0
  catch
    return NaN
  end
  all(isfinite, u) || return NaN
  return dot(u, W, u)                      # e_0ᵀ U⁻¹ W U⁻¹ e_0 = Σ_j N_j ω_j²
end


end # module LocalPolynomial
