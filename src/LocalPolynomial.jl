# LocalPolynomial — the shared local-polynomial jet-smoother: the kernel-weighted local-polynomial
# design solve over a binned 1-D measure, returning the local Taylor coefficients (the jet). This is
# the single *linear* primitive on which the local smoothers are built — `Regression`'s local-linear
# fit is one `local_fit` call, and `MarginalDensity`'s local-likelihood density is `local_fit`'s
# iteratively reweighted (IRLS) inner solve: a local GLM is this linear solve iterated
# (`docs/ARCHITECTURE.md` §2, §11). Pure host-side numerics over recovered arrays. Kernel-agnostic via the
# shared `Kernels._kernel_weight` atom — any symmetric callable at any degree (its normalization cancels);
# an asymmetric boundary kernel has no jet, so it is degree-0 only.

module LocalPolynomial

using LinearAlgebra

using ..Kernels: _kernel_weight, AsymmetricKernel


# `local_fit` (below) is the package-internal shared primitive — not exported and not a manual entry
# (like `Regression._column_sums`); consumers reach it by explicit `using ..LocalPolynomial: local_fit`.


# Kernel-weighted design moments of a binned measure about the query point `x`:
#   moments[k+1] = Σ_j _kernel_weight(kernel, loc_j, x, h) * weights_j * (loc_j - x)^k,   k = 0 … maxorder.
# `weights` is the per-bin weight — the bin counts for the design Gram, a response total for the
# right-hand side, or the IRLS working weight (μ, or the residual count n - μ) for local likelihood.
# The weight comes from the shared `_kernel_weight` atom (the 1/h it folds in cancels in the solve);
# powers are built iteratively in one pass over the bins.
function _design_moments(locations, weights, x, h, kernel, maxorder)
  moments = zeros(Float64, maxorder + 1)
  @inbounds for j in eachindex(locations)
    d = locations[j] - x
    w = _kernel_weight(kernel, locations[j], x, h) * weights[j]
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
# M[i,j] = moments[i+j-1] (the Hankel matrix of T_0 … T_{2d}) against `rhs` (= P_0 … P_d), returning
# the jet β = (β_0, …, β_d). The design matrix is positive semidefinite; a degenerate local design
# (singular — too little local data near x) yields an all-NaN jet, the caller's signal that x is
# unsupported.
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
    local_fit(locations, mass, response, x, h, kernel, degree)

The **local-polynomial jet** of a binned 1-D measure at the query point `x`: the kernel-weighted
least-squares fit of a degree-``d`` polynomial in ``u = \\ell_j - x`` to the bin data, returning the
coefficient vector (jet) ``\\beta = (\\beta_0, \\dots, \\beta_d)`` (with ``d`` the `degree` argument).
From bin locations ``\\ell_j``, design weights (counts) ``w_j`` (`mass`) and response totals ``r_j``
(`response`), and the bandwidth-scaled kernel ``K_h(u) \\propto K(u/h)``, form the design and response
moments

```math
\\begin{aligned}
T_\\ell &= \\sum_j K_h(\\ell_j - x)\\, w_j\\, (\\ell_j - x)^\\ell, \\\\
P_\\ell &= \\sum_j K_h(\\ell_j - x)\\, r_j\\, (\\ell_j - x)^\\ell,
\\end{aligned}
```

and solve the normal equations ``[\\,T_{i+\\ell}\\,]_{i,\\ell=0}^{d}\\,\\beta = (P_0, \\dots, P_d)``. The
jet is the local Taylor expansion: ``\\beta_0`` is the fitted value, ``\\beta_1`` the slope, and
``2\\beta_2`` the curvature — so a regression fit gives ``\\hat m(x) = \\beta_0`` and
``\\hat m''(x) = 2\\beta_2``. The weight is the shared `Kernels._kernel_weight` atom (its normalization
cancels in the solve), so any **symmetric** callable `kernel` works at any `degree`; an **asymmetric**
boundary kernel has no jet and is admitted only at `degree = 0` (the Nadaraya–Watson ratio), rejected
above it. A degenerate local design (too little data near `x`) returns an all-NaN jet.

This is the shared linear primitive of the local smoothers:
``\\hat m(x)`` is ``\\beta_0`` of one `local_fit`, and the local-likelihood density reuses the same
kernel-weighted-moment solve as the iteratively reweighted (IRLS) inner step of its local-Poisson fit.

Arguments:

- `locations` — the bin locations ``\\ell_j`` (geometric centers or recovered centroids).
- `mass` — the per-bin design weight ``w_j`` (the bin counts).
- `response` — the per-bin response total ``r_j`` (e.g. ``\\sum Y`` over the bin); pass `mass` to fit
  the design itself.
- `x` — the query point at which the jet is evaluated.
- `h` — the smoothing bandwidth (``h > 0``); it scales the kernel argument.
- `kernel` — the smoothing kernel: a symmetric callable ``K`` (admitted at any `degree`), or an
  `AsymmetricKernel` (`degree = 0` only). Its normalization cancels in the solve.
- `degree` — the local-polynomial order ``d``: `0` (local-constant / Nadaraya–Watson), `1`
  (local-linear), or `2` (local-quadratic); the returned jet then has `degree + 1` entries (an asymmetric
  `kernel` restricts it to `0`).
"""
function local_fit(locations, mass, response, x, h, kernel, degree)
  (kernel isa AsymmetricKernel && degree ≥ 1) && throw(ArgumentError(
    "an asymmetric kernel is a degree-0 boundary kernel — it has no jet; `degree` must be 0 (got $degree)"))
  T = _design_moments(locations, mass, x, float(h), kernel, 2 * degree)
  P = _design_moments(locations, response, x, float(h), kernel, degree)
  return _hankel_solve(T, P)
end


end # module LocalPolynomial
