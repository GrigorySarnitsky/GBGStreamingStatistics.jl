# Formalization

The typed statements behind the prose. Each operation gets one entry stating its **estimand**, its
**estimator**, and the **consistency claim** — the bias and variance with every normalization
constant pinned — each linked to the test that verifies the convergence rate numerically. The chain
is *prose ↔ formalization ↔ test*: the prose (docstrings, `ARCHITECTURE.md`) agrees with the typed
statement, the statement is verified against the code, and the linked test enforces it.

A nontrivial claim with **no entry here**, or an entry with **no test link**, is a recorded coverage
gap (`CLAUDE.md`, *Definition of done*).

Throughout, ``K`` is the smoothing kernel with ``\int K = 1``, ``\mu_2(K) = \int u^2 K(u)\,\mathrm
du`` its second moment, and ``R(K) = \int K^2`` its roughness; ``h`` is the bandwidth, ``n`` the
number of samples, ``f_X`` the marginal density of the conditioning variable.

The implemented entries are test-linked; a forthcoming operation gets its entry — in the same shape
(estimand / estimator / consistency / test) — when it lands, its estimand recorded in `ARCHITECTURE.md`
§4 meanwhile (the Pawula / Markov diagnostics and the higher Kramers–Moyal coefficients).

## Per-bin power-sums about the bin center (`Histograms`, the summaries)

- **Object.** per bin, the count and the **power-sums about the fixed bin center** ``S_p =
  \sum_{i\in\text{bin}}(X_i - c_j)^p`` for ``p = 1 \dots \texttt{order}`` (and order-2 bivariate ``S_{rs}``
  per cell for `Histogram2D`). Intrinsically finite-sample — no population estimand.
- **Accumulation.** an atomic scatter-add (`observe!` / `observe_pairs!`): one work-item per sample
  adds ``(X_i - c_j)^p`` to its bin's ``S_p``, in the accumulator type (`Float64` by default). `merge!` is elementwise addition (same
  binning ``\Rightarrow`` same centers). The local offsets are well-conditioned, so no compensated
  summation is used.
- **Mergeability (a monoid morphism).** Write ``(\mathcal{M}, \uplus, \varnothing)`` — finite sample
  collections (multisets over the sample space) under disjoint union, ``\varnothing`` the empty
  collection — and ``(\mathcal{S}, \oplus, 0)`` — summaries under `merge!` (elementwise addition of
  counts and power-sums), ``0`` the all-zero summary. Both are **commutative** monoids. Accumulation
  ``T : \mathcal{M} \to \mathcal{S}`` is the unique homomorphism extending the single-sample scatter
  ``t`` (one sample adds ``1`` to its bin's count and ``(X_i - c_j)^p`` to that bin's ``S_p``), so
  ``T(A) = \bigoplus_{x \in A} t(x)`` and the square

  ```math
  \begin{CD}
  \mathcal{M} \times \mathcal{M}  @>{\uplus}>>  \mathcal{M} \\
  @V{T \times T}VV                              @VV{T}V \\
  \mathcal{S} \times \mathcal{S}  @>{\oplus}>>  \mathcal{S}
  \end{CD}
  ```

  commutes: ``T(A \uplus B) = T(A) \oplus T(B)``, which **is** *accumulate-then-merge* ``=``
  *accumulate-together*. Commutativity and associativity of ``\oplus`` then make any reduction order
  or tree (multi-thread, multi-device, time-chunk) yield the same summary. The morphism is the
  **unnormalized** summary only: the normalize-and-smooth readout ``E_h`` (the ``1/n``, the shift
  formula, the kernel) that turns a summary into estimates is **not** a homomorphism — which is exactly
  why `merge!` has no population estimand (pool summaries first, estimate once, never the reverse).
  Exact for the integer counts; the `Float64` power-sums merge up to round-off (non-associative
  floating-point addition), so the test asserts ``=`` on counts and ``\approx`` on the ``S_p``.
- **Recovery.** the global central moments are the **shift formula**

  ```math
  \begin{aligned}
  \sum_i (X_i - \bar X)^p &= \sum_j \sum_{k=0}^{p} \binom{p}{k}\,(c_j - \bar X)^{p-k}\, S_k^{(j)}, \\
  \bar X &= \tfrac1n \sum_j \bigl(c_j n_j + S_1^{(j)}\bigr), \qquad S_0^{(j)} = n_j,
  \end{aligned}
  ```

  which is **exact** (the ``S_p`` use the actual ``X_i``, not the center) and numerically stable for
  the even ``p`` the consumers need; the recovered moments are exposed publicly as `mean`, `var`,
  `std`, and `central_moment` — `var` and `std` apply Bessel's ``n-1`` correction by default (matching
  `Statistics`), `corrected = false` the ``1/n`` form.
- **Test.** `test/histograms.jl` (recovered ``M_2 \dots M_6`` match the direct sum to ``10^{-8}`` and
  the true Gaussian moments to sampling error); `test/histograms2d.jl` (2D marginals + covariance);
  `test/mergeability.jl` (merge ``=`` together); `test/allocations.jl` (bounded `observe!` / `merge!`);
  `test/gpu/suite.jl` (CPU-vs-AMD-GPU parity).

## The cone — `scale!`, weighted `merge!`, hyperstatistics (§10)

- **Object.** the summary monoid ``(\mathcal{S}, \oplus, 0)`` of the previous entry, equipped also with a
  nonnegative scalar action ``w \cdot T`` (`scale!` — multiply every count and power-sum by ``w \ge 0``),
  making ``\mathcal{S}`` a **convex cone**: closed under ``\oplus`` and under nonnegative scaling, hence
  under all nonnegative linear combinations ``\sum_k w_k T_k`` (the weighted
  `merge!(into, from; weight)`, ``\texttt{into} \mathrel{+}= \texttt{weight} \cdot \texttt{from}``; the
  default ``\texttt{weight} = 1`` the plain ``\oplus``). Carried on a **Float-count variant**
  (`count = Float64`): a weighted count ``w\,n`` is not an integer.
- **Estimand (hyperstatistics; Beck & Cohen 2003).** the parameter-mixture (superstatistical) density over
  a slowly varying parameter ``\beta`` with mixing law ``f(\beta)``,

  ```math
  P(x) = \int f(\beta)\,P(x \mid \beta)\,\mathrm{d}\beta \;\approx\; \sum_k w_k\, f_k(x),
  \qquad w_k = f(\beta_k)\,\Delta\beta,
  ```

  with ``f_k = P(\cdot \mid \beta_k)`` the per-parameter density.
- **Estimator.** keep one summary ``T_k`` per ``\beta_k`` (accumulated independently, ``N_k`` samples),
  form the cone combination ``M = \sum_k (w_k / N_k)\, T_k`` by weighted `merge!`, and read out any
  smoothing functional ``E_h`` of ``M``. Since ``E_h`` normalizes by the total mass
  ``N_M = \sum_j m_j = \sum_k w_k``, this yields ``E_h(M) = (\sum_k w_k f_k(x)) / \sum_k w_k`` — the
  mixture (exactly, when the ``w_k`` sum to one); the ``1/N_k`` normalizes each component to a probability
  measure before mixing, and the conditional / transition estimators carry the local weight ``w_k f_k(x)``
  automatically (ratios of linear functionals). **Post hoc**: ``f(\beta)`` is chosen after the runs,
  re-weighting without re-accumulating.
- **Cone laws (the spec; the operations have no population estimand).** `scale!` is the monoid action —
  ``a \cdot (b \cdot T) = (ab) \cdot T``, ``1 \cdot T = T`` — and distributes over ``\oplus``:
  ``w \cdot (A \oplus B) = (w \cdot A) \oplus (w \cdot B)``; the weighted `merge!` is ``A \oplus (w \cdot B)``.
  `scale!` preserves the *normalized* distribution: the mean and the uncorrected central moments are exact
  invariants; only the total mass — and through it the Bessel-corrected variance — scales.
- **Exponential forgetting.** for a non-stationary stream, `scale!(s, λ)` with ``0 < \lambda < 1`` before
  folding the next block gives a geometrically down-weighted summary (effective window ``1/(1-\lambda)``);
  the recovered mean is the exact ``\lambda``-weighted mean of the blocks.
- **Caveat.** ``N_M = \sum_k w_k`` (``\approx 1`` for a normalized law), so ``M``'s *default* bandwidth
  reads ``n \approx 1``; pass an explicit ``h`` (or pre-scale the weights to a representative sample size).
  An effective-sample-size rule (Kish ``(\sum_i w_i)^2 / \sum_i w_i^2``) is a recorded refinement.
- **Test.** `test/hyperstatistics.jl` — the cone laws (`scale!` homogeneity / distributivity, weighted
  merge ``=`` scale-then-merge, ``\texttt{weight} = 1`` ``=`` plain `merge!`, the Float-count guard); the
  discrete Gaussian-mixture oracle; the **χ²-superstatistics → Student-``t``** oracle
  (``X \mid \beta \sim N(0, 1/\beta)``, ``\beta \sim \mathrm{Gamma}(\nu/2, 2/\nu)`` ``\Rightarrow`` the
  Student-``t_\nu`` marginal, recovered into the heavy tail); exact exponential forgetting; the
  `Histogram2D` cone with `conditional_mean` scale-invariance. `test/allocations.jl` (allocation-free
  `scale!` / weighted `merge!`); `test/gpu/suite.jl` (CPU-vs-GPU parity, incl. the atomic-float add of
  Float-count `observe!`).

## Marginal density (§4.1) — `density`

- **Estimand.** ``f_X = \mathrm{d}\mu_X/\mathrm{d}\lambda``, the Lebesgue density of ``X``.
- **Estimator.** the binned KDE ``\hat f_h(x) = \frac{1}{nh}\sum_j n_j K\!\left(\frac{x-\ell_j}{h}\right)``
  (``\int K = 1``, Gaussian ``K``) from the bin counts ``n_j`` at bin locations ``\ell_j`` — the within-bin
  **centroid** ``\ell_j = c_j + S_1^{(j)}/n_j`` (linear binning, the default `binning = :linear`) or the
  geometric center ``\ell_j = c_j`` (`binning = :simple`) — with ``h`` from `bandwidth`
  (Silverman ``h = (4/3)^{1/5}\hat\sigma\,n^{-1/5}``).
- **Consistency claim** (``h \to 0``, ``nh \to \infty``): ``\operatorname{bias}\hat f_h(x) \approx
  \tfrac12 h^2 \mu_2(K)\, f_X''(x)`` and ``\operatorname{Var}\hat f_h(x) \approx R(K)\,f_X(x)/(n h)``;
  the binning error converges to the exact KDE at ``O(\Delta^2)`` in the bin width ``\Delta`` for linear
  binning — the centroid matches the within-bin first moment, ``\sum_{i\in j}(X_i - \ell_j) = 0``, so the
  leading ``O(\Delta)`` discretization term vanishes — and at ``O(\Delta)`` for simple binning; either is
  dominated by the ``O(h^2)`` smoothing bias when ``\Delta`` is small.
- **Test.** `test/oracles.jl` (recovers the ``N(0,1)`` pdf to ``< 3\times10^{-3}`` over ``[-3,3]``,
  integrates to 1). `test/density_convergence.jl` — a representative analytic-PDF suite (Gaussian,
  uniform, exponential, Student-``t`` ``\nu=3`` heavy tail, log-normal, bimodal: each recovered to
  ``\sim 10^{-2}`` and integrating to ``\approx 1``, the variable bandwidth tracking the ``t`` heavy
  tail better than the fixed) and the **convergence-rate test**: the ISE falls at the theoretical
  ``n^{-2p/(2p+1)}`` — slope ``\approx -0.80`` for degree 0/1 (``p=2``), the steeper ``\approx -0.87``
  (``n^{-8/9}``, ``p=4``) for degree 2, the Abramson variable bandwidth sharing the degree-0 order; plus the
  **binning-error convergence** sub-test — the binned estimate approaches the exact KDE at ``O(\Delta^2)``
  for the default linear/centroid `binning` (slope ``\approx 2`` in ``\Delta``) and ``O(\Delta)`` for
  `binning = :simple` (slope ``\approx 1``).

## Local-log-polynomial & adaptive density (§4.1) — `density`'s `degree` / `adaptive`

- **Estimand.** ``f_X``, exactly as for the degree-0 `density` above; the knobs change the *estimator*,
  not the target.
- **`adaptive` (Abramson variable bandwidth).** A fixed-bandwidth pilot ``\tilde f`` at ``h`` sets
  per-bin bandwidths ``h_j = h\,(\tilde f(c_j)/g)^{-1/2}`` (``g`` the geometric mean of ``\tilde f`` over
  the sample), widening the kernel in the sparse tails: ``\operatorname{bias} = O(h^4)`` at interior
  smooth ``x`` (Abramson 1982), with far less tail variance than a fixed ``h``.
- **`degree` ``\ge 1`` (local-log-polynomial).** Fit ``\log f`` by a degree-``d`` polynomial through a binned
  **local Poisson** likelihood — ``\log \mathbb{E}[n_j] = \log(nb) + \sum_{l \le d} \beta_l (c_j-x)^l``,
  weights ``K((c_j-x)/h_j)``, by Newton/IRLS — and report ``\hat f(x) = e^{\hat\beta_0}`` (Loader 1996).
  Degree 1 is exact for an exponential density (tracks exponential tails), degree 2 for a Gaussian. The
  two knobs compose. On a **bounded support** the degree-0 KDE accepts `bounds = (a, b)`, renormalizing
  by the kernel mass ``\int_a^b K_h(x-u)\,\mathrm{d}u`` inside the support (Jones 1993) to remove the
  ``O(1)`` boundary bias that otherwise halves the estimate at the edge.
- **Test.** `test/marginaldensity.jl` (every grid corner recovers the ``N(0,1)`` pdf into the tail and
  integrates to ``\approx 1``; the variable bandwidth resolves the sharp Laplace cusp the fixed KDE
  oversmooths; degree-1 beats the KDE in a Laplace tail).

## Asymmetric kernels (§4.1) — `density`'s `gamma_kernel` / `beta_kernel`

- **Estimand.** ``f_X`` on a bounded or positive support.
- **Estimator.** the varying-shape estimate ``\hat f(x) = \tfrac1n \sum_j n_j K(c_j;\,x,\,b)`` with ``K(\cdot;x,b)`` the modified ``\mathrm{Gamma}(\rho(x),\,b)`` density (gamma kernel, support ``[0,\infty)``; Chen 2000) or the modified ``\mathrm{Beta}`` density (beta kernel, support ``[0,1]``; Chen 1999) — Chen's bias-reduced ``f_2`` shapes (interior shape ``x/b``, kernel mean ``\approx x``, removing the ``f'`` bias of the naive ``x/b+1`` form); the support matches the data, so no mass leaks past the boundary. A reduction over the bin counts (mergeable); default ``b = \hat\sigma\,n^{-2/5}``.
- **Consistency claim.** bias ``O(b)`` uniformly *including at the boundary* (where a symmetric kernel is ``O(h)`` and inconsistent); interior variance ``O(n^{-1}b^{-1/2})``; optimal ``b \sim n^{-2/5}``, MISE ``\sim n^{-4/5}``. Not self-normalizing (an ``O(b)`` mass deviation). The modified ``f_2`` shape lowers the MISE constant (Chen proves ``\mathrm{MISE}(f_2) \le \mathrm{MISE}(f_1)``; only ``f_2`` ships, so the gain is the literature result, not a repo-reproducible number).
- **Test.** `test/marginaldensity.jl` (gamma recovers the exponential at the boundary, beating the symmetric KDE there; beta recovers ``\mathrm{Uniform}(0,1)`` at both edges; integrates to ``\approx 1``; mergeability; degree/adaptive/bounds rejected).

## Parametric start & normalization (§4.1) — `density`'s `start` / `normalize`

- **Estimand.** ``f_X``, exactly as above; `start` and `normalize` change the *estimator*, not the target.
- **Estimator (`start`, Hjort & Glad 1995).** The semiparametric estimate ``\hat f(x) = p(x;\hat\theta)\,\hat r_h(x)``, where ``p(\cdot;\hat\theta)`` is a parametric *start* fit by **method of moments** from the recovered ``\hat\mu``, ``\hat\sigma^2``, and ``\hat r_h(x) = \tfrac1{nh}\sum_j (n_j/p(c_j;\hat\theta))\,K((x-c_j)/h)`` is a KDE of the *correction* ``r = f_X/p`` (the start's local relative error). Five starts are provided — `normal_start`, `lognormal_start`, `gamma_start`, `exponential_start`, `inverse_gaussian_start` — each a density whose parameters are fixed by the first two moments, so ``\hat\theta`` is itself a post-hoc functional of the summary and the whole estimate is **mergeable** (``\hat\theta`` merges through the moments; the reweighted counts add). Degree-0, fixed-``h``, symmetric kernel.
- **Consistency claim.** For a fixed start with bounded, twice-differentiable correction ``r``, consistent for ``f_X`` at ``\operatorname{bias} = O(h^2)``, ``\operatorname{var} = O(1/(nh))`` (MISE ``\sim n^{-4/5}``), but with the **leading bias governed by ``r'' = (f_X/p)''`` in place of ``f_X''``** and the variance ``R(K)f_X(x)/(nh)`` unchanged to leading order. When the start is correct ``r \equiv 1``, ``r'' \equiv 0``: the leading bias vanishes and the estimate collapses to the parametric fit. So the start is a **bias reduction** that pays in the bias-dominated regime (skewed / heavy-tailed shapes, oversmoothed ``h``), not in the large-``n`` limit where the KDE is already nearly unbiased.
- **Estimator (`normalize`).** ``\hat f / Z`` with ``Z = \int \hat f`` over the summary's range intersected with the estimator's support, by the trapezoidal rule on a fine grid — a deterministic, mergeable functional of the estimate making ``\int \hat f = 1`` **exactly** (the asymmetric kernels and a parametric start are not self-normalizing; the symmetric KDE already is, to the grid's ``O(\Delta^2)``).
- **Test.** `test/marginaldensity.jl` (a normal start at a wide ``h`` recovers the Gaussian peak the plain KDE flattens by ``\approx 10\%``, ISE ``> 20\times`` smaller; a log-normal start tracks log-normal data ``> 3\times`` better than the symmetric KDE, and pointwise into the right tail; `normalize` drives ``\int\hat f`` to ``1`` for a start and for the gamma kernel; mergeability; the rejected knob combinations).

## Pointwise standard error (§4.1) — `density_standard_error`

- **Estimand.** the standard error of the density estimator, ``\mathrm{SE}(x) = \sqrt{\operatorname{Var}\hat f_h(x)}`` with the KDE asymptotic variance ``\operatorname{Var}\hat f_h(x) \approx R(K)\,f_X(x)/(nh)``.
- **Estimator.** ``\sqrt{R(K)\,\hat f_h(x)/(nh)}`` from the recovered ``\hat f_h``, ``n``, ``h``, and ``R(K) = \int K^2`` (by quadrature). The ``z``-band ``\hat f_h(x) \pm z\,\mathrm{SE}(x)`` is a variability band about ``\mathbb{E}\hat f_h``, covering ``f_X`` where the ``O(h^2)`` smoothing bias is small. Degree 0/1 (the variable bandwidth and degree ``\ge 2`` are a coverage gap).
- **Test.** `test/marginaldensity.jl` (the SE equals ``\sqrt{R(K)\hat f/(nh)}`` exactly; the 95\% band covers the ``N(0,1)`` truth over ``\ge 75\%`` of the interior; the SE shrinks with ``n``; vector / empty / error forms).

## Bandwidth selection (§4.5) — `bandwidth`

- **Estimand.** the AMISE-optimal bandwidth ``h^\star = \arg\min_h \int \mathbb{E}[(\hat f_h(x) - f_X(x))^2]\,\mathrm{d}x`` for a kernel ``K`` and local-polynomial degree ``d`` — the asymptotically optimal smoothing scale, itself ``O(n^{-1/(2p+1)})`` (a finite-sample selector, not a fixed population constant).
- **Estimator** (`rule = :normal`, the kernel- and degree-aware normal-reference). ``h = C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}`` with bias order ``p = 2`` (degree 0/1) or ``p = 4`` (degree 2/3), and constant ``C(K,d) = [W_{ikk}(K,d)\,2^{p+1}\sqrt{\pi}/(2p-1)!!]^{1/(2p+1)}``. The canonical factor is ``W_{ikk}(K,d) = R(K)/\mu_2(K)^2`` at degree 0/1, and the fourth-order equivalent-kernel analog ``W_{ikk} = (4!^2/8)\,R(K^\star)/\mu_4(K^\star)^2`` with ``K^\star = K\,(n_4 - v^2 n_2)/(n_0 n_4 - n_2^2)`` at degree 2/3 (Loader 1996); here ``R(K) = \int K^2``, ``\mu_2(K) = \int u^2 K`` (by quadrature) and ``n_k = \int v^k K``. ``\hat\sigma`` and ``n`` are read from the one summary.
- **Estimator** (`rule = :sj`, the data-driven Sheather–Jones plug-in; degree 0/1). The curvature ``\psi_4 = R(f'')`` and its ``\psi_6`` pilot are estimated from the bin counts, ``\hat\psi_r(g) = (n^2 g^{r+1})^{-1}\sum_{j,k} n_j n_k\,\varphi^{(r)}((c_j - c_k)/g)`` (Gaussian ``\varphi``, the ``n`` self-pairs removed), and the solve-the-equation bandwidth (Sheather & Jones 1991) is the root of ``h = [R(\varphi)/(\mu_2(\varphi)^2\,\hat\psi_4(\alpha\,h^{5/7})\,n)]^{1/5}``; the eval kernel then multiplies by the canonical-bandwidth ratio ``\delta(K)/\delta(\varphi)`` with ``\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}``. Both functionals are reductions of the bin counts, hence mergeable; degree ``\ge 2`` is a recorded coverage gap. The **improved** variant (`rule = :isj`; Botev, Grotowski & Kroese 2010) instead reads ``\|f^{(s)}\|^2`` off the DCT of the bin counts and solves the AMISE fixed point ``t = \xi\gamma^5(t)`` — also a mergeable functional of the counts, and robust on strongly multimodal data.
- **Estimator** (`rule = :aic`, the fit-based information criterion; degree 0/1). Minimise ``\mathrm{AIC}(h) = -2\sum_j n_j\log\hat f_h(c_j) + 2\,\mathrm{df}(h)`` of the binned KDE over a log-spaced ``h``-grid, with ``\mathrm{df}(h) = (K(0)/h)\sum_j n_j/(n\hat f_h(c_j))``; both terms diverge as ``h\to0``, giving a stable interior optimum. Mergeable (a functional of the counts); the eval kernel multiplies by ``\delta(K)/\delta(\varphi)``. The leave-one-out cross-validation selectors (likelihood-CV, least-squares LSCV) are a recorded coverage gap — binning ties drive their objective to ``h\to0``.
- **Estimator** (`rule = :plugin`, the **asymmetric-kernel** Chen plug-in; `gamma_kernel` / `beta_kernel`). The gamma / beta estimators take a scale ``b`` (rate ``n^{-2/5}``), not ``h``; their AMISE (Chen 1999, 2000; modified ``f_2``) is ``A b^2 + C b^{-1/2}/n`` with ``A = \tfrac14\int w^2 (f'')^2``, ``C = (2\sqrt\pi)^{-1}\mathbb{E}[w(X)^{-1/2}]`` and local-scale weight ``w(x) = x`` (gamma) or ``x(1-x)`` (beta), giving ``b^\star = [\mathbb{E}[w(X)^{-1/2}]/(2\sqrt\pi\,n\int w^2 (f'')^2)]^{2/5}``. ``\mathbb{E}[w(X)^{-1/2}]`` is an exact count moment; the curvature ``\int w^2 (f'')^2`` is a grid quadrature of the pilot KDE's ``\hat f''_g`` (pilot ``g = 1.3\,s\,n^{-1/9}`` at the robust scale ``s = \min(\hat\sigma,\,\mathrm{IQR}/1.349)``, the variance inflation ``R(\varphi'')\,\mathbb{E}[w^2]/(n g^5)`` removed), which — unlike the integration-by-parts ``\psi_4`` U-statistic — stays valid where ``f(0) > 0`` (the ``w \to 0`` boundary weight suppresses the pilot's edge leakage). A deterministic functional of the counts (mergeable). The symmetric `:sj` / `:isj` / `:aic` do not apply; a fully data-adaptive curvature pilot is a recorded refinement.
- **Consistency claim.** with ``\hat\sigma \to \sigma`` and the Gaussian reference ``R(f^{(p)}) = (2p-1)!!/(2^{p+1}\sqrt{\pi}\,\sigma^{2p+1})``, ``C(K,d)`` is the exact AMISE-optimal constant: degree 0/1 gives ``h \sim n^{-1/5}`` (MISE ``\sim n^{-4/5}``), degree 2 the higher-order ``h \sim n^{-1/9}`` (MISE ``\sim n^{-8/9}``); the Gaussian/degree-0 case is Silverman ``(4/3)^{1/5}``, the Epanechnikov ``\approx 2.34``.
- **Test.** `test/bandwidth.jl` (Gaussian degree 0 reproduces Silverman, degree 1 equals degree 0, the Epanechnikov ``h`` is ``2.214\times`` the Gaussian, degree 2 matches the closed-form ``n^{-1/9}`` rule and exceeds degree 0; coverage-gap and sample-size errors; plus the Sheather–Jones tests — SJ ``\approx`` the normal reference on Gaussian data, much tighter on a bimodal mixture, the exact ``\delta(K)`` kernel ratio, mergeability, the degree ``\ge 2`` gap; and the AIC tests — a stable optimum below the normal reference on Gaussian data, far smaller on a bimodal mixture, the exact ``\delta(K)`` ratio, mergeability, and the deliberately-absent CV rules erroring; and the asymmetric Chen plug-in — `:normal` reproduces the ``\hat\sigma\,n^{-2/5}`` rule of thumb, `:plugin` lowers the ISE (``\ge 1.25\times`` on log-normal, ``\ge 1.11\times`` on gamma; ``\approx 1.4\text{–}1.7\times`` typically, seed-audited) and stays ``\approx`` neutral on the near-exponential boundary case, gamma and beta kernels, mergeability, and the inapplicable `:sj`/`:isj`/`:aic` and `degree` ``\ge 1`` rules erroring; the Improved Sheather–Jones tests (`:isj` ``\approx`` the normal reference on Gaussian, ``< 0.5\times`` it on a bimodal mixture, the exact ``\delta(K)`` ratio, mergeability); and the **Marron–Wand claw** stress, where the data-driven `:sj` / `:isj` resolve the five spikes the normal reference oversmooths — ``\approx`` an order of magnitude lower ISE) and `test/oracles.jl` (the kernel-aware ratio on the ``N(0,1)`` oracle).

## Local-polynomial jet primitive (`LocalPolynomial`)

- **Object.** the kernel-weighted local-polynomial **jet** at a query point ``x``: the degree-``d``
  weighted-least-squares fit of a polynomial in ``u = \ell_j - x`` to a binned 1-D measure (locations
  ``\ell_j``, design weights ``w_j``, response totals ``r_j``), returning ``\beta = (\beta_0, \dots,
  \beta_d)`` — ``\beta_0`` the fitted value, ``\beta_1`` the slope, ``2\beta_2`` the curvature. A
  deterministic primitive, not a population estimand; it is the shared linear core of the local smoothers
  (`docs/ARCHITECTURE.md` §11).
- **Computation.** form the design and response moments ``T_\ell = \sum_j K_h(\ell_j - x)\,w_j\,(\ell_j -
  x)^\ell`` (``\ell = 0 \dots 2d``) and ``P_\ell = \sum_j K_h(\ell_j - x)\,r_j\,(\ell_j - x)^\ell``
  (``\ell = 0 \dots d``), and solve the ``(d+1)\times(d+1)`` Hankel normal equations
  ``[\,T_{i+\ell}\,]\,\beta = (P_0, \dots, P_d)`` (a small linear solve, never an explicit inverse). The
  kernel normalization cancels in the solve; a degenerate (singular) design returns an all-NaN jet.
  `conditional_mean` is ``\beta_0`` of one fit; the local-likelihood `density` reuses the same solve as
  its iteratively reweighted (IRLS) inner step.
- **Test.** `test/localpolynomial.jl` — on a binned measure whose per-bin response mean is exactly a
  polynomial ``P``, the jet recovers ``(P(x),\,P'(x),\,P''(x)/2)`` to ``10^{-6}`` (degree 2) and
  ``(P,\,P')`` (degree 1); degree 0 is the Nadaraya–Watson ratio ``P_0/T_0``; a degenerate design gives an
  all-NaN jet; the kernel normalization cancels (Gaussian vs unnormalized agree).

## Local-linear conditional mean (§4.2) — `conditional_mean`

- **Estimand.** ``m(x) = \mathbb{E}[Y \mid X = x] = \int y\, f_{Y\mid X}(y \mid x)\,\mathrm{d}y``,
  the ``L^2(\sigma(X))`` projection of the response onto the conditioning variable.
- **Estimator.** the local-linear fit ``\hat m_h(x) = (T_2 P_0 - T_1 P_1)/(T_0 T_2 - T_1^2)`` from
  the kernel-weighted power-sums ``T_\ell = \sum_i K_h(X_i-x)(X_i-x)^\ell`` and
  ``P_\ell = \sum_i K_h(X_i-x)(X_i-x)^\ell Y_i`` — recovered from the joint `Histogram2D`'s per-X-column
  counts and response totals ``\sum Y`` (exact from the per-cell power-sums).
- **Consistency claim** (interior ``x``, ``h \to 0``, ``nh \to \infty``):

  ```math
  \operatorname{bias}\hat m_h(x) \approx \tfrac12 h^2 \mu_2(K)\, m''(x),
  \qquad
  \operatorname{Var}\hat m_h(x) \approx \frac{R(K)\,\sigma^2(x)}{n\,h\,f_X(x)},
  ```

  giving the optimal rate ``\hat m_{h^\star}(x) - m(x) = O_p(n^{-2/5})`` at ``h^\star \sim n^{-1/5}``.
- **Test.** `test/regression.jl` — the bivariate-Gaussian oracle ``Y = \rho X + \sqrt{1-\rho^2}\,
  \varepsilon`` (so ``\mathbb{E}[Y\mid X=x] = \rho x``): `conditional_mean` recovers the slope to
  ``\lesssim 2\times10^{-2}``, with the `degree` (0/1/2) and `kernel` knobs exercised (degree-2 and
  Epanechnikov both recover the linear mean). The ``n^{-2/5}`` rate test is scheduled.

## Conditional variance (§4.2) — `conditional_variance`

- **Estimand.** ``\sigma^2(x) = \operatorname{Var}(Y \mid X = x) = \mathbb{E}[Y^2 \mid X = x] - m(x)^2``,
  with ``m(x)`` the conditional mean above.
- **Estimator** (`method = :residual`, the default; Fan & Yao 1998). The **two-stage residual-based**
  estimator: fit ``\hat m`` at each X-column center ``c^x_j``, form the per-column squared-residual total
  ``\sum (Y - \hat m_j)^2 = \sum Y^2 - 2\hat m_j \sum Y + \hat m_j^2 N`` (``\hat m_j = \hat m(c^x_j)``,
  exact from the stored column sums), and take a second local fit of it,
  ``\hat\sigma^2(x) = \widehat{\mathbb{E}}[(Y - \hat m(X))^2 \mid x]``. The `:plugin` alternative is the
  simple difference ``\widehat{\mathbb{E}}[Y^2 \mid x] - \hat m_h(x)^2``; both clamped non-negative.
- **Consistency claim.** both are consistent for ``\sigma^2(x)``. The plug-in carries the local-polynomial
  smoothing bias of the curved ``\mathbb{E}[Y^2\mid x]`` — an ``O(h^2)`` term driven by ``m''`` and the
  curvature of ``\mathbb{E}[Y^2\mid x]``, large where ``m`` curves. The residual-based estimator fits the
  far flatter ``\mathbb{E}[(Y - m)^2\mid x] = \sigma^2(x)``, whose ``O(h^2)`` bias ``\propto (\sigma^2)''``
  is much smaller, so it removes the ``m''``-driven bias.
- **Test.** `test/regression.jl` — recovers ``1-\rho^2`` on the (linear-mean) bivariate-Gaussian oracle to
  a few percent, and on a **heteroscedastic curved-mean** oracle (``X\sim N(0,1)``, ``m(x)=x^2``,
  ``\sigma^2(x)=0.2+0.1x^2``) the residual estimator recovers ``\sigma^2(x)`` (degree 2) and beats the
  plug-in in the curved tail (the plug-in inflated ``>1.5\times`` by the ``m''`` curvature); the universal
  contract via `check_point_estimator` (scalar ``=`` vector, functional of the merged summary,
  ``\mathrm{NaN}`` off-support). The ``n^{-2/5}`` rate test is scheduled.

## Response-aware local bandwidth (§4.5) — `bandwidth(::Histogram2D, x; rule = :local)`

- **Estimand.** the AMISE-optimal **local** bandwidth for the conditional mean, balancing the local-linear
  squared bias against the variance at each ``x``:

  ```math
  h^\star(x) = \left[\frac{R(K)\,\sigma^2(x)}{\mu_2(K)^2\,m''(x)^2\,f_X(x)\,n}\right]^{1/5}.
  ```

- **Estimator.** plug in the three ``x``-local quantities from the one joint summary, each at a pilot
  bandwidth: ``\hat\sigma^2(x)`` (the Fan–Yao `conditional_variance`), the curvature
  ``\hat m''(x) = 2\beta_2`` (the degree-2 `LocalPolynomial` jet), and ``\hat f_X(x)`` (the X-marginal
  KDE); ``R(K)``, ``\mu_2(K)`` by quadrature. The default pilot is the curvature-rate rule of thumb
  ``\hat\sigma_X\,n^{-1/9}`` (wider than ``h^\star``, so the second derivative is stable); where the
  regression is locally linear (``m'' \to 0``) the optimal ``h`` is unbounded, capped at the X-range.
- **Consistency claim.** with consistent pilots ``\hat\sigma^2 \to \sigma^2``, ``\hat m'' \to m''``,
  ``\hat f_X \to f_X``, the plug-in is consistent for ``h^\star(x) \sim n^{-1/5}``; each pilot enters
  through a weak (``1/5``-power) factor, so ``h^\star`` is robust to pilot error.
- **Test.** `test/localbandwidth.jl` — the quadratic-mean Gaussian oracle (``X\sim N(0,1)``,
  ``Y = X^2 + \sigma\varepsilon``: ``m''=2``, ``\sigma^2`` constant, ``f_X = \varphi``) where
  ``h^\star(x)`` is closed-form: recovered near the centre to ``\lesssim 12\%``, positive/finite on
  support, ``\mathrm{NaN}`` off it, scalar ``=`` vector, default-pilot sanity, and mergeability.

## Conditional increment moments, drift & diffusion (§4.4) — `increment_moment`, `kramers_moyal`, `drift`, `diffusion`

- **Estimand.** the conditional increment moment and the Kramers–Moyal coefficient

  ```math
  M_n(x, \tau) = \mathbb{E}\!\left[(X_{t+\tau} - X_t)^n \mid X_t = x\right], \qquad
  D^{(n)}(x) = \frac{1}{n!}\lim_{\tau\to0} \frac{M_n(x,\tau)}{\tau},
  ```

  with `kramers_moyal` ``\equiv D^{(n)}``, `drift` ``\equiv D^{(1)}``, and `diffusion` ``\equiv D^{(2)}``. Model-neutral — a conditional
  expectation of an increment power, no diffusion premise.
- **Estimator.** conditioning on ``X_t = x`` (so the increment is ``\Delta X = Y - x`` with
  ``Y = X_{t+\tau}``), the order-2 joint summary resolves

  ```math
  M_1(x) = \hat m_h(x) - x, \quad
  M_2(x) = \widehat{\mathbb{E}}[Y^2\mid x] - 2x\,\hat m_h(x) + x^2, \quad
  \widehat D^{(1)} = \frac{M_1}{\tau}, \quad \widehat D^{(2)} = \frac{M_2}{2\tau},
  ```

  where ``\hat m_h`` and ``\widehat{\mathbb{E}}[Y^2\mid x]`` are the local-linear §4.2 fits (exact from
  the per-cell power-sums). ``n \ge 3`` is a recorded coverage gap — it needs order-``n`` per-cell
  ``Y`` power-sums.
- **Consistency claim.** ``\widehat D^{(1)}`` is consistent for the finite-``\tau`` drift at the
  `conditional_mean` rate, with no smoothing bias where ``m(x)`` is linear in ``x``; ``\widehat D^{(2)}``
  carries the local-polynomial bias from the curvature of the **raw** moment
  ``\mathbb{E}[(\Delta X)^2\mid x]`` (a higher `degree` reduces it; the Fan–Yao *variance* refinement does
  not apply — ``D^{(2)}`` is the raw increment moment, not the conditional variance). The finite-``\tau``
  coefficients approach ``D^{(n)}(x)`` as ``\tau\to0`` up to the ``O(\tau)`` finite-time bias (the
  Bell-polynomial correction, forthcoming).
- **Test.** `test/oracles.jl` — the Ornstein–Uhlenbeck oracle (exact stationary lagged pairs from
  ``\mathrm dX = -\theta X\,\mathrm dt + \sigma\,\mathrm dW``): `drift` recovers
  ``x(e^{-\theta\tau}-1)/\tau`` to ``\le 2\times10^{-2}`` and `diffusion` recovers the closed-form
  ``[s^2(1-e^{-2\theta\tau}) + x^2(e^{-\theta\tau}-1)^2]/(2\tau)`` (``s^2 = \sigma^2/2\theta``) to
  ``8\%``. `test/incrementmoments.jl` — the variance decomposition
  ``M_2 = \operatorname{Var}(Y\mid x) + M_1^2`` (against the matching `:plugin` variance — the exact
  algebraic identity ``M_2 - M_1^2 = \widehat{\mathbb{E}}[Y^2\mid x] - \hat m^2``), the ``1/(n!\,\tau)``
  scaling, the ``n\ge3`` and static-lag coverage-gap errors, and mergeability. The ``n^{-2/5}`` rate test
  is forthcoming.

## Conditional & transition density (§4.3) — `conditional_density`, `transition_density`

- **Estimand.** ``f_{Y\mid X}(y\mid x) = f_{X,Y}(x,y)/f_X(x)``, and the two-time transition density
  ``p_\tau(y\mid x) = f_{X_{t+\tau}\mid X_t}(y\mid x)`` (the same object on a lag-``\tau`` joint).
- **Estimator.** the binned double-kernel (Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996), with a
  symmetric `kernel` ``K`` (``\int K = 1``, default `gaussian_kernel`), conditioning bandwidth
  ``h_1 =`` `hx` and response bandwidth ``h_2 =`` `hy`: from
  the cell counts ``n_{jk}`` at centers ``(c^x_j, c^y_k)`` with ``N_j = \sum_k n_{jk}``,

  ```math
  \hat f(y\mid x) =
    \frac{\sum_{j,k} n_{jk}\, K_{h_1}(c^x_j - x)\, K_{h_2}(c^y_k - y)}{\sum_j N_j\, K_{h_1}(c^x_j - x)} ,
  ```

  computed as ``\hat f(y\mid x) = \big(\sum_k P_k\, K_{h_2}(c^y_k - y)\big) / (\sqrt{2\pi}\,h_2 D)`` with
  the X-kernel-weighted profile ``P_k = \sum_j n_{jk} K_{h_1}(c^x_j - x)`` and ``D = \sum_k P_k``. The
  X-kernel weighting cancels in the ratio; the surviving ``1/h_2`` makes ``\int \hat f(y\mid x)\,\mathrm dy = 1`` exactly (``\int K = 1``).
- **Consistency claim** (``h_1, h_2 \to 0``, ``n h_1 h_2 \to \infty``): the pointwise bias is ``\tfrac12 h_1^2 \mu_2(K)\big[\partial_{xx} f_{Y\mid X} + 2\,\partial_x f_{Y\mid X}\,\partial_x f_X / f_X\big] + \tfrac12 h_2^2 \mu_2(K)\,\partial_{yy} f_{Y\mid X}`` and the variance ``R(K)^2 f_{Y\mid X}/(n h_1 h_2 f_X)``, giving the slower conditional rate ``O(n^{-2/3})`` at the optimal ``h \sim n^{-1/6}`` (the bias is worst at the conditional mode).
- **Test.** `test/conditionaldensity.jl` — the bivariate-Gaussian fixture ``Y\mid X=x \sim N(\rho x, 1-\rho^2)``: ``\int \hat f\,\mathrm dy \approx 1`` and pointwise recovery to ``3\times10^{-2}``; the
  scalar/vector forms, the `transition_density` equivalence, off-support `NaN`, bandwidth checks, and
  mergeability. `test/oracles.jl` — the Ornstein–Uhlenbeck oracle: `transition_density` recovers the
  closed-form ``p_\tau(y\mid x) = N(x e^{-\theta\tau},\, s^2(1-e^{-2\theta\tau}))`` to ``4\times10^{-2}``
  (the ``O(h^2)`` mode bias). The ``n^{-2/3}`` rate test is forthcoming.
