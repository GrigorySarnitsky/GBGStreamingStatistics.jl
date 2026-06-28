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

The implemented entries are test-linked; a forthcoming operation gets its entry in the same shape
(estimand / estimator / consistency / test) before it lands.

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
  \sum_i (X_i - \bar X)^p = \sum_j \sum_{k=0}^{p} \binom{p}{k}\,(c_j - \bar X)^{p-k}\, S_k^{(j)},
  \qquad S_0^{(j)} = n_j, \quad \bar X = \tfrac1n \sum_j \bigl(c_j n_j + S_1^{(j)}\bigr),
  ```

  which is **exact** (the ``S_p`` use the actual ``X_i``, not the center) and numerically stable for
  the even ``p`` the consumers need; the recovered moments are exposed publicly as `mean`, `var`,
  `std`, and `central_moment` — `var` and `std` apply Bessel's ``n-1`` correction by default (matching
  `Statistics`), `corrected = false` the ``1/n`` form.
- **Test.** `test/histograms.jl` (recovered ``M_2 \dots M_6`` match the direct sum to ``10^{-8}`` and
  the true Gaussian moments to sampling error); `test/histograms2d.jl` (2D marginals + covariance);
  `test/mergeability.jl` (merge ``=`` together); `test/allocations.jl` (bounded `observe!` / `merge!`);
  `test/gpu/suite.jl` (CPU-vs-AMD-GPU parity).

## Marginal density (§4.1) — `density`

- **Estimand.** ``f_X = \mathrm{d}\mu_X/\mathrm{d}\lambda``, the Lebesgue density of ``X``.
- **Estimator.** the binned KDE ``\hat f_h(x) = \frac{1}{nh}\sum_j n_j K\!\left(\frac{x-c_j}{h}\right)``
  (``\int K = 1``, Gaussian ``K``) from the bin counts ``n_j`` at centers ``c_j``, with ``h`` from
  `bandwidth` (Silverman ``h = (4/3)^{1/5}\hat\sigma\,n^{-1/5}``).
- **Consistency claim** (``h \to 0``, ``nh \to \infty``): ``\operatorname{bias}\hat f_h(x) \approx
  \tfrac12 h^2 \mu_2(K)\, f_X''(x)`` and ``\operatorname{Var}\hat f_h(x) \approx R(K)\,f_X(x)/(n h)``;
  the binning error ``O(b^2)`` is dominated by the ``O(h^2)`` smoothing bias when ``b \ll h``.
- **Test.** `test/oracles.jl` (recovers the ``N(0,1)`` pdf to ``< 3\times10^{-3}`` over ``[-3,3]``,
  integrates to 1). `test/density_convergence.jl` — a representative analytic-PDF suite (Gaussian,
  uniform, exponential, Student-``t`` ``\nu=3`` heavy tail, log-normal, bimodal: each recovered to
  ``\sim 10^{-2}`` and integrating to ``\approx 1``, the variable bandwidth tracking the ``t`` heavy
  tail better than the fixed) and the **convergence-rate test**: the ISE falls at the theoretical
  ``n^{-2p/(2p+1)}`` — slope ``\approx -0.80`` for degree 0/1 (``p=2``), the steeper ``\approx -0.87``
  (``n^{-8/9}``, ``p=4``) for degree 2, the Abramson variable bandwidth sharing the degree-0 order.

## Local-log-polynomial & adaptive density (§4.1) — `density`'s `degree` / `adaptive`

- **Estimand.** ``f_X``, exactly as for the degree-0 `density` above; the knobs change the *estimator*,
  not the target.
- **`adaptive` (Abramson variable bandwidth).** A fixed-bandwidth pilot ``\tilde f`` at ``h`` sets
  per-bin bandwidths ``h_j = h\,(\tilde f(c_j)/g)^{-1/2}`` (``g`` the geometric mean of ``\tilde f`` over
  the sample), widening the kernel in the sparse tails: ``\operatorname{bias} = O(h^4)`` at interior
  smooth ``x`` (Abramson 1982), with far less tail variance than a fixed ``h``.
- **`degree ≥ 1` (local-log-polynomial).** Fit ``\log f`` by a degree-``d`` polynomial through a binned
  **local Poisson** likelihood — ``\log \mathbb{E}[n_j] = \log(nb) + \sum_{l \le d} \beta_l (c_j-x)^l``,
  weights ``K((c_j-x)/h_j)``, by Newton/IRLS — and report ``\hat f(x) = e^{\hat\beta_0}`` (Loader 1996).
  Degree 1 is exact for an exponential density (tracks exponential tails), degree 2 for a Gaussian. The
  two knobs compose. On a **bounded support** the degree-0 KDE accepts `bounds = (a, b)`, renormalizing
  by the kernel mass ``\int_a^b K_h(x-u)\,\mathrm{d}u`` inside the support (Jones 1993) to remove the
  ``O(1)`` boundary bias that otherwise halves the estimate at the edge.
- **Test.** `test/densityestimation.jl` (every grid corner recovers the ``N(0,1)`` pdf into the tail and
  integrates to ``\approx 1``; the variable bandwidth resolves the sharp Laplace cusp the fixed KDE
  oversmooths; degree-1 beats the KDE in a Laplace tail).

## Asymmetric kernels (§4.1) — `density`'s `gamma_kernel` / `beta_kernel`

- **Estimand.** ``f_X`` on a bounded or positive support.
- **Estimator.** the varying-shape estimate ``\hat f(x) = \tfrac1n \sum_j n_j K(c_j;\,x,\,b)`` with ``K(\cdot;x,b)`` the modified ``\mathrm{Gamma}(\rho(x),\,b)`` density (gamma kernel, support ``[0,\infty)``; Chen 2000) or the modified ``\mathrm{Beta}`` density (beta kernel, support ``[0,1]``; Chen 1999) — Chen's bias-reduced ``f_2`` shapes (interior shape ``x/b``, kernel mean ``\approx x``, removing the ``f'`` bias of the naive ``x/b+1`` form); the support matches the data, so no mass leaks past the boundary. A reduction over the bin counts (mergeable); default ``b = \hat\sigma\,n^{-2/5}``.
- **Consistency claim.** bias ``O(b)`` uniformly *including at the boundary* (where a symmetric kernel is ``O(h)`` and inconsistent); interior variance ``O(n^{-1}b^{-1/2})``; optimal ``b \sim n^{-2/5}``, MISE ``\sim n^{-4/5}``. Not self-normalizing (an ``O(b)`` mass deviation). The modified ``f_2`` shape lowers the MISE constant (Chen proves ``\mathrm{MISE}(f_2) \le \mathrm{MISE}(f_1)``) — measured ``\approx 40\%`` lower ISE than the naive form on a log-normal.
- **Test.** `test/densityestimation.jl` (gamma recovers the exponential at the boundary, beating the symmetric KDE there; beta recovers ``\mathrm{Uniform}(0,1)`` at both edges; integrates to ``\approx 1``; mergeability; degree/adaptive/bounds rejected).

## Parametric start & normalization (§4.1) — `density`'s `start` / `normalize`

- **Estimand.** ``f_X``, exactly as above; `start` and `normalize` change the *estimator*, not the target.
- **Estimator (`start`, Hjort & Glad 1995).** The semiparametric estimate ``\hat f(x) = p(x;\hat\theta)\,\hat r_h(x)``, where ``p(\cdot;\hat\theta)`` is a parametric *start* fit by **method of moments** from the recovered ``\hat\mu``, ``\hat\sigma^2``, and ``\hat r_h(x) = \tfrac1{nh}\sum_j (n_j/p(c_j;\hat\theta))\,K((x-c_j)/h)`` is a KDE of the *correction* ``r = f_X/p`` (the start's local relative error). Five starts are provided — `normal_start`, `lognormal_start`, `gamma_start`, `exponential_start`, `inverse_gaussian_start` — each a density whose parameters are fixed by the first two moments, so ``\hat\theta`` is itself a post-hoc functional of the summary and the whole estimate is **mergeable** (``\hat\theta`` merges through the moments; the reweighted counts add). Degree-0, fixed-``h``, symmetric kernel.
- **Consistency claim.** For a fixed start with bounded, twice-differentiable correction ``r``, consistent for ``f_X`` at ``\operatorname{bias} = O(h^2)``, ``\operatorname{var} = O(1/(nh))`` (MISE ``\sim n^{-4/5}``), but with the **leading bias governed by ``r'' = (f_X/p)''`` in place of ``f_X''``** and the variance ``R(K)f_X(x)/(nh)`` unchanged to leading order. When the start is correct ``r \equiv 1``, ``r'' \equiv 0``: the leading bias vanishes and the estimate collapses to the parametric fit. So the start is a **bias reduction** that pays in the bias-dominated regime (skewed / heavy-tailed shapes, oversmoothed ``h``), not in the large-``n`` limit where the KDE is already nearly unbiased.
- **Estimator (`normalize`).** ``\hat f / Z`` with ``Z = \int \hat f`` over the summary's range intersected with the estimator's support, by the trapezoidal rule on a fine grid — a deterministic, mergeable functional of the estimate making ``\int \hat f = 1`` **exactly** (the asymmetric kernels and a parametric start are not self-normalizing; the symmetric KDE already is, to the grid's ``O(\Delta^2)``).
- **Test.** `test/densityestimation.jl` (a normal start at a wide ``h`` recovers the Gaussian peak the plain KDE flattens by ``\approx 10\%``, ISE ``> 20\times`` smaller; a log-normal start tracks log-normal data ``> 3\times`` better than the symmetric KDE, and pointwise into the right tail; `normalize` drives ``\int\hat f`` to ``1`` for a start and for the gamma kernel; mergeability; the rejected knob combinations).

## Pointwise standard error (§4.1) — `density_standard_error`

- **Estimand.** the standard error of the density estimator, ``\mathrm{SE}(x) = \sqrt{\operatorname{Var}\hat f_h(x)}`` with the KDE asymptotic variance ``\operatorname{Var}\hat f_h(x) \approx R(K)\,f_X(x)/(nh)``.
- **Estimator.** ``\sqrt{R(K)\,\hat f_h(x)/(nh)}`` from the recovered ``\hat f_h``, ``n``, ``h``, and ``R(K) = \int K^2`` (by quadrature). The ``z``-band ``\hat f_h(x) \pm z\,\mathrm{SE}(x)`` is a variability band about ``\mathbb{E}\hat f_h``, covering ``f_X`` where the ``O(h^2)`` smoothing bias is small. Degree 0/1 (the variable bandwidth and degree ``\ge 2`` are a coverage gap).
- **Test.** `test/densityestimation.jl` (the SE equals ``\sqrt{R(K)\hat f/(nh)}`` exactly; the 95\% band covers the ``N(0,1)`` truth over ``\ge 75\%`` of the interior; the SE shrinks with ``n``; vector / empty / error forms).

## Bandwidth selection (§4.5) — `bandwidth`

- **Estimand.** the AMISE-optimal bandwidth ``h^\star = \arg\min_h \int \mathbb{E}[(\hat f_h(x) - f_X(x))^2]\,\mathrm{d}x`` for a kernel ``K`` and local-polynomial degree ``d`` — the asymptotically optimal smoothing scale, itself ``O(n^{-1/(2p+1)})`` (a finite-sample selector, not a fixed population constant).
- **Estimator** (`rule = :normal`, the kernel- and degree-aware normal-reference). ``h = C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}`` with bias order ``p = 2`` (degree 0/1) or ``p = 4`` (degree 2/3), and constant ``C(K,d) = [W_{ikk}(K,d)\,2^{p+1}\sqrt{\pi}/(2p-1)!!]^{1/(2p+1)}``. The canonical factor is ``W_{ikk}(K,d) = R(K)/\mu_2(K)^2`` at degree 0/1, and the fourth-order equivalent-kernel analog ``W_{ikk} = (4!^2/8)\,R(K^\star)/\mu_4(K^\star)^2`` with ``K^\star = K\,(n_4 - v^2 n_2)/(n_0 n_4 - n_2^2)`` at degree 2/3 (Loader 1996); here ``R(K) = \int K^2``, ``\mu_2(K) = \int u^2 K`` (by quadrature) and ``n_k = \int v^k K``. ``\hat\sigma`` and ``n`` are read from the one summary.
- **Estimator** (`rule = :sj`, the data-driven Sheather–Jones plug-in; degree 0/1). The curvature ``\psi_4 = R(f'')`` and its ``\psi_6`` pilot are estimated from the bin counts, ``\hat\psi_r(g) = (n^2 g^{r+1})^{-1}\sum_{j,k} n_j n_k\,\varphi^{(r)}((c_j - c_k)/g)`` (Gaussian ``\varphi``, the ``n`` self-pairs removed), and the solve-the-equation bandwidth (Sheather & Jones 1991) is the root of ``h = [R(\varphi)/(\mu_2(\varphi)^2\,\hat\psi_4(\alpha\,h^{5/7})\,n)]^{1/5}``; the eval kernel then multiplies by the canonical ``\delta(K)/\delta(\varphi)``. Both functionals are reductions of the bin counts, hence mergeable; degree ``\ge 2`` is a recorded coverage gap. The **improved** variant (`rule = :isj`; Botev, Grotowski & Kroese 2010) instead reads ``\|f^{(s)}\|^2`` off the DCT of the bin counts and solves the AMISE fixed point ``t = \xi\gamma^5(t)`` — also a mergeable functional of the counts, and robust on strongly multimodal data.
- **Estimator** (`rule = :aic`, the fit-based information criterion; degree 0/1). Minimise ``\mathrm{AIC}(h) = -2\sum_j n_j\log\hat f_h(c_j) + 2\,\mathrm{df}(h)`` of the binned KDE over a log-spaced ``h``-grid, with ``\mathrm{df}(h) = (K(0)/h)\sum_j n_j/(n\hat f_h(c_j))``; both terms diverge as ``h\to0``, giving a stable interior optimum. Mergeable (a functional of the counts); the eval kernel multiplies by ``\delta(K)/\delta(\varphi)``. The leave-one-out cross-validation selectors (likelihood-CV, least-squares LSCV) are a recorded coverage gap — binning ties drive their objective to ``h\to0``.
- **Consistency claim.** with ``\hat\sigma \to \sigma`` and the Gaussian reference ``R(f^{(p)}) = (2p-1)!!/(2^{p+1}\sqrt{\pi}\,\sigma^{2p+1})``, ``C(K,d)`` is the exact AMISE-optimal constant: degree 0/1 gives ``h \sim n^{-1/5}`` (MISE ``\sim n^{-4/5}``), degree 2 the higher-order ``h \sim n^{-1/9}`` (MISE ``\sim n^{-8/9}``); the Gaussian/degree-0 case is Silverman ``(4/3)^{1/5}``, the Epanechnikov ``\approx 2.34``.
- **Test.** `test/bandwidth.jl` (Gaussian degree 0 reproduces Silverman, degree 1 equals degree 0, the Epanechnikov ``h`` is ``2.214\times`` the Gaussian, degree 2 matches the closed-form ``n^{-1/9}`` rule and exceeds degree 0; coverage-gap and sample-size errors; plus the Sheather–Jones tests — SJ ``\approx`` the normal reference on Gaussian data, much tighter on a bimodal mixture, the exact ``\delta(K)`` kernel ratio, mergeability, the degree ``\ge 2`` gap; and the AIC tests — a stable optimum below the normal reference on Gaussian data, far smaller on a bimodal mixture, the exact ``\delta(K)`` ratio, mergeability, and the deliberately-absent CV rules erroring) and `test/oracles.jl` (the kernel-aware ratio on the ``N(0,1)`` oracle).

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
  ``\lesssim 2\times10^{-2}``; `conditional_variance` recovers ``1-\rho^2`` to a few percent (the
  simple ``\hat{\mathbb E}[Y^2\mid x]-\hat m^2`` estimator carries the local-linear smoothing bias of
  the curved ``E[Y^2\mid x]``; the Fan–Yao residual refinement is forthcoming). The ``n^{-2/5}`` rate
  test is forthcoming.

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
  carries the local-linear bias from the curvature of ``\mathbb{E}[Y^2\mid x]`` (the Fan–Yao residual
  refinement removes it). The finite-``\tau`` coefficients approach ``D^{(n)}(x)`` as ``\tau\to0`` up to
  the ``O(\tau)`` finite-time bias (the Bell-polynomial correction, forthcoming).
- **Test.** `test/oracles.jl` — the Ornstein–Uhlenbeck oracle (exact stationary lagged pairs from
  ``\mathrm dX = -\theta X\,\mathrm dt + \sigma\,\mathrm dW``): `drift` recovers
  ``x(e^{-\theta\tau}-1)/\tau`` to ``\le 2\times10^{-2}`` and `diffusion` recovers the closed-form
  ``[s^2(1-e^{-2\theta\tau}) + x^2(e^{-\theta\tau}-1)^2]/(2\tau)`` (``s^2 = \sigma^2/2\theta``) to
  ``8\%``. `test/markovinference.jl` — the variance decomposition
  ``M_2 = \operatorname{Var}(Y\mid x) + M_1^2``, the ``1/(n!\,\tau)`` scaling, the ``n\ge3`` and
  static-lag coverage-gap errors, and mergeability. The ``n^{-2/5}`` rate test is forthcoming.

## Conditional & transition density (§4.3) — `conditional_density`, `transition_density`

- **Estimand.** ``f_{Y\mid X}(y\mid x) = f_{X,Y}(x,y)/f_X(x)``, and the two-time transition density
  ``p_\tau(y\mid x) = f_{X_{t+\tau}\mid X_t}(y\mid x)`` (the same object on a lag-``\tau`` joint).
- **Estimator.** the binned double-kernel (Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996): from
  the cell counts ``n_{jk}`` at centers ``(c^x_j, c^y_k)`` with ``N_j = \sum_k n_{jk}``,

  ```math
  \hat f(y\mid x) =
    \frac{\sum_{j,k} n_{jk}\, K_{h_1}(c^x_j - x)\, K_{h_2}(c^y_k - y)}{\sum_j N_j\, K_{h_1}(c^x_j - x)} ,
  ```

  computed as ``\hat f(y\mid x) = \big(\sum_k P_k\, K_{h_2}(c^y_k - y)\big) / (\sqrt{2\pi}\,h_2 D)`` with
  the X-kernel-weighted profile ``P_k = \sum_j n_{jk} K_{h_1}(c^x_j - x)`` and ``D = \sum_k P_k``. The
  X-kernel normalization cancels in the ratio; the surviving ``1/(\sqrt{2\pi}\,h_2)`` makes ``\int \hat f(y\mid x)\,\mathrm dy = 1`` exactly.
- **Consistency claim** (``h_1, h_2 \to 0``, ``n h_1 h_2 \to \infty``): the pointwise bias is ``\tfrac12 h_1^2 \mu_2(K)\big[\partial_{xx} f_{Y\mid X} + 2\,\partial_x f_{Y\mid X}\,\partial_x f_X / f_X\big] + \tfrac12 h_2^2 \mu_2(K)\,\partial_{yy} f_{Y\mid X}`` and the variance ``R(K)^2 f_{Y\mid X}/(n h_1 h_2 f_X)``, giving the slower conditional rate ``O(n^{-2/3})`` at the optimal ``h \sim n^{-1/6}`` (the bias is worst at the conditional mode).
- **Test.** `test/conditionaldensity.jl` — the bivariate-Gaussian fixture ``Y\mid X=x \sim N(\rho x, 1-\rho^2)``: ``\int \hat f\,\mathrm dy \approx 1`` and pointwise recovery to ``3\times10^{-2}``; the
  scalar/vector forms, the `transition_density` equivalence, off-support `NaN`, bandwidth checks, and
  mergeability. `test/oracles.jl` — the Ornstein–Uhlenbeck oracle: `transition_density` recovers the
  closed-form ``p_\tau(y\mid x) = N(x e^{-\theta\tau},\, s^2(1-e^{-2\theta\tau}))`` to ``4\times10^{-2}``
  (the ``O(h^2)`` mode bias). The ``n^{-2/3}`` rate test is forthcoming.
