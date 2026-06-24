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

The implemented entries are test-linked; the conditional entry is the format for what lands next.

## Per-bin power-sums about the bin center (`Histograms`, the summaries)

- **Object.** per bin, the count and the **power-sums about the fixed bin center** ``S_p =
  \sum_{i\in\text{bin}}(X_i - c_j)^p`` for ``p = 1 \dots`` `order` (and order-2 bivariate ``S_{rs}``
  per cell for `Histogram2D`). Intrinsically finite-sample — no population estimand.
- **Accumulation.** an atomic scatter-add (`observe!` / `observe_pairs!`): one work-item per sample
  adds ``(X_i - c_j)^p`` to its bin's ``S_p``, in the accumulator type (`Float64` by default). `merge!` is elementwise addition (same
  binning ``\Rightarrow`` same centers). The local offsets are well-conditioned, so no compensated
  summation is used.
- **Recovery.** the global central moments are the **shift formula**

  ```math
  \sum_i (X_i - \bar X)^p = \sum_j \sum_{k=0}^{p} \binom{p}{k}\,(c_j - \bar X)^{p-k}\, S_k^{(j)},
  \qquad S_0^{(j)} = n_j, \quad \bar X = \tfrac1n \sum_j \bigl(c_j n_j + S_1^{(j)}\bigr),
  ```

  which is **exact** (the ``S_p`` use the actual ``X_i``, not the center) and numerically stable for
  the even ``p`` the consumers need.
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
  integrates to 1; the Silverman/Scott bandwidths match their formulas).

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

## Conditional increment moments, drift & diffusion (§4.4) — `increment_moment`, `drift`, `diffusion`

- **Estimand.** the conditional increment moment and the Kramers–Moyal coefficient

  ```math
  M_n(x, \tau) = \mathbb{E}\!\left[(X_{t+\tau} - X_t)^n \mid X_t = x\right], \qquad
  D^{(n)}(x) = \frac{1}{n!}\lim_{\tau\to0} \frac{M_n(x,\tau)}{\tau},
  ```

  with `drift` ``\equiv D^{(1)}`` and `diffusion` ``\equiv D^{(2)}``. Model-neutral — a conditional
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
