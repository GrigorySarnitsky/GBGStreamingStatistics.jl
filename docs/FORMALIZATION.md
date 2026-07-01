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
number of samples, ``f_X`` the marginal density of ``X`` (the conditioning variable in the joint estimators).

The entries are grouped by the lens of `ARCHITECTURE.md` §2: the **atoms** (the power-sums summary, the
cone, the kernel-weighting atom & jet) first, then the two readout families — the **kernel-weighting
readouts** (§4.A) and the **power-sum quadrature readouts** (§4.B). A forthcoming operation gets its
entry — in the same shape (estimand / estimator / consistency / test) — when it lands, its estimand
recorded in `ARCHITECTURE.md` §4 meanwhile.

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
  the even ``p`` the consumers need; the recovered moments are exposed publicly (in `Expectations`, §4.B,
  the ``g = x^p`` quadrature case) as `mean`, `var`, `std`, and `central_moment` — `var` and `std` apply
  Bessel's ``n-1`` correction by default (matching `Statistics`), `corrected = false` the ``1/n`` form.
- **Quadrature readout (the power-sum atom; §4.B).** the same power-sums give any smooth functional by
  Taylor contraction about the bin centers, ``\mathbb{E}[g(X)] = \tfrac1n \sum_j \sum_{p=0}^{\texttt{order}}
  g^{(p)}(c_j)/p!\,S_p^{(j)} + O(\Delta^{\texttt{order}+1})`` — **exact** for polynomials of degree
  ``\le \texttt{order}`` (the moment recovery above is the ``g = x^p`` case), ``O(\Delta^{\texttt{order}+1})``
  for general smooth ``g``. The general ``\mathbb{E}[g(X)]``, its standard error, and the pushforward
  density ``f_{g(X)}`` are **built** in `Expectations` (the §4.B section below), all through this same
  contraction — no new accumulation.
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

## The kernel-weighting atom & local-polynomial jet (`Kernels._kernel_weight`, `LocalPolynomial.local_fit`)

- **Kernel-weighting atom (`_kernel_weight`).** the single per-bin weight ``w_j`` at a query point ``x``
  and bandwidth ``h``, in one of two protocols: **symmetric** ``w_j = K((c_j - x)/h)/h`` (the ``1/h``
  folded in, normalized), or **asymmetric** ``w_j = K(c_j;\,x,\,h)`` (a varying-shape density in the data
  argument, its own normalization). Every smoothing readout consumes ``w_j`` — a count sum (`density`), a
  ratio (the jet below), a ratio of ratios (`conditional_density`) — so the two protocols live in exactly
  one place. Folding the ``1/h`` into the weight is ``\beta``-invariant for the local-polynomial solve
  (both the design matrix and the right-hand side scale by ``1/h``), so it leaves every symmetric estimate
  unchanged.
- **Jet (`local_fit`).** the kernel-weighted local-polynomial **jet** at ``x``: the degree-``d``
  weighted-least-squares fit of a polynomial in ``u = \ell_j - x`` to a binned 1-D measure (locations
  ``\ell_j``, design weights ``w_j``, response totals ``r_j``), returning ``\beta = (\beta_0, \dots,
  \beta_d)`` — ``\beta_0`` the value, ``\beta_1`` the slope, ``2\beta_2`` the curvature. Form the design and
  response moments ``T_\ell = \sum_j w_j\,(\ell_j - x)^\ell`` (``\ell = 0 \dots 2d``) and
  ``P_\ell = \sum_j w_j\,r_j\,(\ell_j - x)^\ell`` (``\ell = 0 \dots d``), the ``w_j`` from `_kernel_weight`,
  and solve the ``(d+1)\times(d+1)`` Hankel normal equations ``[\,T_{i+\ell}\,]\,\beta = (P_0, \dots, P_d)``
  (a small solve, never an explicit inverse); a degenerate (singular) design returns an all-NaN jet.
  `conditional_mean` is ``\beta_0`` of one fit; the local-likelihood `density` reuses the solve as its IRLS
  inner step. **An asymmetric kernel has no jet** (it is not a translation kernel), so `local_fit` admits
  it at `degree = 0` only — the one guard, a clean error above it.
- **Test.** `test/localpolynomial.jl` — on a binned measure whose per-bin response mean is exactly a
  polynomial ``P``, the jet recovers ``(P(x),\,P'(x),\,P''(x)/2)`` to ``10^{-6}`` (degree 2) and
  ``(P,\,P')`` (degree 1); degree 0 is the Nadaraya–Watson ratio ``P_0/T_0``; a degenerate design gives an
  all-NaN jet; the kernel normalization cancels. The ``\beta``-invariance of the ``1/h`` fold and the
  asymmetric-degree-0 admission ride on the unchanged `density` and `conditional_mean` oracle / rate tests.
  The kernel constants the weighting depends on (``\int K = 1``, ``\mu_2(K)``, ``R(K)``) are verified by
  independent quadrature in `test/kernels.jl`.

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
  ``n^{-2p/(2p+1)}`` — slope ``\approx -0.80`` for degree 0/1 (``p=2``), the steeper ``\approx -0.89``
  (``n^{-8/9}``, ``p=4``) for degree 2, the Abramson variable bandwidth sharing the degree-0 order; plus the
  **binning-error convergence** sub-test — the binned estimate approaches the exact KDE at ``O(\Delta^2)``
  for the default linear/centroid `binning` (slope ``\approx 2`` in ``\Delta``) and ``O(\Delta)`` for
  `binning = :simple` (slope ``\approx 1``).

## Local-log-polynomial & adaptive density (§4.1) — `density`'s `degree` / `adaptive`

- **Estimand.** ``f_X``, exactly as for the degree-0 `density` above; the options change the *estimator*,
  not the target.
- **`adaptive` (Abramson variable bandwidth).** A fixed-bandwidth pilot ``\tilde f`` at ``h`` sets
  per-bin bandwidths ``h_j = h\,(\tilde f(c_j)/g)^{-1/2}`` (``g`` the geometric mean of ``\tilde f`` over
  the sample), widening the kernel in the sparse tails: ``\operatorname{bias} = O(h^4)`` at interior
  smooth ``x`` (Abramson 1982), with far less tail variance than a fixed ``h``.
- **`degree` ``\ge 1`` (local-log-polynomial).** Fit ``\log f`` by a degree-``d`` polynomial through a binned
  **local Poisson** likelihood — ``\log \mathbb{E}[n_j] = \log(nb) + \sum_{l \le d} \beta_l (c_j-x)^l``,
  weights ``K((c_j-x)/h_j)``, by Newton/IRLS — and report ``\hat f(x) = e^{\hat\beta_0}`` (Loader 1996).
  Degree 1 is exact for an exponential density (tracks exponential tails), degree 2 for a Gaussian. The
  two options compose. On a **bounded support** the degree-0 KDE accepts `bounds = (a, b)`, renormalizing
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
- **Test.** `test/marginaldensity.jl` (a normal start at a wide ``h`` recovers the Gaussian peak the plain KDE flattens by ``\approx 10\%``, ISE ``> 20\times`` smaller; a log-normal start tracks log-normal data ``> 3\times`` better than the symmetric KDE, and pointwise into the right tail; `normalize` drives ``\int\hat f`` to ``1`` for a start and for the gamma kernel; mergeability; the rejected option combinations).

## Pointwise standard error (§4.1) — `density_standard_error`, `density_with_standard_error`

- **Estimand.** the standard error of the density estimator, ``\mathrm{SE}(x) = \sqrt{\operatorname{Var}\hat f_h(x)}`` with the finite-sample KDE variance ``\operatorname{Var}\hat f_h(x) \approx R(K)\,f_X(x)/(nh) - f_X(x)^2/n`` (Wand & Jones 1995).
- **Estimator.** ``\sqrt{R(K)\,\hat f_h(x)/(nh) - \hat f_h(x)^2/n}`` from the recovered ``\hat f_h``, ``n``, ``h``, and ``R(K) = \int K^2`` (by quadrature) — the ``z``-band ``\hat f_h(x) \pm z\,\mathrm{SE}(x)`` about ``\mathbb{E}\hat f_h``, covering ``f_X`` where the ``O(h^2)`` smoothing bias is small. **It tracks the estimator's options:** `adaptive` uses the local Abramson bandwidth ``h(x)`` in place of ``h``; `bounds` (Jones renormalization) uses the boundary-restricted roughness ``R_a(x) = \int_{\text{support}} K_h(x-u)^2\,\mathrm du`` and divides by the kernel mass ``a(x)`` inside the support, inflating the variance near an edge; a parametric `start` leaves the leading variance unchanged (its ``\hat\theta`` uncertainty is an ``O(1/n)`` omission). `density_with_standard_error` returns the estimate and this band from **one** option specification (the `NamedTuple` ``(\texttt{value}, \texttt{standard\_error})``), so they describe the same estimator by construction. Degree 0/1 and a symmetric kernel (the asymmetric-kernel variance and degree ``\ge 2`` are coverage gaps).
- **Test.** `test/marginaldensity.jl` — the SE equals the two-term ``\sqrt{R(K)\hat f/(nh) - \hat f^2/n}`` exactly; a **Monte-Carlo** check (the predicted SE matches the empirical sampling standard deviation of ``\hat f(0)`` across independent datasets to ``20\%``); `adaptive` lowers the tail variance and `bounds` inflates the edge variance; `density_with_standard_error` matches the separate calls exactly; the 95\% band covers the ``N(0,1)`` truth over ``\ge 75\%`` of the interior; vector / empty / error forms.

## Bandwidth selection (§4.5) — `bandwidth`

- **Estimand.** the AMISE-optimal bandwidth ``h^\star = \arg\min_h \int \mathbb{E}[(\hat f_h(x) - f_X(x))^2]\,\mathrm{d}x`` for a kernel ``K`` and local-polynomial degree ``d`` — the asymptotically optimal smoothing scale, itself ``O(n^{-1/(2p+1)})`` (a finite-sample selector, not a fixed population constant).
- **Estimator** (`rule = :normal`, the kernel- and degree-aware normal-reference). ``h = C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}`` with bias order ``p = 2`` (degree 0/1) or ``p = 4`` (degree 2/3), and constant ``C(K,d) = [W_{ikk}(K,d)\,2^{p+1}\sqrt{\pi}/(2p-1)!!]^{1/(2p+1)}``. The canonical factor is ``W_{ikk}(K,d) = R(K)/\mu_2(K)^2`` at degree 0/1, and the fourth-order equivalent-kernel analog ``W_{ikk} = (4!^2/8)\,R(K^\star)/\mu_4(K^\star)^2`` with ``K^\star = K\,(n_4 - v^2 n_2)/(n_0 n_4 - n_2^2)`` at degree 2/3 (Loader 1996); here ``R(K) = \int K^2``, ``\mu_2(K) = \int u^2 K`` (by quadrature) and ``n_k = \int v^k K``. ``\hat\sigma`` and ``n`` are read from the one summary.
- **Estimator** (`rule = :sj`, the data-driven Sheather–Jones plug-in; degree 0/1). The curvature ``\psi_4 = R(f'')`` and its ``\psi_6`` pilot are estimated from the bin counts, ``\hat\psi_r(g) = (n^2 g^{r+1})^{-1}\sum_{j,k} n_j n_k\,\varphi^{(r)}((c_j - c_k)/g)`` (Gaussian ``\varphi``, the ``n`` self-pairs removed), and the solve-the-equation bandwidth (Sheather & Jones 1991) is the root of ``h = [R(\varphi)/(\mu_2(\varphi)^2\,\hat\psi_4(\alpha\,h^{5/7})\,n)]^{1/5}``; the eval kernel then multiplies by the canonical-bandwidth ratio ``\delta(K)/\delta(\varphi)`` with ``\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}``. Both functionals are reductions of the bin counts, hence mergeable; degree ``\ge 2`` is a recorded coverage gap. The **improved** variant (`rule = :isj`; Botev, Grotowski & Kroese 2010) instead reads ``\|f^{(s)}\|^2`` off the DCT of the bin counts and solves the AMISE fixed point ``t = \xi\gamma^5(t)`` — also a mergeable functional of the counts, and robust on strongly multimodal data.
- **Estimator** (`rule = :aic`, the fit-based information criterion; degree 0/1). Minimise ``\mathrm{AIC}(h) = -2\sum_j n_j\log\hat f_h(c_j) + 2\,\mathrm{df}(h)`` of the binned **Gaussian** KDE over a log-spaced ``h``-grid, with ``\mathrm{df}(h) = (\varphi(0)/h)\sum_j n_j/(n\hat f_h(c_j))`` (``\varphi`` the Gaussian); both terms diverge as ``h\to0``, giving a stable interior optimum. Mergeable (a functional of the counts); the eval kernel multiplies by ``\delta(K)/\delta(\varphi)``. The leave-one-out cross-validation selectors (likelihood-CV, least-squares LSCV) are a recorded coverage gap — binning ties drive their objective to ``h\to0``.
- **Estimator** (`rule = :plugin`, the **asymmetric-kernel** Chen plug-in; `gamma_kernel` / `beta_kernel`). The gamma / beta estimators take a scale ``b`` (rate ``n^{-2/5}``), not ``h``; their AMISE (Chen 1999, 2000; modified ``f_2``) is ``A b^2 + C b^{-1/2}/n`` with ``A = \tfrac14\int w^2 (f'')^2``, ``C = (2\sqrt\pi)^{-1}\mathbb{E}[w(X)^{-1/2}]`` and local-scale weight ``w(x) = x`` (gamma) or ``x(1-x)`` (beta), giving ``b^\star = [\mathbb{E}[w(X)^{-1/2}]/(2\sqrt\pi\,n\int w^2 (f'')^2)]^{2/5}``. ``\mathbb{E}[w(X)^{-1/2}]`` is an exact count moment; the curvature ``\int w^2 (f'')^2`` is a grid quadrature of the pilot KDE's ``\hat f''_g`` (pilot ``g = 1.3\,s\,n^{-1/9}`` at the robust scale ``s = \min(\hat\sigma,\,\mathrm{IQR}/1.349)``, the variance inflation ``R(\varphi'')\,\mathbb{E}[w^2]/(n g^5)`` removed), which — unlike the integration-by-parts ``\psi_4`` U-statistic — stays valid where ``f(0) > 0`` (the ``w \to 0`` boundary weight suppresses the pilot's edge leakage). A deterministic functional of the counts (mergeable). The symmetric `:sj` / `:isj` / `:aic` do not apply; a fully data-adaptive curvature pilot is a recorded refinement.
- **Consistency claim.** with ``\hat\sigma \to \sigma`` and the Gaussian reference ``R(f^{(p)}) = (2p-1)!!/(2^{p+1}\sqrt{\pi}\,\sigma^{2p+1})``, ``C(K,d)`` is the exact AMISE-optimal constant: degree 0/1 gives ``h \sim n^{-1/5}`` (MISE ``\sim n^{-4/5}``), degree 2 the higher-order ``h \sim n^{-1/9}`` (MISE ``\sim n^{-8/9}``); the Gaussian/degree-0 case is Silverman ``(4/3)^{1/5}``, the Epanechnikov ``\approx 2.34``.
- **Test.** `test/bandwidth.jl` (Gaussian degree 0 reproduces Silverman, degree 1 equals degree 0, the Epanechnikov ``h`` is ``2.214\times`` the Gaussian, degree 2 matches the closed-form ``n^{-1/9}`` rule and exceeds degree 0; coverage-gap and sample-size errors; plus the Sheather–Jones tests — SJ ``\approx`` the normal reference on Gaussian data, much tighter on a bimodal mixture, the exact ``\delta(K)`` kernel ratio, mergeability, the degree ``\ge 2`` gap; and the AIC tests — a stable optimum below the normal reference on Gaussian data, far smaller on a bimodal mixture, the exact ``\delta(K)`` ratio, mergeability, and the deliberately-absent CV rules erroring; and the asymmetric Chen plug-in — `:normal` reproduces the ``\hat\sigma\,n^{-2/5}`` rule of thumb, `:plugin` lowers the ISE (``\ge 1.25\times`` on log-normal, ``\ge 1.11\times`` on gamma; ``\approx 1.4\text{–}1.7\times`` typically, seed-audited) and stays ``\approx`` neutral on the near-exponential boundary case, gamma and beta kernels, mergeability, and the inapplicable `:sj`/`:isj`/`:aic` and `degree` ``\ge 1`` rules erroring; the Improved Sheather–Jones tests (`:isj` ``\approx`` the normal reference on Gaussian, ``< 0.5\times`` it on a bimodal mixture, the exact ``\delta(K)`` ratio, mergeability); and the **Marron–Wand claw** stress, where the data-driven `:sj` / `:isj` resolve the five spikes the normal reference oversmooths — ``\approx`` an order of magnitude lower ISE) and `test/oracles.jl` (the kernel-aware ratio on the ``N(0,1)`` oracle).

## Local-linear conditional mean (§4.2) — `conditional_mean`

- **Estimand.** ``m(x) = \mathbb{E}[Y \mid X = x] = \int y\, f_{Y\mid X}(y \mid x)\,\mathrm{d}y``,
  the ``L^2(\sigma(X))`` projection of the response onto the conditioning variable.
- **Estimator.** the local-linear fit ``\hat m_h(x) = (T_2 P_0 - T_1 P_1)/(T_0 T_2 - T_1^2)`` from
  the kernel-weighted power-sums ``T_\ell = \sum_i K_h(X_i-x)(X_i-x)^\ell`` and
  ``P_\ell = \sum_i K_h(X_i-x)(X_i-x)^\ell Y_i`` — recovered from the joint `Histogram2D`'s per-X-column
  counts and response totals ``\sum Y`` (exact from the per-cell power-sums); the ``g = \mathrm{identity}``
  case of `conditional_expectation`.
- **Consistency claim** (interior ``x``, ``h \to 0``, ``nh \to \infty``):

  ```math
  \operatorname{bias}\hat m_h(x) \approx \tfrac12 h^2 \mu_2(K)\, m''(x),
  \qquad
  \operatorname{Var}\hat m_h(x) \approx \frac{R(K)\,\sigma^2(x)}{n\,h\,f_X(x)},
  ```

  giving the optimal rate ``\hat m_{h^\star}(x) - m(x) = O_p(n^{-2/5})`` at ``h^\star \sim n^{-1/5}``.
- **Test.** `test/conditionalexpectations.jl` — the bivariate-Gaussian oracle ``Y = \rho X + \sqrt{1-\rho^2}\,
  \varepsilon`` (so ``\mathbb{E}[Y\mid X=x] = \rho x``): `conditional_mean` recovers the slope to
  ``\lesssim 2\times10^{-2}``, with the `degree` (0/1/2) and `kernel` options exercised (degree-2 and
  Epanechnikov both recover the linear mean). `test/conditional_convergence.jl` — the ``n^{-2/5}`` rate:
  on a curved-mean oracle (``m(x) = \sin x``, so the local-linear fit carries a genuine ``O(h^2)`` bias)
  at ``h \propto n^{-1/5}`` the integrated squared error falls at ``n^{-4/5}`` (the ``n^{-2/5}`` RMSE), the
  fitted log–log slope ``\approx -0.80``.

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
- **Test.** `test/conditionalexpectations.jl` — recovers ``1-\rho^2`` on the (linear-mean) bivariate-Gaussian oracle to
  a few percent, and on a **heteroscedastic curved-mean** oracle (``X\sim N(0,1)``, ``m(x)=x^2``,
  ``\sigma^2(x)=0.2+0.1x^2``) the residual estimator recovers ``\sigma^2(x)`` (degree 2) and beats the
  plug-in in the curved tail (the plug-in inflated ``>1.5\times`` by the ``m''`` curvature); the universal
  contract via `check_point_estimator` (scalar ``=`` vector, functional of the merged summary,
  ``\mathrm{NaN}`` off-support). Its rate is the same ``n^{-2/5}`` local-regression rate (a second
  local-polynomial fit, of the squared residuals); it shares the `conditional_mean` slope test
  (`test/conditional_convergence.jl`) rather than carrying a separate one.

## Conditional expectation of a functional (§4.2) — `conditional_expectation`

- **Estimand.** ``\mathbb{E}[g(Y) \mid X = x] = \int g(y)\, f_{Y\mid X}(y \mid x)\,\mathrm{d}y``, the
  conditional expectation of a smooth functional of the response (``g = \mathrm{identity}`` the regression
  function ``m(x)``; ``g(y) = y^2`` gives ``\mathbb{E}[Y^2\mid x]``).
- **Estimator.** a local-polynomial regression **in ``x``** (as for `conditional_mean`) of a per-``X``-column
  **``y``-quadrature** response ``\sum_i g(Y_i)`` formed from the stored order-2 ``Y``-power-sums by the §4.B
  quadrature (``a_0 n + a_1 S_{01} + a_2 S_{02}`` with ``a_p`` the degree-2 interpolant coefficients of ``g``
  on the ``Y``-bin); the estimate is ``\beta_0``. It **composes the two atoms** — power-sum quadrature in
  ``y`` then kernel-weighting in ``x`` — and reduces to `conditional_mean` at ``g = \mathrm{identity}``
  (folded in exactly).
- **Consistency claim.** the ``x``-side carries `conditional_mean`'s local-regression bias/variance (rate
  ``n^{-2/5}``); the ``y``-quadrature is **exact for ``g`` of degree ``\le 2`` in ``y``** and
  ``O(\Delta_y^2)`` for a general smooth ``g`` (the fixed order-2 ``Y``-power-sums of `Histogram2D`).
- **Test.** `test/conditionalexpectations.jl` — the bivariate-Gaussian oracle ``\mathbb{E}[\cos Y\mid x] =
  \cos(\rho x)\,e^{-(1-\rho^2)/2}`` (and ``\sin`` likewise) to ``\lesssim 2.5\times10^{-2}``; the exact fold
  ``\texttt{conditional\_mean}(x) = \texttt{conditional\_expectation}(\mathrm{identity}, x)``; the `:plugin`
  variance equal to ``\mathbb{E}[Y^2\mid x] - \mathbb{E}[Y\mid x]^2`` via the general reader; scalar ``=``
  vector; the universal contract via `check_point_estimator`.

## Conditional pushforward density (§4.2) — `conditional_pushforward_density`

- **Estimand.** ``f_{g(Y)\mid X}(z \mid x) = \sum_{y:\,g(y)=z} f_{Y\mid X}(y\mid x)/|g'(y)|``, the density of
  the transformed response ``Z = g(Y)`` given ``X = x`` (change of variables on the conditional law).
- **Estimator.** the conditional KDE of ``g(Y)`` written as a conditional expectation,
  ``\hat f_{g(Y)\mid X}(z\mid x) = \widehat{\mathbb{E}}[K_{h_z}(g(Y)-z)\mid x]`` — one `conditional_expectation`
  of the composed test function ``y\mapsto K_{h_z}(g(y)-z)``, no ``g^{-1}`` and no explicit Jacobian; a
  non-monotone ``g`` contributes from every preimage, ``g'=0`` is an integrable caustic. Two smoothing
  directions: the ``x``-conditioning regression (`hx`, `kernelx`, `degree`) and the ``z`` smoothing (`hz`,
  `kernelz`).
- **Consistency claim.** the standard KDE bias–variance in ``z`` (over the ``y``-quadrature), with a
  symmetric ``K_z`` giving ``\int\hat f\,\mathrm dz = 1`` exactly; needs ``g`` smooth on the ``Y``-bin scale
  and ``h_z \gtrsim |g'|\,\Delta_y``. At ``g = \mathrm{identity}``, `degree = 0` it coincides with the
  double-kernel `conditional_density` (up to the within-``Y``-cell ``O(\Delta_y^2)`` correction).
- **Test.** `test/conditionalexpectations.jl` — ``Z = e^Y\mid x \sim \mathrm{Lognormal}(\rho x, 1-\rho^2)``
  recovered pointwise (to ``8\%``) with ``\int\hat f\,\mathrm dz \approx 1``; scalar ``=`` vector; the
  positive-``h_z`` guard; the degree-0 identity-``g`` coincidence with `conditional_density`.

## Conditional expectation standard error (§4.2) — `conditional_expectation_standard_error`

- **Estimand.** the sampling standard error of ``\widehat{\mathbb{E}}[g(Y)\mid x]``,
  ``\mathrm{SE} = \sqrt{\ell(x)\,\operatorname{Var}(g(Y)\mid x)}`` with the local-regression **leverage**
  ``\ell(x) = \sum_j N_j w_j(x)^2`` (``w_j`` the effective weights of the fit's ``\beta_0``) and
  ``\operatorname{Var}(g(Y)\mid x) = \mathbb{E}[g^2\mid x] - \mathbb{E}[g\mid x]^2``.
- **Estimator.** ``\sqrt{\hat\ell(x)\,(\widehat{\mathbb{E}}[g^2\mid x] - \widehat{\mathbb{E}}[g\mid x]^2)}``,
  the **leverage × noise-scale** readout: ``\hat\ell(x) = u^\top Q\, u`` with ``u = M^{-1} e_0`` (``M`` the
  design Hankel ``[T_{i+j}]``) and ``Q`` the Hankel of the **squared-kernel** moments
  ``Q_\ell = \sum_j N_j K_h(c_j-x)^2 (c_j-x)^\ell`` (the `LocalPolynomial._leverage` primitive — the reference
  mechanism for the kernel-family SE), times the plug-in conditional variance of ``g(Y)`` (two
  `conditional_expectation` readouts). It is the exact analog of the marginal `expectation_standard_error`
  ``= \sqrt{\operatorname{Var}(g(X))/n}``, with the global-mean leverage ``1/n`` replaced by
  ``\sum_j N_j w_j^2``; at ``g = \mathrm{identity}`` it is the conditional mean's standard error. The paired
  `conditional_expectation_with_standard_error` returns this error with the estimate
  ``\widehat{\mathbb{E}}[g(Y)\mid x]`` from one shared fit (the ``(\texttt{value}, \texttt{standard\_error})``
  `NamedTuple`).
- **Consistency claim.** it estimates the leading-order sampling variance of the local fit (not the
  ``O(h^2)`` smoothing bias); at ``g = \mathrm{identity}`` the leading term is
  ``R(K)\sigma^2(x)/(n h f_X(x))``.
- **Test.** `test/conditionalexpectations.jl` — a **Monte-Carlo** check: the predicted SE (from one dataset)
  matches the empirical standard deviation of `conditional_mean` across ``K`` independent datasets (to
  ``20\%`` at the density peak, where the finite-``h`` leverage corrections are smallest); the leading-order
  closed form; positivity, tail growth, off-support `NaN`, general ``g``, scalar ``=`` vector.

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

## Conditional & transition density (§4.3) — `conditional_density`, `transition_density`

- **Estimand.** ``f_{Y\mid X}(y\mid x) = f_{X,Y}(x,y)/f_X(x)``, and the two-time transition density
  ``p_\tau(y\mid x) = f_{X_{t+\tau}\mid X_t}(y\mid x)`` (the same object on a lag-``\tau`` joint).
- **Estimator.** the binned double-kernel (Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996), with a
  conditioning kernel ``K_1 =`` `kernelx` (bandwidth ``h_1 =`` `hx`, in ``x``) and a response kernel
  ``K_2 =`` `kernely` (bandwidth ``h_2 =`` `hy`, in ``y``); each ``\int K = 1``, default `gaussian_kernel`:
  from the cell counts ``n_{jk}`` at centers ``(c^x_j, c^y_k)`` with ``N_j = \sum_k n_{jk}``,

  ```math
  \hat f(y\mid x) =
    \frac{\sum_{j,k} n_{jk}\, K_{h_1}(c^x_j - x)\, K_{h_2}(c^y_k - y)}{\sum_j N_j\, K_{h_1}(c^x_j - x)} ,
  ```

  computed as ``\hat f(y\mid x) = \big(\sum_k P_k\, K((c^y_k - y)/h_2)\big) / (h_2 D)`` for a symmetric
  ``K_2`` (the surviving ``1/h_2`` from ``K_{h_2}``), with the conditioning-kernel-weighted profile
  ``P_k = \sum_j n_{jk} K_{h_1}(c^x_j - x)`` and ``D = \sum_k P_k``. The conditioning-kernel weighting
  cancels in the ratio (so `kernelx` only weights the local window); for a symmetric ``K_2`` the
  ``1/h_2`` makes ``\int \hat f(y\mid x)\,\mathrm dy = 1`` exactly (``\int K = 1``). For a **bounded /
  positive response** an **asymmetric** ``K_2`` = `kernely` (gamma / beta; Chen 1999, 2000) is instead a
  density in the data argument, ``\hat f = \big(\sum_k P_k\, K_2(c^y_k;\, y, h_2)\big)/D`` — boundary-bias-free
  in ``y`` but only ``\approx``-normalized, made exact by the `normalize` renormalization.
- **Consistency claim** (``h_1, h_2 \to 0``, ``n h_1 h_2 \to \infty``): for symmetric kernels the pointwise bias is ``\tfrac12 h_1^2 \mu_2(K)\big[\partial_{xx} f_{Y\mid X} + 2\,\partial_x f_{Y\mid X}\,\partial_x f_X / f_X\big] + \tfrac12 h_2^2 \mu_2(K)\,\partial_{yy} f_{Y\mid X}`` and the variance ``R(K)^2 f_{Y\mid X}/(n h_1 h_2 f_X)``, giving the slower conditional rate ``O(n^{-2/3})`` at the optimal ``h \sim n^{-1/6}`` (the bias is worst at the conditional mode). An asymmetric ``K_2`` carries the boundary-uniform ``O(h_2)`` Chen ``y``-bias instead of the ``O(h_2^2)`` interior term — the price of consistency at the ``y = 0`` edge a symmetric kernel cannot reach.
- **Test.** `test/conditionaldensity.jl` — the bivariate-Gaussian fixture ``Y\mid X=x \sim N(\rho x, 1-\rho^2)``: ``\int \hat f\,\mathrm dy \approx 1`` and pointwise recovery to ``3\times10^{-2}``; the
  scalar/vector forms, the `transition_density` equivalence, off-support `NaN`, bandwidth checks, and
  mergeability; plus the **asymmetric-response** oracle (``X\sim N(0,1)``, ``Y\sim\mathrm{Exp}(1)\perp X``
  ``\Rightarrow`` ``f(y\mid x)=e^{-y}``), where the gamma `kernely` recovers the density at the boundary
  ``y\approx0`` that the symmetric Gaussian leaks low, `normalize` driving ``\int\hat f\,\mathrm dy`` to 1.
  `test/oracles.jl` — the Ornstein–Uhlenbeck oracle: `transition_density` recovers the
  closed-form ``p_\tau(y\mid x) = N(x e^{-\theta\tau},\, s^2(1-e^{-2\theta\tau}))`` to ``4\times10^{-2}``
  (the ``O(h^2)`` mode bias). `test/conditional_convergence.jl` — the ``n^{-2/3}`` rate: on the
  bivariate-Gaussian oracle, with ``h_1 = h_2 \propto n^{-1/6}``, the ISE in ``y`` at ``x_0 = 0`` falls
  with fitted log–log slope ``\approx -2/3`` (reaching the 2-D asymptotic regime needs ``n \gtrsim 10^5``,
  so the n-range starts higher than the 1-D density test).


## Conditional density standard error (§4.3) — `conditional_density_standard_error`, `conditional_density_with_standard_error`

- **Estimand.** the pointwise standard error of ``\hat f(y\mid x)``, ``\mathrm{SE}(y\mid x) = \sqrt{\operatorname{Var}\hat f(y\mid x)}`` with the two-term double-kernel variance ``\operatorname{Var}\hat f(y\mid x) \approx R(K_1)R(K_2) f_{Y\mid X}(y\mid x)/(n h_1 h_2 f_X(x)) - f_{Y\mid X}(y\mid x)^2/n`` (Hyndman, Bashtannyk & Grunwald 1996; the ``-f^2/n`` term the double-kernel analog of the marginal Wand–Jones correction).
- **Estimator.** ``\sqrt{R(K_1)R(K_2)\,\hat f(y\mid x)/(h_1 h_2 D) - \hat f(y\mid x)^2/n}`` from the recovered ``\hat f(y\mid x)``, the roughnesses ``R(K) = \int K^2`` (by quadrature), and the ``x``-window denominator ``D = \sum_j N_j K_{h_1}(c^x_j - x) = n\,\hat f_X(x)`` the estimator already forms — so ``n\,\hat f_X(x)`` is free. It is the two-dimensional analog of `density_standard_error`, the ``x``-window effective sample size ``n h_1 \hat f_X(x)`` replacing ``n``. `transition_density_standard_error` is the same on a lag-``\tau`` joint; `conditional_density_with_standard_error` returns the ``(\texttt{value}, \texttt{standard\_error})`` `NamedTuple` from one profile. Symmetric kernels (the asymmetric-``K_2`` variance a recorded gap).
- **Consistency claim.** it estimates the leading double-kernel sampling variance (not the ``O(h^2)`` smoothing bias); the two-term form matches the finite-sample spread, and the relative error ``\mathrm{SE}/\hat f`` is ``O\bigl((n h_1 h_2)^{-1/2}\bigr)`` — the 2-D effective sample size.
- **Test.** `test/conditionaldensity.jl` — the SE equals the two-term closed form (with the true ``f``, ``f_X``) to ``6\%``; a **Monte-Carlo** check (predicted SE ``=`` the empirical spread of ``\hat f(0\mid 0)`` across independent datasets to ``20\%``); the exact ``n^{-1/2}`` copy-merge scaling and the asymptotic ``h_2^{-1/2}`` bandwidth scaling; the paired reader and the `transition_density` aliases match the separate calls; off-support `NaN`, and the asymmetric-kernel gap throws.

## Power-sum quadrature — expectation of a functional (§4.B) — `expectation`

- **Estimand.** the expectation of a smooth functional of ``X``,
  ``\mathbb{E}[g(X)] = \int g(x)\, f_X(x)\,\mathrm dx``.
- **Estimator.** the power-sum quadrature ``\hat{\mathbb{E}}[g] = \tfrac1n \sum_j \sum_{p=0}^{\texttt{order}}
  a_p^{(j)}\, S_p^{(j)}`` (``S_0^{(j)} = n_j``), where ``a_p^{(j)}`` are the coefficients (in ``u = x - c_j``)
  of the degree-``\texttt{order}`` polynomial interpolant of ``g`` on bin ``j`` — a Chebyshev–Lobatto
  interpolant, its inverse Vandermonde precomputed once (uniform bins), so no derivatives of ``g`` are
  formed. For a monomial ``g = x^m`` (``m \le \texttt{order}``) the interpolant is exact and ``a_p^{(j)}``
  is the binomial ``\binom{m}{p} c_j^{m-p}``, so the readout **reduces to the shift formula** — `mean`,
  `var`, `central_moment` are this estimator, and share the `Expectations._quadrature` atom.
- **Consistency claim.** exact (to floating-point) for polynomial ``g`` of degree ``\le \texttt{order}``,
  independent of the bin width; ``O(\Delta^{\,\texttt{order}+1})`` for general smooth ``g`` (the
  interpolation remainder). The estimand needs ``g`` smooth on the bin scale and the mass inside the range.
- **Test.** `test/expectations.jl` — the Gaussian oracles ``\mathbb{E}[\cos X] = e^{-1/2}``,
  ``\mathbb{E}[\sin X] = 0`` (MGF/trig); polynomial exactness (``\hat{\mathbb{E}}[g] = \tfrac1n\sum_i
  g(X_i)`` to ``10^{-10}`` for degree ``\le \texttt{order}``); the binning error vs the raw empirical mean
  ``< 10^{-6}`` for a smooth ``g``; the moment unification (``\mathbb{E}[x^p]`` equals the recovered moment);
  mergeability; empty ``\Rightarrow`` `NaN`.

## Power-sum quadrature — sampling standard error (§4.B) — `expectation_standard_error`

- **Estimand.** the sampling standard error of ``\hat{\mathbb{E}}[g]``,
  ``\mathrm{SE} = \sqrt{\operatorname{Var}(g(X))/n}`` with ``\operatorname{Var}(g(X)) = \mathbb{E}[g^2] -
  \mathbb{E}[g]^2``.
- **Estimator.** ``\sqrt{(\hat{\mathbb{E}}[g^2] - \hat{\mathbb{E}}[g]^2)/n}`` — both moments the same
  quadrature (``\mathbb{E}[g^2]`` is `expectation` with the smooth integrand ``g^2``). It is the plug-in
  sampling error of the empirical mean of ``g(X)``; the ``O(\Delta^{\,\texttt{order}+1})`` quadrature bias
  is separate and not included. The paired `expectation_with_standard_error` returns this error with the
  estimate ``\hat{\mathbb{E}}[g]`` from one shared quadrature pass (the
  ``(\texttt{value}, \texttt{standard\_error})`` `NamedTuple`).
- **Consistency claim.** for ``g = \mathrm{identity}`` it is exactly ``s/\sqrt n`` (``s`` the uncorrected
  standard deviation); in general the plug-in ``\widehat{\operatorname{Var}(g(X))}`` is consistent as
  ``\Delta \to 0``.
- **Test.** `test/expectations.jl` — ``\mathrm{SE}(\mathrm{identity}) = \operatorname{std}_{\text{unc}}/\sqrt
  n`` to ``10^{-10}``; ``\mathrm{SE}(\cos)`` matches the raw empirical ``\sqrt{\operatorname{Var}(\cos X)/n}``;
  positivity; empty ``\Rightarrow`` `NaN`.

## Power-sum quadrature — pushforward density (§4.B) — `pushforward_density`

- **Estimand.** the density of the transformed variable ``Y = g(X)``, the pushforward
  ``f_{g(X)}(y) = \sum_{x:\,g(x)=y} f_X(x)/|g'(x)|`` (change of variables, summed over preimages).
- **Estimator.** the kernel density estimate of the pushforward written as an expectation,
  ``\hat f_{g(X)}(y) = \tfrac1n \sum_i K_h(g(X_i) - y) = \hat{\mathbb{E}}[K_h(g(X) - y)]``, so it is one
  `expectation` of the composed test function ``x \mapsto K_h(g(x) - y)`` — no ``g^{-1}`` and no explicit
  Jacobian. A non-monotone ``g`` contributes from every preimage automatically; ``g' = 0`` is an integrable
  caustic. For a symmetric ``K`` (``\int K = 1``) the surviving ``1/h`` gives ``\int \hat f\,\mathrm dy = 1``
  exactly.
- **Consistency claim.** the standard KDE bias–variance in ``y`` (bias ``\tfrac12 h^2 \mu_2(K) f_{g(X)}''``,
  variance ``R(K) f_{g(X)}/(nh)``), provided ``g`` is smooth on the bin scale and ``h \gtrsim |g'|\,\Delta``
  (so ``K_h\circ g`` is resolved across a bin).
- **Test.** `test/expectations.jl` — ``Y = e^X`` with ``X\sim N(0,1)`` recovers the Lognormal``(0,1)`` pdf
  (pointwise to ``6\%``, ``\int\hat f\,\mathrm dy \approx 1``), and the **non-monotone** ``Y = X^2`` recovers
  the ``\chi^2_1`` pdf in the bulk (both preimages ``\pm\sqrt y`` contributing); scalar ``=`` vector,
  positive-``h`` check, mergeability.
