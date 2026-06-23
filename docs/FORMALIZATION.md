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

## Per-bin central moments — the Pébay recursion (`Accumulators`, `Histogram1D`)

- **Object.** per bin, the count, mean, and unnormalized central moments ``M_p = \sum_{i \in
  \text{bin}}(X_i - \bar X)^p`` for ``p = 2 \dots`` `order` (intrinsically finite-sample — no
  population estimand).
- **Recursion.** the single-sample update (`_observe_moment!`) and the parallel merge
  (`_merge_moment!`) are the Pébay (2008) formulas, accumulated in `Float64`; the single-sample update
  is the ``n_B = 1`` case of the merge, so streaming accumulation and `merge!` are the *same*
  arithmetic.
- **Correctness.** the recursion is *exact* — it computes the same ``M_p`` as a direct two-pass sum,
  to round-off — and the merge gives accumulate-then-merge ``=`` accumulate-together.
- **Test.** `test/histograms.jl` (recovered ``M_2 \dots M_6`` match the direct sum to ``10^{-8}`` and
  the true Gaussian moments to sampling error); `test/mergeability.jl` (merge ``=`` together);
  `test/allocations.jl` (`observe!`/`merge!` allocate 0 bytes).

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
  ``P_\ell = \sum_i K_h(X_i-x)(X_i-x)^\ell Y_i``, read from the per-bin power-sums of the summary.
- **Consistency claim** (interior ``x``, ``h \to 0``, ``nh \to \infty``):

  ```math
  \operatorname{bias}\hat m_h(x) \approx \tfrac12 h^2 \mu_2(K)\, m''(x),
  \qquad
  \operatorname{Var}\hat m_h(x) \approx \frac{R(K)\,\sigma^2(x)}{n\,h\,f_X(x)},
  ```

  giving the optimal rate ``\hat m_{h^\star}(x) - m(x) = O_p(n^{-2/5})`` at ``h^\star \sim n^{-1/5}``.
- **Test.** *(forthcoming with the estimator: an OU oracle, and a rate test that the error falls as
  ``n^{-2/5}``.)*
