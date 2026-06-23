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

!!! note "Scaffold stage"
    The single seeded entry below fixes the format. Entries are added as each module is implemented.

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
- **Test.** *(to be linked: `test/oracles.jl` — Gaussian/OU oracle; a rate test that the error falls
  as ``n^{-2/5}``.)*
