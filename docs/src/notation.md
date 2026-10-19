# Notation

Each symbol below denotes one concept across the documentation, and each concept one symbol; a
docstring that uses a symbol differently is corrected. This rule binds named quantities and
indices with a fixed role; a summation or dummy index bound and defined inside one formula (the binomial
index ``k`` of the shift formula that re-centers power-sums, a measurable set ``A``) may reuse a letter.

## The samples

| Symbol | Meaning |
|:--|:--|
| ``X``, ``Y`` | the variate (conditioning variable of paired samples) and the response |
| ``X_i``, ``Y_i`` | the ``i``-th sample, or the ``i``-th pair ``(X_i, Y_i)``; ``i`` indexes samples only |
| ``n`` | the sample count |
| ``x``, ``y``, ``z`` | the query point, the response point, and a transformed variable's point |
| ``Z`` | a smoothly transformed variable, ``Z = g(X)`` or ``Z = g(Y)`` |
| ``t``, ``\tau`` | the time and lag of a two-time pair ``\bigl(X(t), Y(t+\tau)\bigr)`` |

Axis names are lower-case: an ``x``-column of bins, its ``y``-quadrature of responses.

## The summary

| Symbol | Meaning |
|:--|:--|
| ``\Delta`` | the bin width (``\Delta_x``, ``\Delta_y`` per axis of a joint summary) |
| ``j``, ``l`` | the bin (or ``x``-column) index and the ``y``-cell index; ``j'`` a second bin |
| ``c_j`` | the center of bin ``j`` (``c^x_j``, ``c^y_l`` per axis of a joint summary) |
| ``\tilde c_j`` | the location of bin ``j``'s mass: ``c_j``, or the centroid ``c_j + S_1^{(j)}/n_j`` (center plus mean offset) |
| ``n_j`` | the count of bin ``j``; in the local fit, its design weight |
| ``n_{jl}``, ``N_j`` | the count of cell ``(j, l)`` of a joint summary, and of its ``x``-column ``j``, ``N_j = \sum_l n_{jl}`` |
| ``B`` | the number of bins, as in a cost ``O(B \log B)`` |
| ``S_p^{(j)}`` | the power-sum of order ``p`` about the center of bin ``j``, ``\sum_{i \in j} (X_i - c_j)^p`` |
| ``S_{pp'}`` | the power-sum of a cell, ``\sum (X_i - c^x_j)^p (Y_i - c^y_l)^{p'}`` |
| ``\bar y_j`` | the mean response of bin ``j``; ``n_j \bar y_j`` is its response total |
| ``\mathbb{P}_X``, ``\mathbb{P}_n`` | the law of ``X``, with ``f_X = \mathrm{d}\mathbb{P}_X/\mathrm{d}x``, and the empirical measure, ``\mathbb{P}_n(A) = \#\{i : X_i \in A\}/n`` |
| ``T`` | the summary map and its values: the summary ``T(n\mathbb{P}_n)``; ``n(T)`` its total count |
| ``\mathcal{A}``, ``\mathcal{B}`` | sample multisets, ``T(\mathcal{A} \uplus \mathcal{B}) = T(\mathcal{A}) \oplus T(\mathcal{B})``, where ``\oplus`` is `merge!` |
| ``\lambda`` | a summary's scalar multiplier: rescaling, weighted merging, the forgetting factor |
| ``\mathcal{R}_h`` | the post-hoc readout that turns a summary into an estimate at a bandwidth |
| ``J`` | the number of merged copies of a summary |

## Smoothing kernels

| Symbol | Meaning |
|:--|:--|
| ``K`` | a smoothing kernel, ``\int K = 1`` |
| ``u``, ``v`` | the standardized kernel argument, and a variable in data units |
| ``K_h`` | the scaled kernel, ``K_h(v) = K(v/h)/h`` |
| ``K_x``, ``K_y``, ``K_z`` | a joint estimator's kernels, per direction |
| ``h`` | a bandwidth (``h_x``, ``h_y``, ``h_z`` per direction) |
| ``h_j``, ``h(x)`` | the Abramson bandwidth of bin ``j``, and the balloon bandwidth at ``x`` |
| ``h^\star(x)``, ``h_{\mathrm{AMISE}}`` | the local AMSE-optimal bandwidth, and the global AMISE-optimal one |
| ``\tilde h`` | a pilot bandwidth |
| ``G`` | the geometric mean of the Abramson pilot values |
| ``b`` | the scale of the gamma or beta kernel, and nothing else |
| ``A``, ``V`` | the bias and variance constants of the gamma or beta kernel's AMISE, ``A b^2 + V b^{-1/2}/n`` |
| ``\rho(x)``, ``w(x)`` | the gamma kernel's shape, and the gamma or beta kernel's local-scale weight |
| ``R(K)`` | the kernel roughness, ``\int K^2`` |
| ``\mu_p(K)`` | the kernel moment, ``\int u^p K(u)\,\mathrm{d}u`` |
| ``\delta(K)`` | the canonical factor, ``[R(K)/\mu_2(K)^2]^{1/5}`` |
| ``C(K,d)`` | the normal-reference constant, ``h = C(K,d) \hat\sigma_X n^{-1/(2q+1)}``, with the estimated standard deviation, the degree and the bias order; ``\delta(K)`` is its kernel part |
| ``K_{(4)}`` | the fourth-order equivalent kernel of the local-quadratic fit |

## The local fit

| Symbol | Meaning |
|:--|:--|
| ``d`` | the degree of a local-polynomial fit |
| ``U_p``, ``V_p``, ``W_p`` | the local fit's design, response and squared-kernel moments |
| ``U``, ``W`` | the Hankel matrices ``[U_{r+s}]_{r,s=0}^d`` and ``[W_{r+s}]_{r,s=0}^d`` |
| ``\beta_0, \dots, \beta_d`` | the jet: the local-polynomial coefficients at ``x`` |
| ``\omega_j(x)`` | the fitted value ``\hat\beta_0``'s weight on each sample of column ``j`` |
| ``\ell(x)`` | the local fit's leverage, ``\sum_j N_j \omega_j(x)^2 = e_0^\top U^{-1} W U^{-1} e_0``, with ``e_0`` the first unit vector, which picks out the fitted value |

## Estimands

| Symbol | Meaning |
|:--|:--|
| ``f_X`` | the marginal density of ``X`` |
| ``f_{Y \mid X}`` | the conditional density of ``Y`` given ``X`` |
| ``p_\tau(y \mid x)`` | the transition density at lag ``\tau`` |
| ``m(x)`` | the conditional mean, ``\mathbb{E}[Y \mid X = x]`` |
| ``\sigma^2(x)`` | the conditional variance, ``\operatorname{Var}(Y \mid X = x)`` |
| ``\mu_4(x)`` | the conditional fourth central moment, ``\mathbb{E}[(Y - m(x))^4 \mid X = x]`` |
| ``\mathbb{E}[g(X)]`` | the expectation of a smooth functional ``g``; ``\mathbb{E}`` is the only expectation |
| ``f_{g(X)}``, ``f_{g(Y) \mid X}`` | a transformed variable's density, marginal and conditional |
| ``P_l`` | the conditional density's ``y``-profile of kernel-weighted counts |
| ``D`` | the conditional density's ``x``-window denominator, ``\sum_j N_j K_{h_x}(c^x_j - x) = n \hat f_X(x)`` |
| ``M_n(x, \tau)`` | the conditional increment moment ``\mathbb{E}[(X(t+\tau) - x)^n \mid X(t) = x]``, whose index is an order, not a count; ``M_n/(n! \tau)`` its finite-``\tau`` Kramers–Moyal coefficient |

## Parameters

| Symbol | Meaning |
|:--|:--|
| ``\theta``, ``f_\theta`` | the parameter of a parametric family and its density |
| ``\pi_k``, ``\pi(\theta)`` | the mixing weight of component ``k`` and the mixing law of a parameter mixture |

A density symbol always carries its subscript, and an estimate its hat: ``\hat f_X``, ``\hat m(x)``.
Where the bandwidth varies, as in ``\operatorname{MISE}(h)``, the estimate at bandwidth ``h`` is written ``\hat f_h``, ``\hat m_h``.

## The moments of the data

| Symbol | Meaning |
|:--|:--|
| ``\bar X`` | the sample mean, an estimate of ``\mathbb{E}[X]`` |
| ``\mu_p(X)``, ``\hat\mu_p`` | the central moment ``\mathbb{E}[(X - \mathbb{E}X)^p]``, and its estimate ``\tfrac{1}{n} \sum_i (X_i - \bar X)^p`` |
| ``\sigma_X``, ``\hat\sigma_X`` | the standard deviation of ``X``, estimated by ``\sqrt{\hat\mu_2}`` |

## Orders

| Symbol | Meaning |
|:--|:--|
| ``p`` | a power or an order: of a power-sum, a central moment, a monomial ``x^p``, a kernel moment |
| ``q`` | the bias order of a bandwidth rule, ``h \propto n^{-1/(2q+1)}`` |
| ``r`` | the order of a derivative, as in ``\psi_r = \int f_X^{(r)} f_X`` |

## Spelling

Operators are upright: ``\operatorname{Var}``, ``\operatorname{bias}``, ``\operatorname{SE}``,
``\operatorname{MISE}``, ``\operatorname{AMISE}``, ``\operatorname{AIC}``, ``\mathrm{d}x``. So are named laws:
``\mathrm{Gamma}``, ``\mathrm{Exp}``, ``\mathrm{t}_\nu``, with ``\mathcal{N}`` for the normal law.
Conditioning is written in full, ``\mathbb{E}[\cdot \mid X = x]``.
