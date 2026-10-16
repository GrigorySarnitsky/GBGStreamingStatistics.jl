# Notation

Each symbol below denotes one concept across the documentation, and each concept has one symbol. A
docstring that uses a symbol differently is the one that changes.

## The samples

| Symbol | Meaning |
|:--|:--|
| ``X``, ``Y`` | the variate (the conditioning variable of a joint summary) and the response |
| ``X_i``, ``Y_i`` | the ``i``-th sample, or the ``i``-th pair ``(X_i, Y_i)`` |
| ``n`` | the sample count |
| ``x``, ``y``, ``z`` | the query point, the response point, and the transformed-response point |
| ``\tau`` | the lag of a two-time pair ``\bigl(X(t), Y(t+\tau)\bigr)`` |

## The summary

| Symbol | Meaning |
|:--|:--|
| ``\Delta`` | the bin width (``\Delta_x``, ``\Delta_y`` per axis of a joint summary) |
| ``c_j`` | the center of bin ``j`` (``c^x_j``, ``c^y_k`` per axis of a joint summary) |
| ``n_j`` | the count of bin ``j``; in the local fit, its design weight |
| ``n_{jk}``, ``N_j`` | the count of cell ``(j, k)`` of a joint summary, and of its ``x``-column ``j`` |
| ``S_p^{(j)}`` | the power-sum of order ``p`` about the center of bin ``j``, ``\sum_{i \in j} (X_i - c_j)^p`` |
| ``\tilde c_j`` | the location of bin ``j``'s mass: ``c_j``, or the centroid ``c_j + S_1^{(j)}/n_j`` |
| ``r_j`` | the response total of bin ``j`` in the local fit |

## Smoothing kernels

| Symbol | Meaning |
|:--|:--|
| ``K`` | a smoothing kernel, ``\int K = 1`` |
| ``h`` | the bandwidth of a symmetric kernel (``h_1``, ``h_2`` or ``h_z`` per direction) |
| ``K_h`` | the scaled kernel, ``K_h(u) = K(u/h)/h`` |
| ``b`` | the scale of the gamma or beta kernel, and nothing else |
| ``w(x)`` | the local-scale weight of the gamma or beta kernel |
| ``R(K)`` | the kernel roughness, ``\int K^2`` |

## The local fit

| Symbol | Meaning |
|:--|:--|
| ``d`` | the degree of a local-polynomial fit |
| ``k`` | the power index of the local fit's moments ``T_k``, ``P_k`` |
| ``\beta_0, \dots, \beta_d`` | the jet: the local-polynomial coefficients at ``x`` |
| ``\omega_j(x)`` | the effective weight the fitted value ``\beta_0`` places on column ``j`` |
| ``\ell(x)`` | the leverage of the local fit, ``\sum_j N_j \omega_j(x)^2`` |

## Estimands

| Symbol | Meaning |
|:--|:--|
| ``f_X`` | the marginal density of ``X`` |
| ``f_{Y \mid X}`` | the conditional density of ``Y`` given ``X`` |
| ``p_\tau(y \mid x)`` | the transition density at lag ``\tau`` |
| ``m(x)`` | the conditional mean, ``\mathbb{E}[Y \mid X = x]`` |
| ``\sigma^2(x)`` | the conditional variance, ``\operatorname{Var}(Y \mid X = x)`` |
| ``\mathbb{E}[g(X)]`` | the expectation of a smooth functional ``g`` |
| ``f_{g(X)}`` | the density of the transformed variable ``g(X)`` |
| ``\pi_k`` | the mixing weight of component ``k`` in a parameter mixture |
