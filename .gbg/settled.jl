# .gbg/settled.jl — what the owner has already decided, and agents are told not to report.
#
# Read as data by ConfigIO, which evaluates nothing. Written by `gbg answer`; no build, test
# or check reads it, and only the agents of a run are shown it.
[
    Settled(question = "2026-10-15-documentation-review-01-014-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "docs/src/api.md", node = "API overview", text = "the upper half ([`Expectations`](Expectations/) through [`Bandwidth`](Bandwidth/)) holds the estimands", line = 10),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-018-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "src/Accumulators.jl", node = "merge!", text = "tested as accumulate-then-merge `=` accumulate-together", line = 89),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-071-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "src/MarginalDensity.jl", node = "density", text = "it takes its own scale \$b\$ in place of `h`.", line = 439),
                Location(file = "src/MarginalDensity.jl", node = "density", text = "\$\\hat f(x) \\big/ \\int_a^b K_h(x-u)\\,\\mathrm{d}u\$", line = 466),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-072-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "src/MarginalDensity.jl", node = "density", text = "of the conditioning variable \$X\$, defined by \$P(X \\in A) = \\int_A f_X\\,\\mathrm{d}\\lambda\$.", line = 416),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-077-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting)", text = "Estimand. the AMISE-optimal bandwidth (AMISE = asymptotic mean integrated squared error)", line = 410),
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting)", text = "Estimator. for a marginal, the kernel- and degree-aware normal-reference rule", line = 412),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-080-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator", text = "Accumulation knows nothing of bandwidth, kernel, or which response channel will be asked for.", line = 92),
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half", text = "raise \$P\$ for the heavy-tailed channels", line = 222),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-081-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half", text = "for the OU oracle both marginals are the same stationary Gaussian", line = 282),
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy)", text = "basis projection à la Stochastic Force Inference / ULI", line = 162),
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF)", text = "MD's van Hove function is `ConditionalDensity` with a different sampler", line = 620),
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 11. Decision log", text = "so every readout — count sum, NW ratio, jet — is uniform", line = 897),
                Location(file = "docs/ARCHITECTURE.md", node = "GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half", text = "the conditional increment moments of the Pawula \$Q\$ instead need the `Histogram2D` order-raise", line = 222),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-097-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "docs/src/tutorial.md", node = "Getting started / The workflow", text = "a finely binned empirical measure of a conditioning variable \$X\$ — the canonical summary.", line = 17),
            ]),
    Settled(question = "2026-10-15-documentation-review-01-011-code",
            text = "not a defect",
            date = "2026-10-19",
            locations = [
                Location(file = "src/Bandwidth.jl", node = "bandwidth", text = "R(K) &= \\\\int K^2, \\\\qquad \\\\mu_2(K) = \\\\int u^2 K(u)\\\\,\\\\mathrm{d}u,", line = 337),
                Location(file = "src/Accumulators.jl", node = "merge!", text = "\\\\qquad w_k = f(\\\\beta_k) \\\\Delta\\\\beta,", line = 76),
                Location(file = "src/ConditionalDensity.jl", node = "conditional_density_standard_error", text = "\\\\qquad R(K) = \\\\int K^2 .", line = 256),
            ]),
]
