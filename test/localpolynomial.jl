using Test
using GBGStreamingStatistics

const LP = GBGStreamingStatistics.LocalPolynomial
const gk = gaussian_kernel

@testset "LocalPolynomial.local_fit — the local-polynomial jet" begin
  # When the per-bin response mean r_j/w_j is exactly a polynomial P in the bin location, the
  # kernel-weighted least-squares fit recovers its jet exactly: β = (P(x), P'(x), P''(x)/2),
  # independent of the (positive) weights and the kernel — the model is satisfied with zero residual
  # on a non-degenerate design. This pins the value β₀ = m̂, the slope β₁ = m', and the curvature
  # 2β₂ = m'' that the response-aware local bandwidth h⋆(x) reads.
  centers = collect(-2.0:0.1:2.0)
  mass = fill(3.0, length(centers))             # arbitrary positive design weights
  P(c)   = 1.0 + 2.0 * c + 3.0 * c^2
  dP(c)  = 2.0 + 6.0 * c
  ddP    = 6.0
  response = P.(centers) .* mass

  @testset "degree 2 recovers (value, slope, ½·curvature)" begin
    for x in (-1.0, 0.0, 0.5, 1.3)
      β = LP.local_fit(centers, mass, response, x, gk, 2, 0.8)
      @test β[1] ≈ P(x)       atol = 1e-6     # β₀ = m̂(x)
      @test β[2] ≈ dP(x)      atol = 1e-6     # β₁ = m'(x)
      @test β[3] ≈ ddP / 2    atol = 1e-6     # 2β₂ = m''(x)
    end
  end

  @testset "degree 1 recovers (value, slope) of a linear mean" begin
    Q(c) = -0.5 + 1.5 * c
    resp1 = Q.(centers) .* mass
    for x in (-1.0, 0.0, 1.0)
      β = LP.local_fit(centers, mass, resp1, x, gk, 1, 0.8)
      @test β[1] ≈ Q(x)  atol = 1e-7
      @test β[2] ≈ 1.5   atol = 1e-7
      @test length(β) == 2
    end
  end

  @testset "degree 0 is the kernel-weighted (Nadaraya–Watson) ratio" begin
    β = LP.local_fit(centers, mass, response, 0.4, gk, 0, 0.8)
    nw = sum(gk((centers[j] - 0.4) / 0.8) * response[j] for j in eachindex(centers)) /
         sum(gk((centers[j] - 0.4) / 0.8) * mass[j]     for j in eachindex(centers))
    @test length(β) == 1
    @test β[1] ≈ nw  rtol = 1e-12
  end

  @testset "degenerate design ⇒ all-NaN jet" begin
    # A query point with no local support (all kernel weight underflows) cannot be fit.
    β = LP.local_fit(centers, mass, response, 1.0e6, gk, 1, 0.5)
    @test all(isnan, β)
  end

  @testset "kernel normalization cancels in the solve" begin
    raw(u) = exp(-u^2 / 2)                       # the unnormalised Gaussian (no 1/√(2π))
    for x in (-0.7, 0.3)
      βn = LP.local_fit(centers, mass, response, x, gk, 2, 0.8)
      βr = LP.local_fit(centers, mass, response, x, raw, 2, 0.8)
      @test βn ≈ βr  rtol = 1e-10
    end
  end
end
