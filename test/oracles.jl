using Test
using GBGStreamingStatistics

@testset "analytic oracles (Gaussian, Ornstein–Uhlenbeck)" begin
  # Scaffold stage: estimators will be checked against distributions with closed forms — Gaussian
  # marginals and an Ornstein–Uhlenbeck process (drift, diffusion, transition density known exactly)
  # — and against their theoretical convergence rate (docs/TESTING.md). Pending the estimators.
  @test_skip density_matches_gaussian_oracle()
  @test_skip conditional_mean_matches_ou_oracle()
  @test_skip transition_density_matches_ou_oracle()
  @test_skip estimation_error_falls_as_n_to_the_minus_two_fifths()
end
