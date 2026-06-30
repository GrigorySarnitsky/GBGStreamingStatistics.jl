using Test
using GBGStreamingStatistics

# The shared symmetric smoothing kernels, checked directly against the closed-form constants that set
# every bandwidth and bias/variance factor: the normalization ∫K = 1, the symmetry ∫uK = 0, the second
# moment μ₂(K) = ∫u²K, and the roughness R(K) = ∫K². Integrated here by an independent fine-grid
# trapezoidal quadrature — NOT the package's internal `_quad` — so this is a genuine oracle for the kernel
# definitions, not a self-consistency check. (The δ(K) = [R(K)/μ₂(K)²]^{1/5} canonical factor these feed
# is exercised end-to-end through the Epanechnikov bandwidth ratio in test/bandwidth.jl.)
@testset "Kernels: normalization and moments (∫K=1, μ₂(K), R(K))" begin
  # trapezoidal ∫_{-L}^{L} g on a uniform grid (endpoints halved); independent of the package
  quad(g, L, npts) = let xs = range(-L, L; length = npts)
    (sum(g, xs) - (g(-L) + g(L)) / 2) * step(xs)
  end

  @testset "gaussian_kernel — standard normal density" begin
    @test quad(gaussian_kernel, 12.0, 200_001) ≈ 1.0                         atol = 1e-9   # ∫K = 1
    @test quad(u -> u * gaussian_kernel(u), 12.0, 200_001) ≈ 0.0             atol = 1e-12  # ∫uK = 0
    @test quad(u -> u^2 * gaussian_kernel(u), 12.0, 200_001) ≈ 1.0           atol = 1e-7   # μ₂(K) = 1
    @test quad(u -> gaussian_kernel(u)^2, 12.0, 200_001) ≈ 1 / (2 * sqrt(π)) atol = 1e-8   # R(K)
    @test gaussian_kernel(0.0) ≈ 1 / sqrt(2π)                                               # peak height
  end

  @testset "epanechnikov_kernel — compact support, the 3/4(1−u²) parabola" begin
    @test quad(epanechnikov_kernel, 1.0, 200_001) ≈ 1.0                       atol = 1e-6  # ∫K = 1
    @test quad(u -> u * epanechnikov_kernel(u), 1.0, 200_001) ≈ 0.0           atol = 1e-12 # ∫uK = 0
    @test quad(u -> u^2 * epanechnikov_kernel(u), 1.0, 200_001) ≈ 1 / 5       atol = 1e-6  # μ₂(K) = 1/5
    @test quad(u -> epanechnikov_kernel(u)^2, 1.0, 200_001) ≈ 3 / 5           atol = 1e-6  # R(K) = 3/5
    @test epanechnikov_kernel(0.0) ≈ 0.75                                                   # peak height
    @test epanechnikov_kernel(1.5) == 0.0                                                   # zero off |u| < 1
    @test epanechnikov_kernel(-1.5) == 0.0
  end
end
