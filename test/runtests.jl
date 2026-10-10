using Test
using LinearAlgebra
using ForwardDiff
using GCEEAdaptivity

const GA = GCEEAdaptivity

@testset "metrics vanish on the identity" begin
    T = [1.0 0.0; 0.0 1.0]
    @test mu2(T) ≈ 0 atol = 1e-14
    @test mu7(T) ≈ 0 atol = 1e-14
    @test mu9(T) ≈ 0 atol = 1e-14
    @test mu14(T) ≈ 0 atol = 1e-14
    @test mu55(T) ≈ 0 atol = 1e-14
    @test mu56(T) ≈ 0 atol = 1e-14
    @test mu58(T) ≈ 0 atol = 1e-14
    @test mu77(T) ≈ 0 atol = 1e-14
    @test mu94(T) ≈ 0 atol = 1e-14
    @test mu107(T, T) ≈ 0 atol = 1e-14
end

@testset "shape is scale invariant and size is not" begin
    T = 2.0 * [1.0 0.0; 0.0 1.0]
    @test mu2(T) ≈ 0 atol = 1e-14
    @test mu58(T) ≈ 0 atol = 1e-14
    @test mu7(T) > 1
    @test mu9(T) > 1
    @test mu94(T) ≈ mu2(T) + 1.5 * mu56(T)
end

@testset "unit-square mesh has area one" begin
    patch = unit_square_mesh(4; order=2)
    area = integrate_scalar(patch.xyz0, patch, (x, y) -> 1.0)
    @test area ≈ 1 atol = 1e-12
    @test minimum_area(patch.xyz0, patch) ≈ 1 / 16 atol = 1e-12
end

@testset "element metric vanishes when the target matches the mesh" begin
    patch = unit_square_mesh(2; order=2)
    target = field_target((x, y) -> [0.5 0.0; 0.0 0.5])
    μ = element_metrics(patch.xyz0, patch, 9, target)
    @test all(abs.(μ) .< 1e-10)
    @test tmop_functional(patch.xyz0, patch, 9, target) ≈ 0 atol = 1e-10
end

@testset "equidistributed constant weight reproduces the mesh" begin
    patch = unit_square_mesh(4; order=2)
    build = equidistributed_isotropic((x, y) -> 3.0, patch, patch.xyz0)
    target = field_target(build)
    F = tmop_functional(patch.xyz0, patch, 9, target)
    @test F ≈ 0 atol = 1e-10
end

@testset "sinusoidal band targets" begin
    patch = unit_square_mesh(8; order=4)
    @test sinusoidal_eta(0.0, 0.5) ≈ 2 * tanh(1) atol = 1e-12
    s, Vη, V = GA.two_level_scale(patch, patch.xyz0, sinusoidal_eta, 7.0)
    @test V ≈ 1 atol = 1e-8
    @test 0 < Vη < 1
    @test s > 0
    # Equation (10): the two-level areas fill the element budget.
    @test Vη / s + (V - Vη) / (7 * s) ≈ nelems(patch) atol = 1e-6
    target, metric = sinusoidal_target(patch, patch.xyz0; variant=:both)
    @test metric == 9
    # η = 0 at a corner of the band's complement. y = 0 is far from the layer.
    W = target(1, 1, 0.01, 0.01)
    # Frozen samples ignore the arguments and use the quadrature location.
    # Check the construction directly.
    γ = 3.0
    sγ, _, _ = GA.two_level_scale(patch, patch.xyz0, sinusoidal_eta, γ)
    W0 = [sqrt(γ * sγ) 0.0; 0.0 sqrt(γ * sγ)]
    W1 = [sqrt(γ * sγ) 0.0; 0.0 sqrt(γ * sγ) / γ]
    @test det(W0) ≈ γ * sγ
    @test det(W1) ≈ sγ
    @test W1[1, 1] / W1[2, 2] ≈ γ
end

@testset "ring size" begin
    inside = ring_size_target()(0.5, 0.75)
    outside = ring_size_target()(0.5, 0.5)
    # r = 0.25 is inside [0.2, 0.3]; the origin of the ring is the centre.
    @test det(inside) ≈ 0.001 atol = 1e-6
    # At the centre the tanh tails do not cancel exactly, so ζ is only near `big`.
    @test det(outside) ≈ 0.01 atol = 1e-4
    corner = ring_size_target()(0.0, 0.0)
    @test det(corner) ≈ 0.01 atol = 1e-8
    tid9 = ring_size_target(; r1=0.15, r2=0.35)(0.5, 0.75)
    @test det(tid9) ≈ 0.001 atol = 1e-6
end

@testset "manufactured gradients" begin
    for (x, y) in ((0.2, 0.4), (0.8, 0.1), (0.5, 0.9))
        g = ForwardDiff.gradient(z -> circular_front(z[1], z[2]; sharpness=20.0), [x, y])
        @test [circular_front_grad(x, y; sharpness=20.0)...] ≈ g atol = 1e-10
        gs = ForwardDiff.gradient(z -> oblique_shock(z[1], z[2]; slope=20.0), [x, y])
        @test [oblique_shock_grad(x, y; slope=20.0)...] ≈ gs atol = 1e-10
        gm = ForwardDiff.gradient(z -> begin
                xc = z[1] - 0.5
                yc = z[2] - 0.5
                th = 22.5 * π / 180
                xn = cos(th) * xc + sin(th) * yc
                yn = -sin(th) * xc + cos(th) * yc
                stretch = 1 / cos(th)
                arg = 20 * (yn / stretch) + 3 * sin(3π * (xn / stretch)) + 1
                return tanh(arg)
            end, [x, y])
        @test [material_indicator_grad(x, y)...] ≈ gm atol = 1e-8
    end
end

@testset "Poisson convergence on a smooth solution" begin
    u(x, y) = sin(π * x) * sin(π * y)
    function l2_of(n)
        patch = unit_square_mesh(n; order=3)
        uh = solve_poisson(patch.xyz0, patch, u)
        l2, energy = poisson_errors(patch.xyz0, patch, u, uh)
        return l2, energy
    end
    l2c, ec = l2_of(4)
    l2f, ef = l2_of(8)
    @test l2c / l2f > 3
    @test ec / ef > 1.5
    @test l2f < 0.05
end

@testset "relocation lowers the band objective" begin
    run = run_sinusoidal_band(; variant=:size, n=4, iterations=6, order=2)
    @test run.objective1 < run.objective0
    @test run.min_area > 0
    @test run.iterations >= 1
end

@testset "material interface target keeps positive area" begin
    patch = unit_square_mesh(6; order=3)
    build = material_interface_target(patch, patch.xyz0)
    for (x, y) in ((0.2, 0.3), (0.5, 0.5), (0.8, 0.7), (0.1, 0.9))
        @test det(build(x, y)) > 0
    end
end

@testset "Q4 Gram indicator on a pure square mode" begin
    # u = x² has u_xx = 2. On each element of a 2×2 mesh, h_K = 1/2 and the
    # hourglass vanishes, so the indicator is h_K² / sqrt(3).
    u(x, y) = x^2
    patch = unit_square_mesh(2; order=2)
    η = gcee_indicators(patch.xyz0, patch, u)
    expect = (0.5^2) / sqrt(3)
    @test all(η .≈ expect)
    target = gcee_target(patch.xyz0, patch, u)
    for e in 1:nelems(patch)
        W = target(e, 1, 0.0, 0.0)
        @test det(W) ≈ 0.25 atol = 1e-10
        @test det(W) > 0
    end
end

@testset "Gram indicator peaks on the circular front" begin
    u(x, y) = circular_front(x, y; sharpness=20.0)
    patch = unit_square_mesh(8; order=2)
    η = gcee_indicators(patch.xyz0, patch, u)
    @test maximum(η) / minimum(η) > 5
end

@testset "GCEE relocation of the front stays valid" begin
    run = run_poisson_front(; sharpness=20.0, n=4, iterations=4, order=3, driver=:gcee)
    @test run.objective1 <= run.objective0 * (1 + 1e-8)
    @test run.min_area > 0
    @test run.l2_adapted > 0
    @test run.energy_adapted > 0
end

@testset "GCEE target stays valid on the other shortlist cases" begin
    small = (n = 4, iterations = 3, order = 2)
    runs = [
        run_poisson_front(; sharpness=200.0, driver=:gcee, small...),
        run_analytic_ring(; kind=:paper, driver=:gcee, small...),
        run_analytic_ring(; kind=:tid9, driver=:gcee, small...),
        run_analytic_ring(; kind=:tid10, driver=:gcee, small...),
        run_sinusoidal_band(; variant=:size, driver=:gcee, small...),
        run_sinusoidal_band(; variant=:shape, driver=:gcee, small...),
        run_sinusoidal_band(; variant=:both, driver=:gcee, small...),
        run_oblique_shock(; driver=:gcee, small...),
        run_discrete_size(; driver=:gcee, small...),
        run_material_interface(; driver=:gcee, small...),
    ]
    for run in runs
        @test run.min_area > 0
        @test run.objective1 <= run.objective0 * (1 + 1e-6)
    end
end

@testset "Poisson front relocation stays valid" begin
    run = run_poisson_front(; sharpness=20.0, n=4, iterations=4, order=3)
    @test run.objective1 <= run.objective0 * (1 + 1e-8)
    @test run.min_area > 0
    @test run.l2_uniform > 0
    @test run.l2_adapted > 0
    @test run.energy_uniform > 0
end

include("test_studies.jl")
