# Recovered-indicator study: elements, meshes, fields, recovery, estimators,
# elasticity and relocation. Included from `runtests.jl`, which already loads
# the packages. Every tolerance is a structural identity, so the suite stays
# meaningful if the quadrature or the solver changes.

@testset "reference shape functions reproduce the nodes" begin
    x, w = gauss_legendre(7)
    @test sum(w) ≈ 2 atol = 1e-14
    @test sum(w .* x .^ 12) ≈ 2 / 13 atol = 1e-13
    for kind in ("Q4", "Q8", "Q9")
        nodes = basis(kind)[1]
        N = shapes(kind, nodes)[1]
        @test N ≈ Matrix(I, size(nodes, 1), size(nodes, 1)) atol = 1e-12
        q = [0.3 -0.2; -0.7 0.5]
        @test vec(sum(shapes(kind, q)[1]; dims=2)) ≈ ones(2) atol = 1e-12
        # Shape derivatives reproduce a linear field exactly.
        D = shapes(kind, q)[2]
        @test D[1, :, 1]' * nodes[:, 1] ≈ 1 atol = 1e-12
        @test D[1, :, 2]' * nodes[:, 1] ≈ 0 atol = 1e-12
    end
    @test length(EXP) == 10
    @test HIGH == collect(4:10)
end

@testset "structured meshes have unit area and a positive Jacobian" begin
    for kind in ("Q4", "Q8", "Q9")
        m = structured_mesh(4; kind=kind, family="regular")
        @test sum(mesh_geometry(m, 5)[4]) ≈ 1 atol = 1e-12
        @test ncells(m) == 16
        @test all(corner_determinants(m) .> 0)
    end
    # The jitter depends only on its arguments.
    a = structured_mesh(8; kind="Q4", family="random", amplitude=0.2, seed=7)
    b = structured_mesh(8; kind="Q4", family="random", amplitude=0.2, seed=7)
    other = structured_mesh(8; kind="Q4", family="random", amplitude=0.2, seed=19)
    @test a.points == b.points
    @test other.points != a.points
    # Curved sides need a Bernstein certificate, and that certificate is positive.
    m = structured_mesh(4; kind="Q8", family="curved_A", amplitude=0.15)
    lower, tag = certify(m)
    @test tag == "Bernstein-degree3"
    @test lower > 0
    @test_throws ErrorException structured_mesh(4; kind="Q4", family="curved_A")
end

@testset "manufactured cubic field matches its Taylor coefficients" begin
    u = manufactured("cubic")
    c = field_coefficients(u, 0.0, 0.0)
    @test c ≈ [0, 1, -1, 0, 0, 0, 1, 2, -0.6, 1] atol = 1e-10
    # The load is minus the Laplacian.
    H = field_hessian(u, 0.3, 0.7)
    @test field_force(u, 0.3, 0.7) ≈ -(H[1, 1] + H[2, 2]) atol = 1e-12
end

@testset "patch test reproduces polynomials in the element space" begin
    for (kind, family) in (("Q4", "random"), ("Q8", "curved_B"), ("Q9", "checkerboard"))
        m = structured_mesh(4; kind=kind, family=family, amplitude=0.2, seed=7)
        uh, res = solve_field(m, "affine"; order=7)
        EI, Eh, _ = field_errors(m, uh, "affine"; order=9)
        @test res < 1e-10
        @test sqrt(sum(Eh)) < 1e-10
        @test sqrt(sum(EI)) < 1e-10
    end
    # Q9 reproduces a quadratic exactly. Q4 does not, so this check stays on Q9.
    m9 = structured_mesh(4; kind="Q9", family="checkerboard", amplitude=0.4)
    uh, _ = solve_field(m9, "quadratic"; order=7)
    @test sqrt(sum(field_errors(m9, uh, "quadratic"; order=9)[2])) < 1e-10
end

@testset "recovered Gram is the quadratic form of the high modes" begin
    m = structured_mesh(4; kind="Q4", family="random", amplitude=0.2, seed=7)
    u = manufactured("cubic")
    vals = [u(m.points[i, 1], m.points[i, 2]) for i in 1:npoints(m)]
    rec = element_recovery(m, vals)
    @test physical_coefficients(rec) ≈ coefficients_at(u, rec.centers) atol = 1e-10
    # The Gram operator is symmetric and positive semidefinite.
    G, R, w = gram(m, rec.centers; order=7)
    for e in 1:ncells(m)
        Ge = G[e, :, :]
        @test Ge ≈ Ge' atol = 1e-12
        @test minimum(eigvals(Symmetric(Ge))) > -1e-12
    end
    # The contraction of a coefficient set equals the quadratic form it defines.
    C = physical_coefficients(rec)
    direct = contraction(C, R, w)
    quad = [dot(C[e, HIGH], G[e, :, :] * C[e, HIGH]) for e in 1:ncells(m)]
    @test direct ≈ quad atol = 1e-10
    # The budget is a seminorm of the coefficient difference, so it bounds the gap.
    uh, _ = solve_field(m, "paper1"; order=7)
    s, _, _ = diagnostics(m, uh, "paper1"; order=7)
    @test s["budget_inequality_slack"] >= -1e-12
end

@testset "residual and ZZ indicators vanish on an affine solution" begin
    m = structured_mesh(4; kind="Q4", family="regular")
    uh, _ = solve_field(m, "affine"; order=7)
    @test sqrt(sum(zz_indicator(m, uh; order=7))) < 1e-10
    @test sqrt(sum(residual_indicator(m, uh, "affine"; order=7))) < 1e-10
end

@testset "elastic probes satisfy Navier and the patch Pythagorean identity" begin
    for name in ("affine", "quadratic_cos", "quadratic_sin", "cubic_cos", "cubic_sin")
        r = 0.0
        for x in (0.2, 0.8), y in (0.3, 0.9)
            r = max(r, maximum(abs, navier_residual(name, x, y)))
        end
        @test r < 1e-12
    end
    C = navier_basis()
    @test size(C) == (2, 10, 8)
    μ, λ = 5 / 13, 15 / 26
    for mode in 1:8
        f1(z) = sum(C[1, k, mode] * z[1]^EXP[k][1] * z[2]^EXP[k][2] for k in 1:10)
        f2(z) = sum(C[2, k, mode] * z[1]^EXP[k][1] * z[2]^EXP[k][2] for k in 1:10)
        H1 = ForwardDiff.hessian(f1, [0.3, 0.7])
        H2 = ForwardDiff.hessian(f2, [0.3, 0.7])
        @test abs(μ * (H1[1, 1] + H1[2, 2]) + (λ + μ) * (H1[1, 1] + H2[1, 2])) < 1e-12
        @test abs(μ * (H2[1, 1] + H2[2, 2]) + (λ + μ) * (H1[1, 2] + H2[2, 2])) < 1e-12
    end
    m = structured_mesh(4; kind="Q4", family="random", amplitude=0.2, seed=7)
    sys = PatchSystem(m; order=11)
    uh, info = patch_solve(sys, "affine")
    s, _ = assess_patch(m, uh, "affine"; order=13)
    @test info["algebraic_residual"] < 1e-10
    @test s["E_h"] < 1e-10
    # The oracle Gram value is the interpolation error of a non-affine probe.
    uh2, _ = patch_solve(sys, "cubic_sin")
    s2, _ = assess_patch(m, uh2, "cubic_sin"; order=13)
    @test s2["oracle"] ≈ s2["E_I"] rtol = 1e-8
    @test s2["pythagoras_absolute_defect"] < 1e-12
end

@testset "critical mode is the worst geometric response of its space" begin
    kinds = ("Q4", "Q8", "Q9")
    meshes = Dict(
        k => structured_mesh(4; kind=k, family="checkerboard", amplitude=0.8, seed=7) for k in kinds
    )
    modes = Dict{String,Any}()
    for k in kinds
        spec = geometric_response(meshes[k]; order=13)
        mode, s, _, _ = critical_mode(
            meshes[k], 6, spec["G"][6, :, :], spec["H"][6, :, :]; order=21,
        )
        modes[k] = (mode, s)
        # The eigenvalue and a direct integration of the same mode agree.
        @test s["q_direct"] ≈ s["q_eigen"] rtol = 1e-7
        @test s["normalization"] ≈ 1 atol = 1e-10
        @test s["eigen_residual"] < 1e-12
        @test s["relative_top_gap"] > 0
        # The degree decomposition is exact.
        @test s["error2"] ≈ s["quadratic_error2"] + s["cubic_error2"] + s["cross_error2"] atol = 1e-12
    end
    # A mode taken from one space cannot beat the worst case of the space it is evaluated in.
    for src in kinds, tgt in kinds
        @test direct_response(meshes[tgt], 6, modes[src][1]) <= modes[tgt][2]["q_direct"] + 1e-7
    end
end

@testset "node relocation descends and keeps the boundary fixed" begin
    m = curved_beam(6; kind="Q4", amplitude=0.3, seed=11)
    @test all(corner_determinants(m) .> 0)
    uh, _ = solve_field(m, "beam"; order=7)
    _, _, rec = diagnostics(m, uh, "beam"; order=7)
    obj, _ = relocation_objective(m, rec, "gcee"; order=5)
    z = zeros(2 * count(.!m.boundary))
    # The AD gradient matches a central difference.
    g = ForwardDiff.gradient(obj, z)
    h = 1e-6
    for k in (1, 7, length(z))
        e = zeros(length(z))
        e[k] = h
        @test (obj(z + e) - obj(z - e)) / (2h) ≈ g[k] rtol = 1e-5
    end
    # Descent is monotone and every accepted mesh stays valid.
    m2, hist, _ = relocate_nodes(m, rec; mode="gcee", iterations=6, order=5)
    @test hist[end].J <= hist[1].J
    @test all(step.det_min > 0 for step in hist)
    @test all(corner_determinants(m2) .> 0)
    # Boundary nodes never move.
    @test m2.points[m.boundary, :] ≈ m.points[m.boundary, :] atol = 1e-14
end
