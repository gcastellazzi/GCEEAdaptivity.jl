# The six shortlist benchmarks. Each run relocates a uniform quad mesh and,
# where the benchmark has an exact solution, solves Poisson before and after.

"""One relocation, with optional Poisson errors on the uniform and adapted meshes."""
struct BenchmarkRun
    name::String
    metric::Int
    n::Int
    objective0::Float64
    objective1::Float64
    iterations::Int
    min_area::Float64
    l2_uniform::Union{Nothing,Float64}
    l2_adapted::Union{Nothing,Float64}
    energy_uniform::Union{Nothing,Float64}
    energy_adapted::Union{Nothing,Float64}
end

function _finish(name, metric, patch, xyz0, xyz, history, solution)
    l2u = l2a = eu = ea = nothing
    if solution !== nothing
        uh0 = solve_poisson(xyz0, patch, solution)
        uh1 = solve_poisson(xyz, patch, solution)
        l2u, eu = poisson_errors(xyz0, patch, solution, uh0)
        l2a, ea = poisson_errors(xyz, patch, solution, uh1)
    end
    return BenchmarkRun(
        name, metric, patch.n,
        history[1], history[end], length(history) - 1,
        minimum_area(xyz, patch),
        l2u, l2a, eu, ea,
    )
end

"""
    run_sinusoidal_band(; variant=:size, n=8, iterations=12, order=4, driver=:analytic)

IMR 2019, §6.1. `variant` is `:size` (μ₇, α=7), `:shape` (μ₅₈, β=3) or
`:both` (μ₉, γ=3). Geometry is bilinear, not the published Q₃ mesh, and the
solver is gradient descent.

`driver = :analytic` uses the published target built from `η`.
`driver = :gcee` keeps that metric and replaces `η` with the Q4 Gram indicator
of the same band field.
"""
function run_sinusoidal_band(; variant::Symbol=:size, n::Integer=8, iterations::Integer=12,
        order::Integer=4, driver::Symbol=:analytic)
    patch = unit_square_mesh(n; order=order)
    target, metric = sinusoidal_target(patch, patch.xyz0; variant=variant)
    label = "sinusoidal-$variant"
    if driver == :gcee
        target = gcee_target(patch.xyz0, patch, sinusoidal_eta)
        label = label * "-gcee"
    elseif driver != :analytic
        error("unknown band driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, metric, target; iterations=iterations)
    return _finish(label, metric, patch, patch.xyz0, xyz, history, nothing)
end

"""
    run_analytic_ring(; kind=:paper, n=8, iterations=12, order=4, driver=:analytic)

hr-adaptivity §4.2, r-branch only.

- `:paper` — radii 0.2 and 0.3, μ₇, size 0.001 / 0.01
- `:tid9` — MFEM `-tid 9`, radii 0.15 and 0.35, same sizes, μ₇
- `:tid10` — MFEM `-tid 10`, size plus angular aspect, μ₉

`driver = :gcee` keeps the metric and builds `W` from the Gram indicator of
the unclamped ring field.
"""
function run_analytic_ring(; kind::Symbol=:paper, n::Integer=8, iterations::Integer=12,
        order::Integer=4, driver::Symbol=:analytic)
    patch = unit_square_mesh(n; order=order)
    if kind == :paper
        build = ring_size_target()
        metric = 7
        r1, r2 = 0.2, 0.3
    elseif kind == :tid9
        build = ring_size_target(; r1=0.15, r2=0.35)
        metric = 7
        r1, r2 = 0.15, 0.35
    elseif kind == :tid10
        build = ring_aspect_target()
        metric = 9
        r1, r2 = 0.2, 0.3
    else
        error("unknown ring kind $kind")
    end
    label = "ring-$kind"
    if driver == :analytic
        target = field_target(build)
    elseif driver == :gcee
        target = gcee_target(patch.xyz0, patch, (x, y) -> ring_feature(x, y; r1=r1, r2=r2, β=30.0))
        label = label * "-gcee"
    else
        error("unknown ring driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, metric, target; iterations=iterations)
    return _finish(label, metric, patch, patch.xyz0, xyz, history, nothing)
end

"""
    run_poisson_front(; sharpness=20.0, n=8, iterations=12, order=4, driver=:gradient)

Circular front on `[0, 1]²`. `sharpness = 20` is Kolev et al. 2026, §8.1;
`sharpness = 200` is the hr paper, §4.1.

`driver = :gradient` uses the gradient shape+size matrix (15).
`driver = :gcee` uses the scalar Q4 Gram indicator as the size and aspect
field, frozen on the initial mesh. Both keep metric μ₉.

Uniform squares already have μ₂ = 0 for `W = I`, so a pure shape target does
not move this mesh. The 2026 paper reduces the error with a PDE objective;
that objective is not this relocation. The Poisson errors reported here are
the uniform mesh and the mesh after relocation.
"""
function run_poisson_front(; sharpness::Real=20.0, n::Integer=8, iterations::Integer=12,
        order::Integer=4, driver::Symbol=:gradient)
    patch = unit_square_mesh(n; order=order)
    α = Float64(sharpness)
    solution(x, y) = circular_front(x, y; sharpness=α)
    if driver == :gradient
        grad(x, y) = circular_front_grad(x, y; sharpness=α)
        target = field_target(gradient_shape_size(grad, patch, patch.xyz0))
        label = "poisson-front-α$α"
    elseif driver == :gcee
        target = gcee_target(patch.xyz0, patch, solution)
        label = "poisson-front-gcee-α$α"
    else
        error("unknown Poisson driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, 9, target; iterations=iterations)
    return _finish(label, 9, patch, patch.xyz0, xyz, history, solution)
end

"""
    run_oblique_shock(; slope=20.0, n=8, iterations=12, order=4, driver=:analytic)

Kolev et al. 2026, §8.2, without the PDE-constrained weight. Metric 107.

`driver = :analytic` uses the orientation field `θ = π y (1 − y) cos(2π x)`
at the uniform element area. `driver = :gcee` uses the Gram indicator of the
exact shock. `slope` scales that solution; §8.2 does not give its value.
The published `L²` table is for a different objective and is not reproduced.
"""
function run_oblique_shock(; slope::Real=20.0, n::Integer=8, iterations::Integer=12,
        order::Integer=4, driver::Symbol=:analytic)
    patch = unit_square_mesh(n; order=order)
    α = Float64(slope)
    solution(x, y) = oblique_shock(x, y; slope=α)
    label = "oblique-shock"
    if driver == :analytic
        target = field_target(orientation_target(patch.n))
    elseif driver == :gcee
        target = gcee_target(patch.xyz0, patch, solution)
        label = label * "-gcee"
    else
        error("unknown shock driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, 107, target; iterations=iterations)
    return _finish(label, 107, patch, patch.xyz0, xyz, history, solution)
end

"""
    run_discrete_size(; n=8, iterations=12, order=4, driver=:analytic)

MFEM `-tid 5` / `-mid 94`: ideal shape, size from `size_indicator`, with the
2D small/big ratio 9 of `ConstructSizeGF`.

`driver = :gcee` keeps metric 94 and builds `W` from the Gram indicator of
the unclamped annulus.
"""
function run_discrete_size(; n::Integer=8, iterations::Integer=12, order::Integer=4,
        driver::Symbol=:analytic)
    patch = unit_square_mesh(n; order=order)
    label = "discrete-size"
    if driver == :analytic
        target = field_target(discrete_size_target(patch, patch.xyz0))
    elseif driver == :gcee
        target = gcee_target(patch.xyz0, patch, size_feature)
        label = label * "-gcee"
    else
        error("unknown size driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, 94, target; iterations=iterations)
    return _finish(label, 94, patch, patch.xyz0, xyz, history, nothing)
end

"""
    run_material_interface(; n=8, iterations=12, order=4, ratio=9.0, driver=:analytic)

MFEM `material_indicator_2d` (`-tid 6` indicator), metric μ₇.

The analytic target area uses the IMR two-level rule with the given `ratio`.
The blend weight is `‖∇η‖` normalised by its maximum on the initial mesh. The
aspect ratio grows from 1 to 3 across the interface, and the short axis is the
gradient direction. This is the gradient construction of the hr paper applied
to the miniapp indicator, not a line-by-line port of `DiscreteAdaptTC`.

`driver = :gcee` keeps metric μ₇ and builds `W` from the Gram indicator of
the unclamped interface.
"""
function run_material_interface(; n::Integer=8, iterations::Integer=12, order::Integer=4,
        ratio::Float64=9.0, driver::Symbol=:analytic)
    patch = unit_square_mesh(n; order=order)
    label = "material-interface"
    if driver == :analytic
        target = field_target(material_interface_target(patch, patch.xyz0; ratio=ratio))
    elseif driver == :gcee
        target = gcee_target(patch.xyz0, patch, material_feature)
        label = label * "-gcee"
    else
        error("unknown interface driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, 7, target; iterations=iterations)
    return _finish(label, 7, patch, patch.xyz0, xyz, history, nothing)
end

"""
    run_shortlist(; n=8, iterations=8, order=4)

Run the six shortlist cases at one resolution, each with the published target
and with the GCEE target, and print a table. The Poisson front is run at both
published sharpness values.
"""
function run_shortlist(; n::Integer=8, iterations::Integer=8, order::Integer=4)
    kw = (; n=n, iterations=iterations, order=order)
    runs = BenchmarkRun[
        run_poisson_front(; sharpness=20.0, driver=:gradient, kw...),
        run_poisson_front(; sharpness=20.0, driver=:gcee, kw...),
        run_poisson_front(; sharpness=200.0, driver=:gradient, kw...),
        run_poisson_front(; sharpness=200.0, driver=:gcee, kw...),
        run_analytic_ring(; kind=:paper, driver=:analytic, kw...),
        run_analytic_ring(; kind=:paper, driver=:gcee, kw...),
        run_sinusoidal_band(; variant=:size, driver=:analytic, kw...),
        run_sinusoidal_band(; variant=:size, driver=:gcee, kw...),
        run_oblique_shock(; driver=:analytic, kw...),
        run_oblique_shock(; driver=:gcee, kw...),
        run_material_interface(; driver=:analytic, kw...),
        run_material_interface(; driver=:gcee, kw...),
    ]
    print_runs(runs)
    return runs
end

"""Print a benchmark table."""
function print_runs(runs::AbstractVector{BenchmarkRun})
    println("name                    metric   n    F0          F1          steps   L2 uniform   L2 adapted")
    for r in runs
        l2u = r.l2_uniform === nothing ? "—" : @sprintf("%12.4e", r.l2_uniform)
        l2a = r.l2_adapted === nothing ? "—" : @sprintf("%12.4e", r.l2_adapted)
        @printf("%-22s  %6d  %3d  %10.4e  %10.4e  %5d  %12s  %12s\n",
            r.name, r.metric, r.n, r.objective0, r.objective1, r.iterations, l2u, l2a)
    end
    return nothing
end
