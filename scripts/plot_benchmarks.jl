# Mesh figures and convergence for the shortlist.
#
# Each case keeps its published metric. TMOP uses the published target.
# GCEE replaces that target with the scalar Q4 Gram indicator of the same
# feature. The circular front and the oblique shock also report the Poisson
# error on the uniform mesh.
#
#   julia --project=. scripts/plot_benchmarks.jl
#
# Writes figures/ and data/*.csv.

using GCEEAdaptivity
using DelimitedFiles
using Plots
using Printf

gr()
default(;
    guidefontsize = 10,
    tickfontsize = 8,
    legendfontsize = 8,
    titlefontsize = 10,
    linewidth = 1.4,
)

const ROOT = normpath(joinpath(@__DIR__, ".."))
const FIG = joinpath(ROOT, "figures")
const DATA = joinpath(ROOT, "data")
mkpath(FIG)
mkpath(DATA)

function mesh_panel(xyz, conn, nodal;
        title="", vmin=minimum(nodal), vmax=maximum(nodal), cmap=:viridis)
    p = plot(;
        aspect_ratio = :equal,
        legend = false,
        xlims = (-0.02, 1.02),
        ylims = (-0.02, 1.02),
        title = title,
        xlabel = "x",
        ylabel = "y",
        size = (420, 400),
        grid = false,
        framestyle = :box,
    )
    span = max(vmax - vmin, 1e-12)
    colors = cgrad(cmap)
    for e in axes(conn, 1)
        ids = conn[e, :]
        value = sum(nodal[ids]) / length(ids)
        t = clamp((value - vmin) / span, 0, 1)
        plot!(p, Shape(xyz[ids, 1], xyz[ids, 2]);
            fillcolor = colors[t],
            linecolor = :black,
            linewidth = 0.35,
            colorbar = false,
            label = "",
        )
    end
    scatter!(p, [NaN, NaN], [NaN, NaN];
        marker_z = [vmin, vmax],
        seriescolor = colors,
        colorbar = true,
        markersize = 0,
        label = "",
    )
    return p
end

function geometry_panel(xyz, conn; title="")
    p = plot(;
        aspect_ratio = :equal,
        legend = false,
        xlims = (-0.02, 1.02),
        ylims = (-0.02, 1.02),
        title = title,
        xlabel = "x",
        ylabel = "y",
        size = (420, 400),
        grid = false,
        framestyle = :box,
    )
    for e in axes(conn, 1)
        ids = [conn[e, :]; conn[e, 1]]
        plot!(p, xyz[ids, 1], xyz[ids, 2]; color = :black, linewidth = 0.45, label = "")
    end
    return p
end

function metric_title(name, id, functional, values)
    return @sprintf("%s\nmetric %d    F = %.3e    max = %.3e", name, id, functional, maximum(values))
end

"""Common color limits of a TMOP / GCEE pair."""
function shared_limits(left, right)
    lo = min(minimum(left), minimum(right))
    hi = max(maximum(left), maximum(right))
    if hi <= lo
        hi = lo + 1e-14
    end
    return lo, hi
end

"""Element color is the quadrature mean of the TMOP metric on that element."""
function metric_panel(xyz, conn, values; title="", cmap=:YlOrRd, vmin=nothing, vmax=nothing)
    lo = vmin === nothing ? minimum(values) : vmin
    hi = vmax === nothing ? maximum(values) : vmax
    if hi <= lo
        hi = lo + 1e-14
    end
    p = plot(;
        aspect_ratio = :equal,
        legend = false,
        xlims = (-0.02, 1.02),
        ylims = (-0.02, 1.02),
        title = title,
        xlabel = "x",
        ylabel = "y",
        size = (460, 430),
        grid = false,
        framestyle = :box,
        titlefontsize = 9,
    )
    span = hi - lo
    colors = cgrad(cmap)
    for e in axes(conn, 1)
        ids = conn[e, :]
        t = clamp((values[e] - lo) / span, 0, 1)
        plot!(p, Shape(xyz[ids, 1], xyz[ids, 2]);
            fillcolor = colors[t],
            linecolor = :black,
            linewidth = 0.35,
            colorbar = false,
            label = "",
        )
    end
    scatter!(p, [NaN, NaN], [NaN, NaN];
        marker_z = [lo, hi],
        seriescolor = colors,
        colorbar = true,
        colorbar_title = "mu",
        markersize = 0,
        label = "",
    )
    return p
end

"""Uniform mesh plus one relocation. `driver` is `:tmop` or `:gcee`."""
function poisson_pair(n; sharpness=20.0, iterations=8, order=3, driver::Symbol=:tmop)
    patch = unit_square_mesh(n; order=order)
    α = Float64(sharpness)
    solution(x, y) = circular_front(x, y; sharpness=α)
    if driver == :tmop
        grad(x, y) = circular_front_grad(x, y; sharpness=α)
        target = field_target(gradient_shape_size(grad, patch, patch.xyz0))
    elseif driver == :gcee
        target = gcee_target(patch.xyz0, patch, solution)
    else
        error("unknown driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, 9, target; iterations=iterations)
    mu = element_metrics(xyz, patch, 9, target)
    uh0 = solve_poisson(patch.xyz0, patch, solution)
    uh1 = solve_poisson(xyz, patch, solution)
    l2u, eu = poisson_errors(patch.xyz0, patch, solution, uh0)
    l2a, ea = poisson_errors(xyz, patch, solution, uh1)
    exact = [solution(patch.xyz0[i, 1], patch.xyz0[i, 2]) for i in axes(patch.xyz0, 1)]
    return (;
        patch, xyz, uh0, uh1, exact, solution,
        l2u, eu, l2a, ea, history, mu, n, driver,
        ndofs = (n - 1)^2,
        h = 1 / n,
    )
end

"""Oblique shock: uniform mesh and one relocation. `driver` is `:analytic` or `:gcee`."""
function shock_pair(n; slope=20.0, iterations=8, order=3, driver::Symbol=:analytic)
    patch = unit_square_mesh(n; order=order)
    solution(x, y) = oblique_shock(x, y; slope=Float64(slope))
    if driver == :analytic
        target = field_target(orientation_target(n))
    elseif driver == :gcee
        target = gcee_target(patch.xyz0, patch, solution)
    else
        error("unknown driver $driver")
    end
    xyz, history = relocate(patch.xyz0, patch, 107, target; iterations=iterations)
    mu = element_metrics(xyz, patch, 107, target)
    uh0 = solve_poisson(patch.xyz0, patch, solution)
    uh1 = solve_poisson(xyz, patch, solution)
    l2u, eu = poisson_errors(patch.xyz0, patch, solution, uh0)
    l2a, ea = poisson_errors(xyz, patch, solution, uh1)
    exact = [solution(patch.xyz0[i, 1], patch.xyz0[i, 2]) for i in axes(patch.xyz0, 1)]
    return (;
        patch, xyz, uh0, uh1, exact,
        l2u, eu, l2a, ea, history, mu, n, driver,
        ndofs = (n - 1)^2,
        h = 1 / n,
    )
end

"""One TMOP relocation and one GCEE relocation, same metric, same feature field."""
function pair_panels(metric, tmop_target, field, title_tmop, title_gcee; n=12, iterations=8, order=3)
    patch = unit_square_mesh(n; order=order)
    gtarget = gcee_target(patch.xyz0, patch, field)
    xyz_t, hist_t = relocate(patch.xyz0, patch, metric, tmop_target; iterations=iterations)
    xyz_g, hist_g = relocate(patch.xyz0, patch, metric, gtarget; iterations=iterations)
    area_t = minimum_area(xyz_t, patch)
    area_g = minimum_area(xyz_g, patch)
    mu_t = element_metrics(xyz_t, patch, metric, tmop_target)
    mu_g = element_metrics(xyz_g, patch, metric, gtarget)
    lo, hi = shared_limits(mu_t, mu_g)
    @printf("%-22s  F TMOP %.4e -> %.4e  max mu %.3e | F GCEE %.4e -> %.4e  max mu %.3e\n",
        title_tmop, hist_t[1], hist_t[end], maximum(mu_t), hist_g[1], hist_g[end], maximum(mu_g))
    return (
        geometry_panel(xyz_t, patch.conn; title=title_tmop),
        geometry_panel(xyz_g, patch.conn; title=title_gcee),
        hist_t[1], hist_t[end], hist_g[1], hist_g[end], area_t, area_g, metric,
        metric_panel(xyz_t, patch.conn, mu_t;
            title = metric_title(title_tmop, metric, hist_t[end], mu_t), vmin=lo, vmax=hi),
        metric_panel(xyz_g, patch.conn, mu_g;
            title = metric_title(title_gcee, metric, hist_g[end], mu_g), vmin=lo, vmax=hi),
        maximum(mu_t), maximum(mu_g),
    )
end

function savefig_pair(path, plotobj)
    savefig(plotobj, path * ".pdf")
    savefig(plotobj, path * ".png")
    return path
end

function observed_rate(h, err)
    length(h) >= 2 || return NaN
    return log(err[end] / err[end - 1]) / log(h[end] / h[end - 1])
end

function save_convergence(path, runs, adapted; stem)
    h = [r.h for r in runs]
    l2u = [r.l2u for r in runs]
    l2a = [r.l2a for r in runs]
    l2g = [g.l2a for g in adapted]
    eu = [r.eu for r in runs]
    ea = [r.ea for r in runs]
    eg = [g.ea for g in adapted]
    @printf("%s rates: L2 u %.2f, TMOP %.2f, GCEE %.2f; energia u %.2f, TMOP %.2f, GCEE %.2f\n",
        stem,
        observed_rate(h, l2u), observed_rate(h, l2a), observed_rate(h, l2g),
        observed_rate(h, eu), observed_rate(h, ea), observed_rate(h, eg))
    guide_l2 = l2u[1] .* (h ./ h[1]) .^ 2
    guide_e = eu[1] .* (h ./ h[1])
    conv = plot(;
        layout = (1, 2),
        size = (860, 400),
        xscale = :log10,
        yscale = :log10,
        xlabel = "h",
        legend = :bottomleft,
    )
    plot!(conv, h, l2u; subplot = 1, marker = :circle, label = "uniform", color = :black)
    plot!(conv, h, l2a; subplot = 1, marker = :square, label = "TMOP", color = :firebrick)
    plot!(conv, h, l2g; subplot = 1, marker = :diamond, label = "GCEE", color = :steelblue)
    plot!(conv, h, guide_l2; subplot = 1, linestyle = :dash, color = :gray, label = "O(h^2)")
    plot!(conv, h, eu; subplot = 2, marker = :circle, label = "uniform", color = :black)
    plot!(conv, h, ea; subplot = 2, marker = :square, label = "TMOP", color = :firebrick)
    plot!(conv, h, eg; subplot = 2, marker = :diamond, label = "GCEE", color = :steelblue)
    plot!(conv, h, guide_e; subplot = 2, linestyle = :dash, color = :gray, label = "O(h)")
    plot!(conv; subplot = 1, ylabel = "L2 error", title = "L2 norm")
    plot!(conv; subplot = 2, ylabel = "H1 seminorm", title = "energy norm")
    savefig_pair(joinpath(FIG, path), conv)
    table = hcat([r.n for r in runs], [r.ndofs for r in runs], h, l2u, l2a, l2g, eu, ea, eg)
    header = ["n", "interior_dofs", "h", "L2_uniform", "L2_TMOP", "L2_GCEE", "energy_uniform", "energy_TMOP", "energy_GCEE"]
    open(joinpath(DATA, stem * ".csv"), "w") do io
        println(io, join(header, ","))
        writedlm(io, table, ',')
    end
    return nothing
end

function solution_figure(path, uniform_run, adapted_run)
    vmin = min(
        minimum(uniform_run.exact), minimum(uniform_run.uh0),
        minimum(uniform_run.uh1), minimum(adapted_run.uh1),
    )
    vmax = max(
        maximum(uniform_run.exact), maximum(uniform_run.uh0),
        maximum(uniform_run.uh1), maximum(adapted_run.uh1),
    )
    conn = uniform_run.patch.conn
    n = uniform_run.n
    panels = [
        mesh_panel(uniform_run.patch.xyz0, conn, uniform_run.exact;
            title = "exact, $(n)x$(n)", vmin=vmin, vmax=vmax),
        mesh_panel(uniform_run.patch.xyz0, conn, uniform_run.uh0;
            title = "uniform", vmin=vmin, vmax=vmax),
        mesh_panel(uniform_run.xyz, conn, uniform_run.uh1;
            title = "TMOP", vmin=vmin, vmax=vmax),
        mesh_panel(adapted_run.xyz, conn, adapted_run.uh1;
            title = "GCEE", vmin=vmin, vmax=vmax),
    ]
    savefig_pair(joinpath(FIG, path), plot(panels...; layout = (2, 2), size = (860, 780)))
    return nothing
end

function metric_pair_figure(path, tmop_run, gcee_run, id)
    lo, hi = shared_limits(tmop_run.mu, gcee_run.mu)
    panels = [
        metric_panel(tmop_run.xyz, tmop_run.patch.conn, tmop_run.mu;
            title = metric_title("TMOP", id, tmop_run.history[end], tmop_run.mu), vmin=lo, vmax=hi),
        metric_panel(gcee_run.xyz, gcee_run.patch.conn, gcee_run.mu;
            title = metric_title("GCEE", id, gcee_run.history[end], gcee_run.mu), vmin=lo, vmax=hi),
    ]
    savefig_pair(joinpath(FIG, path), plot(panels...; layout = (1, 2), size = (920, 460)))
    return nothing
end

function report_pde(label, runs, adapted)
    println(label)
    for (r, g) in zip(runs, adapted)
        @printf("n=%2d  dofs=%4d  L2 u %.4e  TMOP %.4e  GCEE %.4e  E u %.4e  TMOP %.4e  GCEE %.4e\n",
            r.n, r.ndofs, r.l2u, r.l2a, g.l2a, r.eu, r.ea, g.ea)
        @printf("         F TMOP %.4e -> %.4e    F GCEE %.4e -> %.4e\n",
            r.history[1], r.history[end], g.history[1], g.history[end])
    end
    return nothing
end

function main()
    levels = [4, 6, 8, 12]
    iterations = 8
    ngeo = 12

    runs = poisson_pair.(levels; sharpness=20.0, iterations=iterations, order=3, driver=:tmop)
    gcee = poisson_pair.(levels; sharpness=20.0, iterations=iterations, order=3, driver=:gcee)
    report_pde("Poisson front, sharpness 20", runs, gcee)
    solution_figure("poisson_solution", runs[end], gcee[end])
    metric_pair_figure("poisson_metric", runs[end], gcee[end], 9)
    save_convergence("convergence_poisson", runs, gcee; stem="convergence_poisson")

    sharp = poisson_pair.(levels; sharpness=200.0, iterations=iterations, order=3, driver=:tmop)
    sharp_g = poisson_pair.(levels; sharpness=200.0, iterations=iterations, order=3, driver=:gcee)
    report_pde("Poisson front, sharpness 200", sharp, sharp_g)
    solution_figure("poisson_solution_sharp", sharp[end], sharp_g[end])
    metric_pair_figure("poisson_sharp_metric", sharp[end], sharp_g[end], 9)
    save_convergence("convergence_poisson_sharp", sharp, sharp_g; stem="convergence_poisson_sharp")

    shocks = shock_pair.(levels; iterations=iterations, order=3, driver=:analytic)
    shocks_g = shock_pair.(levels; iterations=iterations, order=3, driver=:gcee)
    report_pde("Oblique shock, slope 20, metric 107", shocks, shocks_g)
    solution_figure("shock_solution", shocks[end], shocks_g[end])
    metric_pair_figure("shock_metric", shocks[end], shocks_g[end], 107)
    save_convergence("convergence_shock", shocks, shocks_g; stem="convergence_shock")

    ring_field(x, y) = ring_feature(x, y)
    tid9_field(x, y) = ring_feature(x, y; r1=0.15, r2=0.35)
    band_patch = unit_square_mesh(ngeo; order=3)
    band_size, metric_size = sinusoidal_target(band_patch, band_patch.xyz0; variant=:size)
    band_shape, metric_shape = sinusoidal_target(band_patch, band_patch.xyz0; variant=:shape)
    band_both, metric_both = sinusoidal_target(band_patch, band_patch.xyz0; variant=:both)

    rows = [
        pair_panels(7, field_target(ring_size_target()), ring_field,
            "annulus TMOP", "annulus GCEE"; n=ngeo, iterations=iterations),
        pair_panels(7, field_target(ring_size_target(; r1=0.15, r2=0.35)), tid9_field,
            "tid 9 TMOP", "tid 9 GCEE"; n=ngeo, iterations=iterations),
        pair_panels(9, field_target(ring_aspect_target()), ring_field,
            "tid 10 TMOP", "tid 10 GCEE"; n=ngeo, iterations=iterations),
        pair_panels(metric_size, band_size, sinusoidal_eta,
            "band size TMOP", "band size GCEE"; n=ngeo, iterations=iterations),
        pair_panels(metric_shape, band_shape, sinusoidal_eta,
            "band shape TMOP", "band shape GCEE"; n=ngeo, iterations=iterations),
        pair_panels(metric_both, band_both, sinusoidal_eta,
            "band both TMOP", "band both GCEE"; n=ngeo, iterations=iterations),
        pair_panels(7, field_target(material_interface_target(band_patch, band_patch.xyz0)), material_feature,
            "interface TMOP", "interface GCEE"; n=ngeo, iterations=iterations),
        pair_panels(94, field_target(discrete_size_target(band_patch, band_patch.xyz0)), size_feature,
            "discrete size TMOP", "discrete size GCEE"; n=ngeo, iterations=iterations),
    ]
    savefig_pair(joinpath(FIG, "meshes_tmop"), plot(
        rows[1][1], rows[1][2], rows[4][1], rows[4][2];
        layout = (2, 2), size = (860, 780),
    ))
    savefig_pair(joinpath(FIG, "meshes_ring"), plot(
        rows[2][1], rows[2][2], rows[3][1], rows[3][2];
        layout = (2, 2), size = (860, 780),
    ))
    savefig_pair(joinpath(FIG, "meshes_band"), plot(
        rows[5][1], rows[5][2], rows[6][1], rows[6][2];
        layout = (2, 2), size = (860, 780),
    ))
    savefig_pair(joinpath(FIG, "meshes_indicators"), plot(
        rows[7][1], rows[7][2], rows[8][1], rows[8][2];
        layout = (2, 2), size = (860, 780),
    ))
    savefig_pair(joinpath(FIG, "meshes_tmop_metric"), plot(
        rows[1][10], rows[1][11], rows[4][10], rows[4][11];
        layout = (2, 2), size = (920, 860),
    ))
    savefig_pair(joinpath(FIG, "meshes_ring_metric"), plot(
        rows[2][10], rows[2][11], rows[3][10], rows[3][11];
        layout = (2, 2), size = (920, 860),
    ))
    savefig_pair(joinpath(FIG, "meshes_band_metric"), plot(
        rows[5][10], rows[5][11], rows[6][10], rows[6][11];
        layout = (2, 2), size = (920, 860),
    ))
    savefig_pair(joinpath(FIG, "meshes_indicators_metric"), plot(
        rows[7][10], rows[7][11], rows[8][10], rows[8][11];
        layout = (2, 2), size = (920, 860),
    ))

    names = ["ring-paper", "ring-tid9", "ring-tid10", "band-size", "band-shape", "band-both", "material", "discrete-size"]
    open(joinpath(DATA, "cases_gcee.csv"), "w") do io
        println(io, "case,metric,F0_TMOP,F1_TMOP,F0_GCEE,F1_GCEE,min_area_TMOP,min_area_GCEE,max_mu_TMOP,max_mu_GCEE")
        for (name, row) in zip(names, rows)
            println(io, join((name, row[9], row[3], row[4], row[5], row[6], row[7], row[8], row[12], row[13]), ","))
        end
    end
    println("scritto ", FIG)
    return runs
end

main()
