# Writes every table and figure that Paper A includes from `generated/`.
# The recovered-indicator study replaces the Python prototype. The shortlist
# figures use the same metric, solver and mesh as `scripts/plot_benchmarks.jl`,
# with the Gram target in place of the published target.
#
#   julia --project=. scripts/export_paper.jl --paper /path/to/PaperA
#   julia --project=. scripts/export_paper.jl --paper /path/to/PaperA --quick
#
# `--quick` drops the n = 16 recovery level and the finer shortlist meshes.
# The paper tables need the full run.

using GCEEAdaptivity
using DelimitedFiles
using LinearAlgebra
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
const DATA = joinpath(ROOT, "data")
const FIG = joinpath(ROOT, "figures")
mkpath(DATA); mkpath(FIG)

# Paper A colour convention.
const C_GCEE = RGB(0.659, 0.216, 0.173)      # #a8372c
const C_TMOPU = RGB(0.141, 0.384, 0.545)     # #24628b
const C_TMOPH = RGB(0.220, 0.525, 0.357)     # #38865b
const C_UNIFORM = RGB(0.267, 0.267, 0.267)  # #444444
const MODE_COLOR = Dict("gcee" => C_GCEE, "tmop_uniform" => C_TMOPU,
                        "tmop_hessian" => C_TMOPH, "tmop" => C_TMOPU)
const MODE_LABEL = Dict("gcee" => "GCEE", "tmop_uniform" => "TMOP, uniform target",
                        "tmop_hessian" => "TMOP, Hessian shape", "tmop" => "TMOP, Hessian shape")
const SHORT = Dict("gcee" => "GCEE", "tmop_uniform" => "TMOP-U",
                   "tmop_hessian" => "TMOP-H", "tmop" => "TMOP")
const MESH_NAME = Dict("regular" => "regular", "checkerboard" => "alternating",
                       "random" => "random", "curved_A" => "curved A",
                       "curved_B" => "curved B", "curved_C" => "curved C")

"""Write a complete `tabular`, so the paper can `\\input` it inside a table float."""
function write_table(path, headers, rows, alignment)
    body = ["\\begin{tabular}{$alignment}", "\\toprule",
            join(headers, " & ") * "\\\\", "\\midrule"]
    append!(body, [join(r, " & ") * "\\\\" for r in rows])
    append!(body, ["\\bottomrule", "\\end{tabular}"])
    write(path, join(body, "\n") * "\n")
    return path
end

"""Outline of a quadrilateral mesh, in the style of the benchmark figures."""
function mesh_outline(xyz, conn; title="", lims=nothing, color=:black)
    p = plot(; aspect_ratio = :equal, legend = false, title = title,
             xlabel = "x", ylabel = "y", size = (420, 400), grid = false,
             framestyle = :box)
    if lims !== nothing
        plot!(p; xlims = lims[1], ylims = lims[2])
    end
    for e in axes(conn, 1)
        ids = [conn[e, 1:4]; conn[e, 1]]
        plot!(p, xyz[ids, 1], xyz[ids, 2]; color = color, linewidth = 0.45, label = "")
    end
    return p
end

"""Element means of a nodal field, the values that colour the cells."""
function cell_means(conn, nodal)
    return [sum(nodal[conn[e, 1:4]]) / 4 for e in axes(conn, 1)]
end

"""Two-line title: the panel name, then the plotted quantity and its maximum."""
function field_heading(name, values)
    return @sprintf("%s\nu    max = %.3e", name, maximum(values))
end

"""Filled mesh coloured by an element value, with a colour bar."""
function mesh_field(xyz, conn, nodal; title="", lims=nothing, cmap=:YlOrRd,
        vmin=minimum(nodal), vmax=maximum(nodal), linewidth=0.35, bar="u")
    p = plot(; aspect_ratio = :equal, legend = false, title = title,
             xlabel = "x", ylabel = "y", size = (420, 400), grid = false,
             framestyle = :box, titlefontsize = 9)
    if lims !== nothing
        plot!(p; xlims = lims[1], ylims = lims[2])
    end
    span = max(vmax - vmin, 1e-12)
    colors = cgrad(cmap)
    for e in axes(conn, 1)
        ids = conn[e, 1:4]
        value = sum(nodal[ids]) / 4
        t = clamp((value - vmin) / span, 0, 1)
        plot!(p, Shape(xyz[ids, 1], xyz[ids, 2]); fillcolor = colors[t],
              linecolor = :black, linewidth = linewidth, colorbar = false, label = "")
    end
    scatter!(p, [NaN, NaN], [NaN, NaN]; marker_z = [vmin, vmax],
             seriescolor = colors, colorbar = true, colorbar_title = bar,
             markersize = 0, label = "")
    return p
end

savepair(path, p) = (savefig(p, path * ".pdf"); savefig(p, path * ".png"); path)

"""Physical image of the four reference edges, so a curved side stays curved."""
function cell_edges(kind, nodes; samples = 9)
    t = collect(range(-1.0, 1.0; length = samples))
    refs = (
        (t, fill(-1.0, samples)),
        (fill(1.0, samples), t),
        (reverse(t), fill(1.0, samples)),
        (fill(-1.0, samples), reverse(t)),
    )
    edges = Matrix{Float64}[]
    for (ξ, η) in refs
        N, _, _ = shapes(kind, [ξ η])
        push!(edges, N * nodes)
    end
    return edges
end

"""Outline of one structured mesh, corners and midside nodes included."""
function family_outline(m; title = "")
    p = plot(; aspect_ratio = :equal, legend = false, title = title,
             xlabel = "x", ylabel = "y", size = (420, 400), grid = false,
             framestyle = :box, xlims = (-0.04, 1.04), ylims = (-0.04, 1.04))
    for e in 1:ncells(m)
        nodes = m.points[m.cells[e, :], :]
        for edge in cell_edges(m.kind, nodes)
            plot!(p, edge[:, 1], edge[:, 2]; color = :black, linewidth = 0.45, label = "")
        end
    end
    return p
end

"""Element colour is the quadrature mean of the TMOP metric."""
function metric_fill(xyz, conn, values; title = "", lo = minimum(values), hi = maximum(values))
    hi <= lo && (hi = lo + 1e-14)
    p = plot(; aspect_ratio = :equal, legend = false, title = title,
             xlabel = "x", ylabel = "y", size = (460, 430), grid = false,
             framestyle = :box, titlefontsize = 9,
             xlims = (-0.02, 1.02), ylims = (-0.02, 1.02))
    span = hi - lo
    colors = cgrad(:YlOrRd)
    for e in axes(conn, 1)
        ids = conn[e, :]
        t = clamp((values[e] - lo) / span, 0, 1)
        plot!(p, Shape(xyz[ids, 1], xyz[ids, 2]); fillcolor = colors[t],
              linecolor = :black, linewidth = 0.35, colorbar = false, label = "")
    end
    scatter!(p, [NaN, NaN], [NaN, NaN]; marker_z = [lo, hi], seriescolor = colors,
             colorbar = true, colorbar_title = "mu", markersize = 0, label = "")
    return p
end

function metric_heading(name, id, functional, values)
    return @sprintf("%s\nmetric %d    F = %.3e    max = %.3e",
                    name, id, functional, maximum(values))
end

"""One published target and the Gram target on the same metric."""
function target_pair(metric, tmop_target, field; n = 12, iterations = 8, order = 3)
    patch = unit_square_mesh(n; order = order)
    gtarget = gcee_target(patch.xyz0, patch, field)
    xyz_t, hist_t = relocate(patch.xyz0, patch, metric, tmop_target; iterations = iterations)
    xyz_g, hist_g = relocate(patch.xyz0, patch, metric, gtarget; iterations = iterations)
    mu_t = element_metrics(xyz_t, patch, metric, tmop_target)
    mu_g = element_metrics(xyz_g, patch, metric, gtarget)
    lo = min(minimum(mu_t), minimum(mu_g))
    hi = max(maximum(mu_t), maximum(mu_g))
    return (;
        patch, xyz_t, xyz_g, mu_t, mu_g, lo, hi, metric,
        f0t = hist_t[1], f1t = hist_t[end], f0g = hist_g[1], f1g = hist_g[end],
        area_t = minimum_area(xyz_t, patch), area_g = minimum_area(xyz_g, patch),
    )
end

function metric_pair_plot(pair, title_t, title_g)
    return (
        metric_fill(pair.xyz_t, pair.patch.conn, pair.mu_t;
            title = metric_heading(title_t, pair.metric, pair.f1t, pair.mu_t),
            lo = pair.lo, hi = pair.hi),
        metric_fill(pair.xyz_g, pair.patch.conn, pair.mu_g;
            title = metric_heading(title_g, pair.metric, pair.f1g, pair.mu_g),
            lo = pair.lo, hi = pair.hi),
    )
end

# Published size ratio admitted by each geometric target, in table order.
const QUALITY_ROWS = (
    ("ring-paper", "Annulus, hr \\S4.2", "\$10\\!:\\!1\$"),
    ("ring-tid9", "Annulus, \\texttt{tid 9}", "\$10\\!:\\!1\$"),
    ("ring-tid10", "Annulus, \\texttt{tid 10}", "\$10\\!:\\!1\$"),
    ("band-size", "Band \$W_1\$, IMR \\S6.1", "\$7\\!:\\!1\$"),
    ("band-shape", "Band \$W_2\$, IMR \\S6.1", "\$3\\!:\\!1\$"),
    ("band-both", "Band \$W_3\$, IMR \\S6.1", "\$3\\!:\\!1\$"),
    ("material", "Interface, \\texttt{tid 6}", "\$9\\!:\\!1\$"),
    ("discrete-size", "Size, \\texttt{tid 5}", "\$9\\!:\\!1\$"),
)

fmt_mu(v) = v >= 100 ? @sprintf("%.0f", v) : @sprintf("%.2f", v)

"""Convergence, quality and the three metric figures of the shortlist section."""
function export_shortlist(out; levels = (4, 6, 8, 12), ngeo = 12, iterations = 8, order = 3)
    println("shortlist ...")
    fronts = []
    for α in (20.0, 200.0)
        tmop = [run_poisson_front(; sharpness = α, n = n, iterations = iterations,
                                  order = order, driver = :gradient) for n in levels]
        gcee = [run_poisson_front(; sharpness = α, n = n, iterations = iterations,
                                  order = order, driver = :gcee) for n in levels]
        push!(fronts, (tmop, gcee))
    end
    shock_t = [run_oblique_shock(; n = n, iterations = iterations, order = order,
                                 driver = :analytic) for n in levels]
    shock_g = [run_oblique_shock(; n = n, iterations = iterations, order = order,
                                 driver = :gcee) for n in levels]
    sweeps = (
        ("Circular front, \$\\beta=20\$", "circular front, beta = 20",
         "convergence_poisson", fronts[1]...),
        ("Circular front, \$\\beta=200\$", "circular front, beta = 200",
         "convergence_poisson_sharp", fronts[2]...),
        ("Oblique shock, \$\\mu_{107}\$", "oblique shock",
         "convergence_shock", shock_t, shock_g),
    )

    lines = [raw"\begin{tabular}{llrrrrrr}", raw"\toprule",
             raw"& & \multicolumn{3}{c}{$\|u-u_h\|_{L^2}$} & \multicolumn{3}{c}{$|u-u_h|_{H^1}$}" * "\\\\",
             raw"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
             raw"Case & $n$ & uniform & TMOP & GCEE & uniform & TMOP & GCEE" * "\\\\",
             raw"\midrule"]
    for (k, (label, _, _, tmop, gcee)) in enumerate(sweeps)
        for (j, (r, g)) in enumerate(zip(tmop, gcee))
            head = j == 1 ? label : ""
            push!(lines, @sprintf("%s & %d & %.4f & %.4f & %.4f & %.3f & %.3f & %.3f\\\\",
                                  head, r.n, r.l2_uniform, r.l2_adapted, g.l2_adapted,
                                  r.energy_uniform, r.energy_adapted, g.energy_adapted))
        end
        k < length(sweeps) && push!(lines, raw"\midrule")
    end
    append!(lines, [raw"\bottomrule", raw"\end{tabular}"])
    write(joinpath(out, "shortlist_table.tex"), join(lines, "\n") * "\n")

    for (_, _, stem, tmop, gcee) in sweeps
        open(joinpath(DATA, stem * ".csv"), "w") do io
            println(io, "n,interior_dofs,h,L2_uniform,L2_TMOP,L2_GCEE,energy_uniform,energy_TMOP,energy_GCEE")
            for (r, g) in zip(tmop, gcee)
                @printf(io, "%d,%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                        r.n, (r.n - 1)^2, 1 / r.n,
                        r.l2_uniform, r.l2_adapted, g.l2_adapted,
                        r.energy_uniform, r.energy_adapted, g.energy_adapted)
            end
        end
    end

    conv = plot(; layout = (2, 3), size = (940, 520), legend = false)
    for (col, (_, title, _, tmop, gcee)) in enumerate(sweeps)
        h = [1 / r.n for r in tmop]
        blocks = (
            ([r.l2_uniform for r in tmop], [r.l2_adapted for r in tmop],
             [g.l2_adapted for g in gcee], 2, "L2 error"),
            ([r.energy_uniform for r in tmop], [r.energy_adapted for r in tmop],
             [g.energy_adapted for g in gcee], 1, "H1 seminorm"),
        )
        for (row, (yu, yt, yg, power, ylab)) in enumerate(blocks)
            sp = (row - 1) * 3 + col
            show = sp == 1
            plot!(conv, h, yu; subplot = sp, marker = :circle, markersize = 3,
                  color = C_UNIFORM, label = show ? "uniform" : "",
                  xscale = :log10, yscale = :log10, xlabel = "h",
                  ylabel = col == 1 ? ylab : "", title = row == 1 ? title : "")
            plot!(conv, h, yt; subplot = sp, marker = :square, markersize = 3,
                  color = C_TMOPU, label = show ? "TMOP target" : "")
            plot!(conv, h, yg; subplot = sp, marker = :diamond, markersize = 3,
                  color = C_GCEE, label = show ? "GCEE target" : "")
            guide = yu[1] .* (h ./ h[1]) .^ power
            plot!(conv, h, guide; subplot = sp, linestyle = :dot, color = :gray,
                  label = show ? (power == 2 ? "O(h^2)" : "O(h)") : "")
        end
    end
    plot!(conv; subplot = 1, legend = :topright)
    savepair(joinpath(out, "shortlist_convergence"), conv)

    ring(x, y) = ring_feature(x, y)
    tid9(x, y) = ring_feature(x, y; r1 = 0.15, r2 = 0.35)
    band = unit_square_mesh(ngeo; order = order)
    band_size, metric_size = sinusoidal_target(band, band.xyz0; variant = :size)
    band_shape, metric_shape = sinusoidal_target(band, band.xyz0; variant = :shape)
    band_both, metric_both = sinusoidal_target(band, band.xyz0; variant = :both)
    cases = [
        ("ring-paper", target_pair(7, field_target(ring_size_target()), ring;
                                    n = ngeo, iterations = iterations, order = order)),
        ("ring-tid9", target_pair(7, field_target(ring_size_target(; r1 = 0.15, r2 = 0.35)), tid9;
                                  n = ngeo, iterations = iterations, order = order)),
        ("ring-tid10", target_pair(9, field_target(ring_aspect_target()), ring;
                                   n = ngeo, iterations = iterations, order = order)),
        ("band-size", target_pair(metric_size, band_size, sinusoidal_eta;
                                  n = ngeo, iterations = iterations, order = order)),
        ("band-shape", target_pair(metric_shape, band_shape, sinusoidal_eta;
                                   n = ngeo, iterations = iterations, order = order)),
        ("band-both", target_pair(metric_both, band_both, sinusoidal_eta;
                                  n = ngeo, iterations = iterations, order = order)),
        ("material", target_pair(7, field_target(material_interface_target(band, band.xyz0)),
                                 material_feature; n = ngeo, iterations = iterations, order = order)),
        ("discrete-size", target_pair(94, field_target(discrete_size_target(band, band.xyz0)),
                                      size_feature; n = ngeo, iterations = iterations, order = order)),
    ]
    byname = Dict(name => pair for (name, pair) in cases)

    shock = target_pair(107, field_target(orientation_target(ngeo)),
                        (x, y) -> oblique_shock(x, y; slope = 20.0);
                        n = ngeo, iterations = iterations, order = order)
    left, right = metric_pair_plot(shock, "TMOP", "GCEE")
    savepair(joinpath(out, "shortlist_shock_metric"),
             plot(left, right; layout = (1, 2), size = (920, 460)))

    r9 = metric_pair_plot(byname["ring-tid9"], "tid 9 TMOP", "tid 9 GCEE")
    r10 = metric_pair_plot(byname["ring-tid10"], "tid 10 TMOP", "tid 10 GCEE")
    savepair(joinpath(out, "shortlist_ring_metric"),
             plot(r9..., r10...; layout = (2, 2), size = (920, 860)))

    iface = metric_pair_plot(byname["material"], "interface TMOP", "interface GCEE")
    disc = metric_pair_plot(byname["discrete-size"], "discrete size TMOP", "discrete size GCEE")
    savepair(joinpath(out, "shortlist_indicators_metric"),
             plot(iface..., disc...; layout = (2, 2), size = (920, 860)))

    area0 = 1 / ngeo^2
    qlines = [raw"\begin{tabular}{llcrrrrrr}", raw"\toprule",
              raw"& & & \multicolumn{2}{c}{$\min_K |K|\ (\times 10^{-3})$}" *
              raw" & \multicolumn{2}{c}{$|K_{\mathrm{unif}}|/\min_K|K|$}" *
              raw" & \multicolumn{2}{c}{$\max_K \mu$}" * "\\\\",
              raw"\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}",
              raw"Target & $\mu$ & ratio & TMOP & GCEE & TMOP & GCEE & TMOP & GCEE" * "\\\\",
              raw"\midrule"]
    open(joinpath(DATA, "cases_gcee.csv"), "w") do io
        println(io, "case,metric,F0_TMOP,F1_TMOP,F0_GCEE,F1_GCEE,min_area_TMOP,min_area_GCEE,max_mu_TMOP,max_mu_GCEE")
        for (key, label, ratio) in QUALITY_ROWS
            pair = byname[key]
            mt, mg = maximum(pair.mu_t), maximum(pair.mu_g)
            println(io, join((key, pair.metric, pair.f0t, pair.f1t, pair.f0g, pair.f1g,
                              pair.area_t, pair.area_g, mt, mg), ","))
            push!(qlines, @sprintf("%s & \$\\mu_{%d}\$ & %s & %.2f & %.2f & %.1f & %.1f & %s & %s\\\\",
                                   label, pair.metric, ratio, 1e3 * pair.area_t, 1e3 * pair.area_g,
                                   area0 / pair.area_t, area0 / pair.area_g, fmt_mu(mt), fmt_mu(mg)))
        end
    end
    append!(qlines, [raw"\bottomrule", raw"\end{tabular}"])
    write(joinpath(out, "shortlist_quality.tex"), join(qlines, "\n") * "\n")
    return nothing
end

# ---------------------------------------------------------------------------

function main()
    paper = nothing
    quick = false
    for (i, a) in enumerate(ARGS)
        a == "--paper" && (paper = ARGS[i + 1])
        a == "--quick" && (quick = true)
    end
    control_only = any(a == "--only" && ARGS[i + 1] == "control" for (i, a) in enumerate(ARGS[1:end - 1]))
    paper === nothing && error("pass --paper /path/to/PaperA")
    out = joinpath(paper, "generated")
    mkpath(out)
    levels = quick ? (4, 8) : (4, 8, 16)

    # --- Experiments 1-3 ----------------------------------------------------
    if !control_only
    println("recovery study ...")
    rec = run_recovery_study(; levels = levels, verbose = false)
    @printf("  %d cases\n", length(rec))
    top = quick ? 8 : 16
    selected = [r for r in rec if r["n"] == top && r["seed"] == 7 &&
                (r["family"] == "regular" ||
                 (r["family"] == "checkerboard" && r["amplitude"] == 0.4) ||
                 (r["family"] == "random" && r["amplitude"] == 0.25) ||
                 (startswith(r["family"], "curved") && r["amplitude"] in (0.3, 0.5)))]
    write_table(joinpath(out, "recovery_table.tex"),
        ["Space", "Mesh", "\$e_{2,I}\$", "\$e_{2,h}\$", "\$e_{3,I}\$", "\$e_{3,h}\$",
         "\$\\delta_G/\\vartheta_o\$"],
        [[r["kind"], MESH_NAME[r["family"]],
          @sprintf("%.3f", r["exact_nodal_degree2"]), @sprintf("%.3f", r["recovery_degree2"]),
          @sprintf("%.3f", r["exact_nodal_degree3"]), @sprintf("%.3f", r["recovery_degree3"]),
          @sprintf("%.3f", r["recovery_budget"] / r["oracle"])] for r in selected],
        "llrrrrr")
    write_table(joinpath(out, "estimator_table.tex"),
        ["Space", "Mesh", "\$\\vartheta_o/E_I\$", "\$\\vartheta_r/E_I\$",
         "\$\\vartheta_r/E_h\$", "\$\\eta_{ZZ}/E_h\$", "\$\\eta_R/E_h\$"],
        [[r["kind"], MESH_NAME[r["family"]],
          @sprintf("%.3f", r["oracle_I"]), @sprintf("%.3f", r["rec_I"]),
          @sprintf("%.3f", r["rec_h"]), @sprintf("%.3f", r["ZZ_h"]),
          @sprintf("%.3f", r["residual_h"])] for r in selected],
        "llrrrrr")

    conv = plot(; layout = (1, 3), size = (900, 290), legend = :topleft)
    for (j, kind) in enumerate(("Q4", "Q8", "Q9"))
        for (family, amp, color) in (("regular", 0.0, :black),
                                     ("checkerboard", 0.4, C_TMOPU),
                                     ("random", 0.25, C_GCEE))
            rr = sort([r for r in rec if r["kind"] == kind && r["family"] == family &&
                       r["amplitude"] == amp && r["seed"] == 7], by = r -> r["n"])
            isempty(rr) && continue
            plot!(conv, [r["n"] for r in rr], [r["recovery_degree3"] for r in rr];
                  subplot = j, marker = :circle, markersize = 3, color = color,
                  label = j == 1 ? MESH_NAME[family] : "",
                  xscale = :log10, yscale = :log10, xlabel = "n",
                  ylabel = j == 1 ? "degree-3 recovery error" : "", title = kind)
        end
    end
    savepair(joinpath(out, "recovery_effectivity"), conv)

    families = (
        ("regular", 0.0, "regular"),
        ("checkerboard", 0.4, "alternating"),
        ("random", 0.25, "random"),
        ("curved_A", 0.3, "curved A"),
        ("curved_B", 0.3, "curved B"),
        ("curved_C", 0.5, "curved C"),
    )
    outlines = [family_outline(structured_mesh(6; kind = "Q8", family = f,
                                              amplitude = a, seed = 7); title = title)
                for (f, a, title) in families]
    savepair(joinpath(out, "recovery_meshes"),
             plot(outlines...; layout = (2, 3), size = (980, 640)))

    end

    # --- Experiments 4-5 ----------------------------------------------------
    println("relocation study ...")
    rel, frozen, relmesh, relcells = run_relocation_study(; verbose = false)
    # The solid GCEE curve is the relocation without the shape sub-step.
    rel_free, frozen_free, relmesh_free, _ = run_relocation_study(;
        modes = ("gcee",), ratio_growth = Inf, verbose = false)
    rel_report = vcat(rel_free, [r for r in rel if r["mode"] != "gcee"])
    frozen_report = vcat(frozen_free, [a for a in frozen if a["mode"] != "gcee"])
    rows = String[]
    table_rows = Vector{String}[]
    for r in rel_report
        r["cycle"] == 4 || continue
        f = first(a for a in frozen_report if a["n"] == r["n"] && a["field"] == r["field"] &&
                  a["mode"] == r["mode"] && a["cycle"] == 0)
        c1 = first(a for a in rel_report if a["n"] == r["n"] && a["field"] == r["field"] &&
                   a["mode"] == r["mode"] && a["cycle"] == 1)
        push!(table_rows, [string(r["n"]), r["field"], SHORT[r["mode"]],
                           @sprintf("%.3f", f["J_after"] / f["J_before"]),
                           @sprintf("%.3f", c1["E_h_ratio"]),
                           @sprintf("%.3f", r["E_h_ratio"]),
                           @sprintf("%.3f", r["E_I_ratio"])])
    end
    write_table(joinpath(out, "relocation_table.tex"),
        ["\$n\$", "Field", "Objective", "\$J_1/J_0\$", "\$E_h^1/E_h^0\$",
         "\$E_h^4/E_h^0\$", "\$E_I^4/E_I^0\$"], table_rows, "rllrrrr")

    errs = plot(; layout = (2, 2), size = (860, 560))
    panel = 0
    for n in (8, 12), fieldname in ("layer", "bump")
        panel += 1
        for mode in ("gcee", "tmop_uniform", "tmop_hessian")
            source = mode == "gcee" ? rel_free : rel
            rr = sort([r for r in source if r["n"] == n && r["field"] == fieldname &&
                       r["mode"] == mode], by = r -> r["cycle"])
            plot!(errs, [r["cycle"] for r in rr], [r["E_h_ratio"] for r in rr];
                  subplot = panel, marker = :circle, markersize = 3,
                  color = MODE_COLOR[mode], label = panel == 1 ? MODE_LABEL[mode] : "",
                  xlabel = "cycle", ylabel = "\$E_h/E_h^0\$",
                  title = "$fieldname, n = $n", legend = panel == 1 ? :topright : false)
        end
        relaxed = sort([r for r in rel if r["n"] == n && r["field"] == fieldname &&
                        r["mode"] == "gcee"], by = r -> r["cycle"])
        plot!(errs, [r["cycle"] for r in relaxed], [r["E_h_ratio"] for r in relaxed];
              subplot = panel, marker = :diamond, markersize = 3, linestyle = :dash,
              color = C_GCEE, label = panel == 1 ? "GCEE, relaxed" : "")
    end
    savepair(joinpath(out, "relocation_errors"), errs)

    # Layer at n = 12, the case named in the caption. The field panel is a fine
    # regular mesh so the colour is the feature; the other three share its scale.
    u_layer = manufactured("layer")
    fieldmesh = structured_mesh(40; kind = "Q4", family = "regular")
    square = ((-0.02, 1.02), (-0.02, 1.02))
    reloc_meshes = (
        fieldmesh.points,
        relmesh["layer_gcee_12"][1],
        relmesh["layer_tmop_uniform_12"][end],
        relmesh["layer_gcee_12"][end],
    )
    reloc_conns = (fieldmesh.cells, relcells["layer_gcee_12"],
                   relcells["layer_tmop_uniform_12"], relcells["layer_gcee_12"])
    layer_ratio = first(r["det_ratio"] for r in rel if r["n"] == 12 &&
                        r["field"] == "layer" && r["mode"] == "gcee" && r["cycle"] == 4)
    reloc_titles = ("manufactured field", "initial", "TMOP-U", "GCEE, relaxed")
    reloc_ratios = (nothing, nothing, nothing, layer_ratio)
    reloc_nodal = [[u_layer(pts[i, 1], pts[i, 2]) for i in axes(pts, 1)] for pts in reloc_meshes]
    reloc_paint = [cell_means(conn, val) for (conn, val) in zip(reloc_conns, reloc_nodal)]
    rvmin, rvmax = minimum(minimum.(reloc_paint)), maximum(maximum.(reloc_paint))
    reloc_panels = [mesh_field(pts, conn, val;
                               title = ratio === nothing ? field_heading(title, paint) :
                                       field_heading(title, paint) *
                                       @sprintf("    ratio = %.1f", ratio),
                               lims = square, cmap = :YlOrRd, vmin = rvmin, vmax = rvmax,
                               linewidth = pts === fieldmesh.points ? 0.0 : 0.35)
                    for (pts, conn, val, paint, title, ratio) in
                        zip(reloc_meshes, reloc_conns, reloc_nodal, reloc_paint,
                            reloc_titles, reloc_ratios)]
    savepair(joinpath(out, "relocation_meshes"),
             plot(reloc_panels...; layout = (2, 2), size = (920, 860)))

    # --- Displacement-only elastic control ----------------------------------
    if !control_only
    println("elastic patch ...")
    elastic_rows = Vector{String}[]
    response = Dict{String,Vector{NTuple{4,Float64}}}()
    for kind in ("Q4", "Q8", "Q9")
        response[kind] = NTuple{4,Float64}[]
        for delta in (0.0, 0.4, 0.8)
            m = structured_mesh(8; kind = kind, family = "checkerboard", amplitude = delta, seed = 7)
            sys = PatchSystem(m; order = 11)
            uh, _ = patch_solve(sys, "cubic_cos")
            s, _ = assess_patch(m, uh, "cubic_cos"; order = 13)
            g = geometric_response(m; order = 13)
            qmax = maximum(g["qmax"])
            push!(elastic_rows, [kind, @sprintf("%.1f", delta), @sprintf("%.3f", qmax),
                                 @sprintf("%.4f", s["E_I"] / s["exact_energy_norm"]),
                                 @sprintf("%.4f", s["E_h"] / s["exact_energy_norm"]),
                                 @sprintf("%.3f", s["rec_I"]), @sprintf("%.3f", s["rec_h"])])
            push!(response[kind], (delta, qmax, s["E_I"] / s["exact_energy_norm"], s["rec_h"]))
        end
    end
    write_table(joinpath(out, "elastic_patch_table.tex"),
        ["Space", "\$\\delta\$", "\$q_{\\max}\$", "\$E_I/\\|u\\|_a\$",
         "\$E_h/\\|u\\|_a\$", "\$\\vartheta_r/E_I\$", "\$\\vartheta_r/E_h\$"],
        elastic_rows, "lrrrrrr")

    resp = plot(; layout = (1, 3), size = (900, 290))
    for (j, (ylab, idx)) in enumerate((("\$q_{\\max}\$", 2),
                                       ("\$E_I/\\|u\\|_a\$", 3),
                                       ("\$\\vartheta_r/E_h\$", 4)))
        for (kind, color) in (("Q4", C_GCEE), ("Q8", C_TMOPU), ("Q9", C_TMOPH))
            d = response[kind]
            plot!(resp, [t[1] for t in d], [t[idx] for t in d]; subplot = j,
                  marker = :circle, markersize = 3, color = color,
                  label = j == 1 ? kind : "", xlabel = "\$\\delta\$", ylabel = ylab,
                  legend = j == 1 ? :topleft : false)
        end
    end
    savepair(joinpath(out, "elastic_patch_response"), resp)

    patch_panels = [mesh_outline(structured_mesh(8; kind = "Q4", family = "checkerboard",
                                                amplitude = d, seed = 7).points,
                                 structured_mesh(8; kind = "Q4", family = "checkerboard",
                                                 amplitude = d, seed = 7).cells;
                                 title = "\$\\delta = $d\$",
                                 lims = ((-0.02, 1.02), (-0.02, 1.02)))
                    for d in (0.0, 0.4, 0.8)]
    savepair(joinpath(out, "elastic_patch_meshes"),
             plot(patch_panels...; layout = (1, 3), size = (900, 300)))

    # --- Critical modes ------------------------------------------------------
    println("critical modes ...")
    kinds = ("Q4", "Q8", "Q9")
    meshes = Dict(k => structured_mesh(4; kind = k, family = "checkerboard",
                                       amplitude = 0.8, seed = 7) for k in kinds)
    modes = Dict{String,Any}()
    stats = Dict{String,Any}()
    for k in kinds
        spec = geometric_response(meshes[k]; order = 13)
        mode, s, _, _ = critical_mode(meshes[k], 6, spec["G"][6, :, :], spec["H"][6, :, :]; order = 21)
        modes[k] = mode
        stats[k] = s
    end
    crit_rows = [[src,
                  (@sprintf("%.3f", direct_response(meshes["Q4"], 6, modes[src]))),
                  (@sprintf("%.3f", direct_response(meshes["Q8"], 6, modes[src]))),
                  (@sprintf("%.3f", direct_response(meshes["Q9"], 6, modes[src]))),
                  (@sprintf("%.3f", stats[src]["cross_error2"] / stats[src]["error2"]))]
                 for src in kinds]
    write_table(joinpath(out, "critical_modes_table.tex"),
        ["Worst mode of", "Q4 response", "Q8 response", "Q9 response",
         "\$2(R_2,R_3)_a/\\|R\\|_a^2\$"], crit_rows, "lrrrr")

    crit = plot(; layout = (1, 3), size = (960, 320))
    for (j, k) in enumerate(kinds)
        m = meshes[k]
        nodes = m.points[m.cells[6, 1:4], :]
        g = range(-1, 1; length = 26)
        pts = reduce(vcat, [[a b] for a in g, b in g])
        N = shapes("Q4", pts)[1]
        xy = N * nodes
        v = mode_value(modes[k], xy)
        mag = sqrt.(v[:, 1] .^ 2 .+ v[:, 2] .^ 2)
        scatter!(crit, xy[:, 1], xy[:, 2]; subplot = j, marker_z = mag, markersize = 2.4,
                 markerstrokewidth = 0, seriescolor = cgrad(:viridis), label = "",
                 aspect_ratio = :equal, framestyle = :box, grid = false,
                 xlabel = "x", ylabel = j == 1 ? "y" : "",
                 title = @sprintf("%s, q = %.3f", k, stats[k]["q_direct"]), colorbar = j == 3)
        ids = [m.cells[6, 1:4]; m.cells[6, 1]]
        plot!(crit, m.points[ids, 1], m.points[ids, 2]; subplot = j, color = :black,
              linewidth = 0.8, label = "")
    end
    savepair(joinpath(out, "critical_modes"), crit)

    end

    # --- Experiment 6: the curved beam ---------------------------------------
    println("curved beam ...")
    part, partmesh, partcells = run_part_study(; ratio_growth = Inf, verbose = true)
    part_relax, partmesh_relax, _ = run_part_study(; modes = ("gcee",), verbose = false)
    write_table(joinpath(out, "part_table.tex"),
        ["Objective", "Cycle", "\$|u-u_h|_{H^1}\$", "GCEE rec.", "ZZ", "Residual",
         "\$\\det\$ ratio"],
        [[uppercase(r["mode"]) == "GCEE" ? "GCEE" : "TMOP", string(r["cycle"]),
          @sprintf("%.4f", r["E_h"]), @sprintf("%.4f", r["gcee_recovered"]),
          @sprintf("%.4f", r["zz"]), @sprintf("%.3f", r["residual_estimator"]),
          @sprintf("%.1f", r["det_ratio"])] for r in part],
        "llrrrrr")

    u = manufactured("beam")
    lims = ((-0.05, 2.05), (-0.05, 2.05))
    # The field panel is the exact annulus, fine enough that the cells disappear
    # and the colour is the manufactured feature that drives the relocation.
    fieldmesh = curved_beam(40; kind = "Q4", amplitude = 0.0)
    meshes = (
        fieldmesh.points,
        partmesh["gcee"][1],
        partmesh["tmop"][end],
        partmesh_relax["gcee"][end],
    )
    conns = (fieldmesh.cells, partcells, partcells, partcells)
    titles = ("manufactured field", "initial", "TMOP", "GCEE, relaxed")
    nodal = [[u(pts[i, 1], pts[i, 2]) for i in axes(pts, 1)] for pts in meshes]
    painted = [cell_means(conn, val) for (conn, val) in zip(conns, nodal)]
    vmin, vmax = minimum(minimum.(painted)), maximum(maximum.(painted))
    beam_relaxed_ratio = last([r["det_ratio"] for r in part_relax if r["mode"] == "gcee"])
    beam_ratios = (nothing, nothing, nothing, beam_relaxed_ratio)
    panels = [mesh_field(pts, conn, val;
                         title = ratio === nothing ? field_heading(title, paint) :
                                 field_heading(title, paint) *
                                 @sprintf("    ratio = %.1f", ratio),
                         lims = lims, cmap = :YlOrRd, vmin = vmin, vmax = vmax,
                         linewidth = pts === fieldmesh.points ? 0.0 : 0.35)
              for (pts, conn, val, paint, title, ratio) in
                  zip(meshes, conns, nodal, painted, titles, beam_ratios)]
    savepair(joinpath(out, "part_meshes"),
             plot(panels...; layout = (2, 2), size = (920, 860)))

    hist = plot(; layout = (1, 2), size = (900, 320))
    for mode in ("gcee", "tmop")
        rr = [r for r in part if r["mode"] == mode]
        plot!(hist, [r["cycle"] for r in rr], [r["E_h"] for r in rr]; subplot = 1,
              marker = :circle, markersize = 3, color = MODE_COLOR[mode],
              label = SHORT[mode], xlabel = "cycle", ylabel = "\$|u-u_h|_{H^1}\$",
              title = "true error", legend = :topleft)
    end
    relaxed_beam = [r for r in part_relax if r["mode"] == "gcee"]
    plot!(hist, [r["cycle"] for r in relaxed_beam], [r["E_h"] for r in relaxed_beam];
          subplot = 1, marker = :diamond, markersize = 3, linestyle = :dash,
          color = C_GCEE, label = "GCEE, relaxed")
    gg = [r for r in part if r["mode"] == "gcee"]
    for (key, lab, color) in (("E_h", "true error", :black),
                              ("gcee_recovered", "GCEE recovered", C_GCEE),
                              ("zz", "ZZ", C_TMOPU),
                              ("residual_estimator", "residual", C_TMOPH))
        plot!(hist, [r["cycle"] for r in gg], [r[key] / gg[1][key] for r in gg];
              subplot = 2, marker = :circle, markersize = 3, color = color, label = lab,
              xlabel = "cycle", ylabel = "normalised by cycle 0",
              title = "GCEE run, indicators", legend = :bottomleft)
    end
    savepair(joinpath(out, "part_history"), hist)

    control = plot(; layout = (2, 2), size = (860, 560))
    free_beam = [r for r in part if r["mode"] == "gcee"]
    tmop_beam = [r for r in part if r["mode"] == "tmop"]
    for (sp, key, ylab, title) in ((1, "E_h", "\$|u-u_h|_{H^1}\$", "curved beam, error"),
                                   (2, "det_ratio", "Jacobian ratio", "curved beam, ratio"))
        for (rows, style, marker, lab, color) in (
                (free_beam, :solid, :circle, "GCEE", C_GCEE),
                (relaxed_beam, :dash, :diamond, "GCEE, relaxed", C_GCEE),
                (tmop_beam, :solid, :square, "TMOP", C_TMOPU))
            plot!(control, [r["cycle"] for r in rows], [r[key] for r in rows];
                  subplot = sp, marker = marker, markersize = 3, linestyle = style,
                  color = color, label = sp == 1 ? lab : "",
                  xlabel = "cycle", ylabel = ylab, title = title,
                  legend = sp == 1 ? :topright : false)
        end
    end
    layer_free = sort([r for r in rel_free if r["n"] == 12 && r["field"] == "layer"],
                      by = r -> r["cycle"])
    layer_relax = sort([r for r in rel if r["n"] == 12 && r["field"] == "layer" &&
                        r["mode"] == "gcee"], by = r -> r["cycle"])
    for (sp, key, ylab, title) in ((3, "E_h_ratio", "\$E_h/E_h^0\$", "layer, n = 12, error"),
                                   (4, "det_ratio", "Jacobian ratio", "layer, n = 12, ratio"))
        for (rows, style, marker, lab) in (
                (layer_free, :solid, :circle, "GCEE"),
                (layer_relax, :dash, :diamond, "GCEE, relaxed"))
            plot!(control, [r["cycle"] for r in rows], [r[key] for r in rows];
                  subplot = sp, marker = marker, markersize = 3, linestyle = style,
                  color = C_GCEE, label = sp == 3 ? lab : "",
                  xlabel = "cycle", ylabel = ylab, title = title,
                  legend = sp == 3 ? :topright : false)
        end
    end
    savepair(joinpath(out, "gcee_control"), control)

    u_beam = manufactured("beam")
    beam_lims = ((-0.05, 2.05), (-0.05, 2.05))
    before_pts, after_pts = partmesh["gcee"][end], partmesh_relax["gcee"][end]
    before_row = free_beam[end]
    after_row = relaxed_beam[end]
    beam_nodal = [[u_beam(pts[i, 1], pts[i, 2]) for i in axes(pts, 1)]
                  for pts in (before_pts, after_pts)]
    beam_paint = [cell_means(partcells, val) for val in beam_nodal]
    bmin, bmax = minimum(minimum.(beam_paint)), maximum(maximum.(beam_paint))
    beam_titles = (
        @sprintf("GCEE\nu max = %.3e    ratio = %.1f", maximum(beam_paint[1]), before_row["det_ratio"]),
        @sprintf("GCEE, relaxed\nu max = %.3e    ratio = %.1f", maximum(beam_paint[2]), after_row["det_ratio"]),
    )
    beam_panels = [mesh_field(pts, partcells, val; title = title, lims = beam_lims,
                              cmap = :YlOrRd, vmin = bmin, vmax = bmax)
                   for (pts, val, title) in zip((before_pts, after_pts), beam_nodal, beam_titles)]
    savepair(joinpath(out, "part_control_meshes"),
             plot(beam_panels...; layout = (1, 2), size = (920, 460)))

    # --- raw data ------------------------------------------------------------
    open(joinpath(DATA, "part.csv"), "w") do io
        println(io, "mode,cycle,E_h,E_I,relative,gcee_recovered,gcee_oracle,zz,residual_estimator,det_min,det_ratio")
        for r in part
            @printf(io, "%s,%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                r["mode"], r["cycle"], r["E_h"], r["E_I"], r["relative"],
                r["gcee_recovered"], r["gcee_oracle"], r["zz"],
                r["residual_estimator"], r["det_min"], r["det_ratio"])
        end
    end
    if !control_only
        if quick
            export_shortlist(out; levels = (4, 8), ngeo = 6, iterations = 3, order = 2)
        else
            export_shortlist(out)
        end
    end
    println("wrote tables and figures to ", out)
end

main()
