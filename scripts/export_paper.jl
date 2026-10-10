# Runs the recovered-indicator study and writes every table and figure that
# Paper A includes from `generated/`. Replaces the Python prototype.
#
#   julia --project=. scripts/export_paper.jl --paper /path/to/PaperA
#   julia --project=. scripts/export_paper.jl --paper /path/to/PaperA --quick
#
# `--quick` drops the n = 16 recovery level; the paper tables need the full run.

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

"""Filled mesh coloured by an element value, with a colour bar."""
function mesh_field(xyz, conn, nodal; title="", lims=nothing, cmap=:viridis,
        vmin=minimum(nodal), vmax=maximum(nodal))
    p = plot(; aspect_ratio = :equal, legend = false, title = title,
             xlabel = "x", ylabel = "y", size = (420, 400), grid = false,
             framestyle = :box)
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
              linecolor = :black, linewidth = 0.35, colorbar = false, label = "")
    end
    scatter!(p, [NaN, NaN], [NaN, NaN]; marker_z = [vmin, vmax],
             seriescolor = colors, colorbar = true, markersize = 0, label = "")
    return p
end

savepair(path, p) = (savefig(p, path * ".pdf"); savefig(p, path * ".png"); path)

# ---------------------------------------------------------------------------

function main()
    paper = nothing
    quick = false
    for (i, a) in enumerate(ARGS)
        a == "--paper" && (paper = ARGS[i + 1])
        a == "--quick" && (quick = true)
    end
    paper === nothing && error("pass --paper /path/to/PaperA")
    out = joinpath(paper, "generated")
    mkpath(out)
    levels = quick ? (4, 8) : (4, 8, 16)

    # --- Experiments 1-3 ----------------------------------------------------
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

    # --- Experiments 4-5 ----------------------------------------------------
    println("relocation study ...")
    rel, frozen, relmesh, relcells = run_relocation_study(; verbose = false)
    rows = String[]
    table_rows = Vector{String}[]
    for r in rel
        r["cycle"] == 4 || continue
        f = first(a for a in frozen if a["n"] == r["n"] && a["field"] == r["field"] &&
                  a["mode"] == r["mode"] && a["cycle"] == 0)
        c1 = first(a for a in rel if a["n"] == r["n"] && a["field"] == r["field"] &&
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
            rr = sort([r for r in rel if r["n"] == n && r["field"] == fieldname &&
                       r["mode"] == mode], by = r -> r["cycle"])
            plot!(errs, [r["cycle"] for r in rr], [r["E_h_ratio"] for r in rr];
                  subplot = panel, marker = :circle, markersize = 3,
                  color = MODE_COLOR[mode], label = panel == 1 ? MODE_LABEL[mode] : "",
                  xlabel = "cycle", ylabel = "\$E_h/E_h^0\$",
                  title = "$fieldname, n = $n", legend = panel == 1 ? :topright : false)
        end
    end
    savepair(joinpath(out, "relocation_errors"), errs)

    key = "layer_gcee_8"
    seq = relmesh[key]
    panels = [mesh_outline(seq[k], relcells[key];
                           title = k == 1 ? "initial" : "cycle $(k-1)",
                           lims = ((-0.02, 1.02), (-0.02, 1.02)))
              for k in (1, 2, 3, length(seq))]
    savepair(joinpath(out, "relocation_meshes"), plot(panels...; layout = (1, 4), size = (1050, 280)))

    # --- Displacement-only elastic control ----------------------------------
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

    # --- Experiment 6: the curved beam ---------------------------------------
    println("curved beam ...")
    part, partmesh, partcells = run_part_study(; verbose = true)
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
    function beam_panel(points, title)
        nodal = [u(points[i, 1], points[i, 2]) for i in axes(points, 1)]
        return mesh_field(points, partcells, nodal; title = title, lims = lims,
                          cmap = :viridis)
    end
    savepair(joinpath(out, "part_meshes"), plot(
        beam_panel(partmesh["gcee"][1], "initial"),
        beam_panel(partmesh["gcee"][end], "GCEE, four cycles"),
        beam_panel(partmesh["tmop"][end], "geometric, four cycles");
        layout = (1, 3), size = (1120, 380)))

    hist = plot(; layout = (1, 2), size = (900, 320))
    for mode in ("gcee", "tmop")
        rr = [r for r in part if r["mode"] == mode]
        plot!(hist, [r["cycle"] for r in rr], [r["E_h"] for r in rr]; subplot = 1,
              marker = :circle, markersize = 3, color = MODE_COLOR[mode],
              label = SHORT[mode], xlabel = "cycle", ylabel = "\$|u-u_h|_{H^1}\$",
              title = "true error", legend = :topleft)
    end
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
    println("wrote tables and figures to ", out)
end

main()
