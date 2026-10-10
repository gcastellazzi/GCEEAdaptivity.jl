# The experiments of the recovered-indicator study, ported from the Python
# prototype so that the whole programme runs from one package.
#
# Experiments 1-3 measure recovery, oracle survival and the two competing
# estimators on the mesh families of the first paper. Experiments 4-5 relocate
# nodes with the frozen recovered field and with two target-matrix controls.
# Experiment 6 does the same on a real part, the quarter-annulus curved beam.

"""One row of the recovery study."""
const RecoveryRow = Dict{String,Any}

"""
    recovery_cases(; levels=(4, 8, 16))

The mesh specifications of Experiments 1-3: three element kinds, the regular,
checkerboard and random families, and the three curved families for `Q8`/`Q9`.
"""
function recovery_cases(; levels=(4, 8, 16))
    cases = Tuple{Int,String,String,Float64,Int}[]
    for kind in ("Q4", "Q8", "Q9")
        specs = Tuple{String,Float64,Int}[("regular", 0.0, 7), ("checkerboard", 0.2, 7),
                                          ("checkerboard", 0.4, 7)]
        for a in (0.1, 0.25), seed in (7, 19)
            push!(specs, ("random", a, seed))
        end
        if kind != "Q4"
            for f in ("curved_A", "curved_B"), a in (0.15, 0.3)
                push!(specs, (f, a, 7))
            end
            push!(specs, ("curved_C", 0.5, 7))
        end
        for n in levels, (family, a, seed) in specs
            push!(cases, (n, kind, family, a, seed))
        end
    end
    return cases
end

"""
    run_recovery_study(; levels=(4, 8, 16), name="paper1", verbose=true)

Experiments 1-3. One row per mesh with the recovery budget, both effectivity
normalisations, the two reference estimators and the Jacobian certificate.
"""
function run_recovery_study(; levels=(4, 8, 16), name::AbstractString="paper1",
        verbose::Bool=true)
    cases = recovery_cases(; levels=levels)
    rows = RecoveryRow[]
    u = manufactured(name)
    for (idx, (n, kind, family, a, seed)) in enumerate(cases)
        m = structured_mesh(n; kind=kind, family=family, amplitude=a, seed=seed)
        lower, certificate = certify(m)
        uh, algebraic = solve_field(m, name; order=7)
        EI, Eh, norm2 = field_errors(m, uh, name; order=9)
        summary, _, _ = diagnostics(m, uh, name; order=7)
        ez = sqrt(sum(zz_indicator(m, uh; order=7)))
        er = sqrt(sum(residual_indicator(m, uh, name; order=7)))
        ei, eh = sqrt(sum(EI)), sqrt(sum(Eh))
        _, B, _, w, _ = mesh_geometry(m, 9)
        ui = [u(m.points[i, 1], m.points[i, 2]) for i in 1:npoints(m)]
        distance2 = 0.0
        for e in 1:ncells(m), q in axes(w, 2)
            d1 = d2 = 0.0
            for aa in axes(m.cells, 2)
                δ = uh[m.cells[e, aa]] - ui[m.cells[e, aa]]
                d1 += B[e, q, aa, 1] * δ
                d2 += B[e, q, aa, 2] * δ
            end
            distance2 += w[e, q] * (d1 * d1 + d2 * d2)
        end
        row = RecoveryRow(
            "n" => n, "kind" => kind, "family" => family, "amplitude" => a,
            "seed" => seed, "dofs" => npoints(m), "E_I" => ei, "E_h" => eh,
            "norm_exact" => sqrt(norm2),
            "galerkin_interpolant_distance" => sqrt(distance2),
            "pythagoras_relative_defect" => abs(ei^2 - eh^2 - distance2) / ei^2,
            "ZZ" => ez, "residual_estimator" => er,
            "rec_I" => summary["recovered"] / ei, "rec_h" => summary["recovered"] / eh,
            "oracle_I" => summary["oracle"] / ei, "oracle_h" => summary["oracle"] / eh,
            "ZZ_I" => ez / ei, "ZZ_h" => ez / eh,
            "residual_I" => er / ei, "residual_h" => er / eh,
            "det_lower" => lower, "certificate" => certificate,
            "algebraic_residual" => algebraic)
        for (k, v) in summary
            row[k] = v
        end
        push!(rows, row)
        verbose && @printf("%3d/%3d  %s %-12s a=%-4g n=%2d  rec/oracle=%.3f  rec/Eh=%.3f\n",
                           idx, length(cases), kind, family, a, n,
                           summary["recovered"] / summary["oracle"], row["rec_h"])
    end
    return rows
end

"""
    relocation_experiment(m, name; mode="gcee", cycles=4, inner=25)

Repeated relocate-and-resolve. Returns the per-cycle records, the frozen-step
summaries and the node coordinates of every mesh in the sequence.
"""
function relocation_experiment(m::Mesh, name::AbstractString="layer";
        mode::AbstractString="gcee", cycles::Integer=4, inner::Integer=25,
        ratio_growth::Real=4)
    records = Dict{String,Any}[]
    frozen = Dict{String,Any}[]
    meshes = [copy(m.points)]
    Dc = shapes("Q4", CORN)[2]
    # `Inf` reproduces the relocation with no shape sub-step and no ratio cap.
    # Anchor and cap stay on the mesh that entered the experiment.
    anchor = jacobian_ratio(m)
    cap = mode == "gcee" && isfinite(ratio_growth) ? ratio_growth * anchor : Inf
    for cycle in 0:cycles
        uh, res = solve_field(m, name; order=9)
        EI, Eh, norm2 = field_errors(m, uh, name; order=11)
        summary, _, rec = diagnostics(m, uh, name; order=9)
        push!(records, Dict{String,Any}(
            "cycle" => cycle, "mode" => mode, "field" => name, "n" => m.n,
            "E_I" => sqrt(sum(EI)), "E_h" => sqrt(sum(Eh)),
            "relative_error" => sqrt(sum(Eh) / norm2),
            "recovered" => summary["recovered"], "oracle" => summary["oracle"],
            "recovery_budget" => summary["recovery_budget"], "residual" => res,
            "det_min" => _corner_minimum(Dc, m.cells, m.points),
            "det_ratio" => jacobian_ratio(m)))
        cycle == cycles && break
        new, history, rejected = relocate_nodes(m, rec; mode=mode, iterations=inner, order=9,
                                                ratio_cap=cap, ratio_anchor=anchor)
        push!(frozen, Dict{String,Any}(
            "cycle" => cycle, "J_before" => history[1].J, "J_after" => history[end].J,
            "rejected_trials" => rejected, "accepted_steps" => length(history) - 1,
            "minimum_det" => minimum(h.det_min for h in history),
            "reverted" => mode == "gcee" && jacobian_ratio(new) > cap))
        m = new
        push!(meshes, copy(m.points))
    end
    return records, frozen, meshes
end

"""
    run_relocation_study(; levels=(8, 12), fields=("layer", "bump"), verbose=true)

Experiments 4-5: the frozen recovered field against the two target-matrix
controls, on random meshes, with four relocate-and-resolve cycles each.
"""
function run_relocation_study(; levels=(8, 12), fields=("layer", "bump"),
        modes=("gcee", "tmop_uniform", "tmop_hessian"), cycles::Integer=4,
        inner::Integer=25, ratio_growth::Real=4, verbose::Bool=true)
    rows = Dict{String,Any}[]
    frozen = Dict{String,Any}[]
    meshes = Dict{String,Vector{Matrix{Float64}}}()
    cells = Dict{String,Matrix{Int}}()
    for n in levels, name in fields, mode in modes
        m = structured_mesh(n; kind="Q4", family="random", amplitude=0.2, seed=7)
        records, steps, seq = relocation_experiment(m, name; mode=mode, cycles=cycles,
                                                    inner=inner, ratio_growth=ratio_growth)
        for r in records
            r["E_h_ratio"] = r["E_h"] / records[1]["E_h"]
            r["E_I_ratio"] = r["E_I"] / records[1]["E_I"]
        end
        append!(rows, records)
        for s in steps
            push!(frozen, merge(Dict{String,Any}("n" => n, "field" => name, "mode" => mode), s))
        end
        key = "$(name)_$(mode)_$(n)"
        meshes[key] = seq
        cells[key] = m.cells
        verbose && @printf("n=%2d %-6s %-13s: Eh ratio=%.3f  EI ratio=%.3f\n",
                           n, name, mode, records[end]["E_h_ratio"], records[end]["E_I_ratio"])
    end
    return rows, frozen, meshes, cells
end

"""
    run_part_study(; n=10, amplitude=0.35, seed=11, cycles=4, inner=25)

Experiment 6: relocation on the quarter-annulus curved beam, with the true
error and three indicators recorded after every re-solve.
"""
function run_part_study(; n::Integer=10, amplitude::Real=0.35, seed::Integer=11,
        cycles::Integer=4, inner::Integer=25, field::AbstractString="beam",
        modes=("gcee", "tmop"), ratio_growth::Real=4, verbose::Bool=true)
    rows = Dict{String,Any}[]
    meshes = Dict{String,Vector{Matrix{Float64}}}()
    cells = Matrix{Int}(undef, 0, 0)
    for mode in modes
        m = curved_beam(n; kind="Q4", amplitude=amplitude, seed=seed)
        cells = m.cells
        seq = [copy(m.points)]
        anchor = jacobian_ratio(m)
        cap = mode == "gcee" && isfinite(ratio_growth) ? ratio_growth * anchor : Inf
        for cycle in 0:cycles
            uh, res = solve_field(m, field; order=9)
            EI, Eh, norm2 = field_errors(m, uh, field; order=11)
            summary, _, rec = diagnostics(m, uh, field; order=9)
            zz = sqrt(sum(zz_indicator(m, uh; order=7)))
            resid = sqrt(sum(residual_indicator(m, uh, field; order=7)))
            det = corner_determinants(m)
            push!(rows, Dict{String,Any}(
                "cycle" => cycle, "mode" => mode,
                "E_h" => sqrt(sum(Eh)), "E_I" => sqrt(sum(EI)),
                "relative" => sqrt(sum(Eh) / norm2),
                "gcee_recovered" => summary["recovered"],
                "gcee_oracle" => summary["oracle"],
                "zz" => zz, "residual_estimator" => resid, "solver_residual" => res,
                "det_min" => minimum(det), "det_ratio" => maximum(det) / minimum(det)))
            cycle == cycles && break
            m, _, _ = relocate_nodes(m, rec; mode=mode, iterations=inner, order=9,
                                     ratio_cap=cap, ratio_anchor=anchor)
            push!(seq, copy(m.points))
        end
        meshes[mode] = seq
        if verbose
            first_row = findfirst(r -> r["mode"] == mode, rows)
            e0 = rows[first_row]["E_h"]
            e1 = rows[end]["E_h"]
            @printf("%-13s: E_h %.6e -> %.6e   (%+.1f%%)\n", mode, e0, e1, 100 * (e1 / e0 - 1))
        end
    end
    return rows, meshes, cells
end
