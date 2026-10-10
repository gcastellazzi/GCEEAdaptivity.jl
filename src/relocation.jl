# Fixed-topology node relocation on a `Mesh`, driven either by the frozen
# recovered Gram field or by a target-matrix shape objective.
#
# The recovered cubic is frozen at the start of a relocation, so the objective
# is a function of the node positions alone and its gradient is exact under
# forward-mode AD. Every trial mesh is tested for a positive Jacobian before
# the energy is evaluated, so no invalid mesh is ever an accepted iterate.
# A GCEE step is then followed by a short shape relaxation. The relaxation
# squares cells without changing the area the GCEE step chose. A step that
# still leaves the corner-Jacobian ratio above the cap is undone.

"""Value and physical gradient of the frozen recovered cubic of element `e` at `X`."""
function recovered_at(rec::Recovery, e::Integer, X1, X2)
    s = rec.scales[e]
    q1 = (X1 - rec.centers[e, 1]) / s
    q2 = (X2 - rec.centers[e, 2]) / s
    T = promote_type(typeof(q1), typeof(q2), Float64)
    val = zero(T)
    d1 = zero(T)
    d2 = zero(T)
    for (k, (i, j)) in enumerate(EXP)
        c = rec.coefficients[e, k]
        val += c * q1^i * q2^j
        i > 0 && (d1 += c * i * q1^(i - 1) * q2^j / s)
        j > 0 && (d2 += c * j * q1^i * q2^(j - 1) / s)
    end
    return val, d1, d2
end

"""
    target_inverse(m, rec, mode)

Inverse target Jacobian of each element. `"tmop_uniform"` asks for the ideal
square at the mesh scale; `"tmop_hessian"` asks for the shape of the recovered
Hessian, normalised to unit determinant, so it adapts shape and not size.
"""
function target_inverse(m::Mesh, rec::Recovery, mode::AbstractString)
    ne = ncells(m)
    out = Array{Float64,3}(undef, ne, 2, 2)
    if mode == "tmop_uniform"
        for e in 1:ne
            out[e, :, :] = 2 * m.n * [1.0 0.0; 0.0 1.0]
        end
        return out
    end
    a = physical_coefficients(rec)
    for e in 1:ne
        H = [2a[e, 4] a[e, 5]; a[e, 5] 2a[e, 6]]
        F = eigen(Symmetric(H))
        floor = 0.05 * max(maximum(abs, F.values), 1e-12)
        λ = max.(abs.(F.values), floor)
        M = F.vectors * Diagonal(sqrt.(λ)) * F.vectors'
        M ./= sqrt(det(M))
        out[e, :, :] = 2 * m.n * M
    end
    return out
end

"""Smallest and largest corner Jacobian determinants of a bilinear mesh."""
function _corner_extent(Dc::Array{Float64,3}, cells::AbstractMatrix, P::AbstractMatrix)
    lo = Inf
    hi = -Inf
    for e in axes(cells, 1), q in 1:4
        J11 = J12 = J21 = J22 = 0.0
        for a in axes(cells, 2)
            px, py = P[cells[e, a], 1], P[cells[e, a], 2]
            J11 += Dc[q, a, 1] * px
            J12 += Dc[q, a, 2] * px
            J21 += Dc[q, a, 1] * py
            J22 += Dc[q, a, 2] * py
        end
        d = J11 * J22 - J12 * J21
        lo = min(lo, d)
        hi = max(hi, d)
    end
    return lo, hi
end

"""Smallest corner Jacobian determinant of a bilinear mesh with the given points."""
function _corner_minimum(Dc::Array{Float64,3}, cells::AbstractMatrix, P::AbstractMatrix)
    lo, _ = _corner_extent(Dc, cells, P)
    return lo
end

"""
    jacobian_ratio(m)

Largest corner Jacobian determinant divided by the smallest. A non-positive
corner returns `Inf`.
"""
function jacobian_ratio(m::Mesh)
    d = corner_determinants(m)
    lo = minimum(d)
    return lo > 0 ? maximum(d) / lo : Inf
end

"""Cap `growth * jacobian_ratio(m)`, anchored on the mesh passed in."""
jacobian_cap(m::Mesh; growth::Real=4) = growth * jacobian_ratio(m)

"""
    relocation_objective(m, rec, mode; order=7)

Returns `(objective, build)`, where `objective(z)` is the energy of the
displaced mesh and `build(z)` its node coordinates. `z` holds the two
displacement components of every free node.
"""
function relocation_objective(m::Mesh, rec::Recovery, mode::AbstractString="gcee";
        order::Integer=7)
    m.kind == "Q4" || error("relocation is restricted to bilinear Q4 geometry")
    N, D, _, weights, _ = reference(m.kind, order)
    free = findall(.!m.boundary)
    ne, na, nq = ncells(m), size(m.cells, 2), length(weights)
    target = startswith(mode, "tmop") ? target_inverse(m, rec, mode) : nothing

    function build(z)
        T = promote_type(eltype(z), Float64)
        P = Matrix{T}(undef, npoints(m), 2)
        P .= m.points
        for (k, i) in enumerate(free)
            P[i, 1] += z[2k - 1]
            P[i, 2] += z[2k]
        end
        return P
    end

    function objective(z)
        P = build(z)
        T = eltype(P)
        acc = zero(T)
        for e in 1:ne
            ids = @view m.cells[e, :]
            for q in 1:nq
                J11 = J12 = J21 = J22 = zero(T)
                x1 = x2 = zero(T)
                for a in 1:na
                    px, py = P[ids[a], 1], P[ids[a], 2]
                    J11 += D[q, a, 1] * px
                    J12 += D[q, a, 2] * px
                    J21 += D[q, a, 1] * py
                    J22 += D[q, a, 2] * py
                    x1 += N[q, a] * px
                    x2 += N[q, a] * py
                end
                dJ = J11 * J22 - J12 * J21
                if target !== nothing
                    t11, t12 = target[e, 1, 1], target[e, 1, 2]
                    t21, t22 = target[e, 2, 1], target[e, 2, 2]
                    T11 = J11 * t11 + J12 * t21
                    T12 = J11 * t12 + J12 * t22
                    T21 = J21 * t11 + J22 * t21
                    T22 = J21 * t12 + J22 * t22
                    dT = T11 * T22 - T12 * T21
                    mu = 0.5 * (T11^2 + T12^2 + T21^2 + T22^2) / dT - 1
                    acc += weights[q] * mu
                else
                    i11, i12 = J22 / dJ, -J12 / dJ
                    i21, i22 = -J21 / dJ, J11 / dJ
                    # Gradient of the frozen cubic at the quadrature point,
                    # against the gradient of its nodal interpolant.
                    _, g1, g2 = recovered_at(rec, e, x1, x2)
                    s1 = zero(T)
                    s2 = zero(T)
                    for a in 1:na
                        v, _, _ = recovered_at(rec, e, P[ids[a], 1], P[ids[a], 2])
                        b1 = D[q, a, 1] * i11 + D[q, a, 2] * i21
                        b2 = D[q, a, 1] * i12 + D[q, a, 2] * i22
                        s1 += b1 * v
                        s2 += b2 * v
                    end
                    acc += weights[q] * dJ * ((g1 - s1)^2 + (g2 - s2)^2)
                end
            end
        end
        return target !== nothing ? acc / (4 * m.n * m.n) : acc
    end

    return objective, build
end

"""
One Armijo step inside the trust box around the frozen mesh. `step_scale`
shortens the proposed step. A trial with a corner determinant below `floor`,
or a Jacobian ratio above `ratio_limit`, is rejected.
"""
function _line_search(z, g, objective, value, detmin, ratiofn;
        n, trust, floor, step_scale=1.0, ratio_limit=Inf)
    direction = -g .* ((0.12 / n) * step_scale) ./ max(maximum(abs, g), 1e-30)
    box = trust / n
    candidate = clamp.(z .+ direction, -box, box)
    direction = candidate .- z
    slope = dot(g, direction)
    slope >= -1e-18 && return z, value, 0, false
    α = 1.0
    rejected = 0
    for _ in 1:28
        trial = z .+ α .* direction
        if detmin(trial) > floor && ratiofn(trial) <= ratio_limit
            v = objective(trial)
            if isfinite(v) && v <= value + 1e-4 * α * slope
                return trial, v, rejected, true
            end
        end
        rejected += 1
        α *= 0.5
    end
    return z, value, rejected, false
end

"""
    relocate_nodes(m, rec; mode="gcee", iterations=25, order=7, trust=0.3,
                   ratio_cap=Inf, ratio_anchor=NaN, relax_steps=5)

Scaled steepest descent with a trust box and an Armijo line search. Returns the
relocated mesh, the iteration history and the number of rejected trials.

For `mode="gcee"` and a finite `ratio_cap`, every accepted step that leaves the
corner-Jacobian ratio above `ratio_anchor` is followed by at most `relax_steps`
shape steps on the uniform-square target. Each of those steps is a quarter of
a GCEE step and may not raise the ratio, so the nodes move: the GCEE step
chooses the area and the shape step squares the cells. The sub-steps stop
once the ratio is back to the anchor. A step that is still above `ratio_cap`
afterwards is undone, and the descent stops. Pass the anchor and the cap from
the mesh the budget is measured on. An unset anchor uses the mesh at the
start of this call.
"""
function relocate_nodes(m::Mesh, rec::Recovery; mode::AbstractString="gcee",
        iterations::Integer=25, order::Integer=7, trust::Real=0.3,
        ratio_cap::Real=Inf, ratio_anchor::Real=NaN, relax_steps::Integer=5)
    objective, build = relocation_objective(m, rec, mode; order=order)
    relax = mode == "gcee" && isfinite(ratio_cap)
    shape = relax ? relocation_objective(m, rec, "tmop_uniform"; order=order)[1] : nothing
    nfree = count(.!m.boundary)
    z = zeros(2 * nfree)
    Dc = shapes("Q4", CORN)[2]
    floor = 0.05 / (4 * m.n * m.n)
    history = NamedTuple{(:iteration, :J, :det_min, :gradient_norm),
                         Tuple{Int,Float64,Float64,Float64}}[]
    rejected = 0
    detmin(zz) = _corner_minimum(Dc, m.cells, build(zz))
    ratiofn(zz) = begin
        lo, hi = _corner_extent(Dc, m.cells, build(zz))
        return lo > 0 ? hi / lo : Inf
    end
    detmin(z) > floor || error("initial mesh violates the chosen Jacobian floor")
    anchor = isnan(ratio_anchor) ? ratiofn(z) : ratio_anchor
    value = objective(z)
    for k in 0:(iterations - 1)
        g = ForwardDiff.gradient(objective, z)
        push!(history, (iteration=k, J=value, det_min=detmin(z),
                        gradient_norm=norm(g)))
        maximum(abs, g) < 1e-12 && break
        z_before, value_before = copy(z), value
        z, value, nrej, accepted = _line_search(z, g, objective, value, detmin, ratiofn;
                                                 n=m.n, trust=trust, floor=floor)
        rejected += nrej
        accepted || break
        if relax && ratiofn(z) > anchor
            for _ in 1:relax_steps
                ratiofn(z) <= anchor && break
                gs = ForwardDiff.gradient(shape, z)
                z, _, nrej, ok = _line_search(z, gs, shape, shape(z), detmin, ratiofn;
                                              n=m.n, trust=trust, floor=floor,
                                              step_scale=0.25, ratio_limit=ratiofn(z))
                rejected += nrej
                ok || break
            end
            value = objective(z)
        end
        if relax && ratiofn(z) > ratio_cap
            z, value = z_before, value_before
            break
        end
    end
    final = objective(z)
    push!(history, (iteration=length(history), J=final, det_min=detmin(z),
                    gradient_norm=NaN))
    return with_points(m, build(z)), history, rejected
end
