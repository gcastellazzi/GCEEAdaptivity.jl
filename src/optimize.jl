# Node relocation with fixed corners and edges that slide on the boundary
# of the unit square. The step is a scaled gradient with an Armijo line search.
# This is not the Newton / MINRES solver used by MFEM, so published objective
# values are not a regression target.

"""Free degrees of freedom: corners fixed, boundary nodes slide, interior nodes free."""
function boundary_dofs(xyz::AbstractMatrix)
    nodes = Int[]
    comps = Int[]
    for i in axes(xyz, 1)
        x, y = xyz[i, 1], xyz[i, 2]
        onx = x <= 1e-12 || x >= 1 - 1e-12
        ony = y <= 1e-12 || y >= 1 - 1e-12
        if onx && ony
            continue
        elseif onx
            push!(nodes, i)
            push!(comps, 2)
        elseif ony
            push!(nodes, i)
            push!(comps, 1)
        else
            push!(nodes, i)
            push!(comps, 1)
            push!(nodes, i)
            push!(comps, 2)
        end
    end
    return nodes, comps
end

"""Coordinates after a displacement `z` of the free dofs. Generic in the element type of `z`."""
function displaced(xyz0::AbstractMatrix, nodes, comps, z)
    T = promote_type(eltype(xyz0), eltype(z))
    xyz = Matrix{T}(undef, size(xyz0, 1), 2)
    xyz .= xyz0
    for k in eachindex(z)
        xyz[nodes[k], comps[k]] += z[k]
    end
    return xyz
end

"""True when every node lies in `[0, 1]²` and every quadrature Jacobian has positive area."""
function mesh_feasible(xyz::AbstractMatrix, patch::QuadPatch; area_floor=0.0)
    for i in axes(xyz, 1)
        (0 - 1e-10 <= xyz[i, 1] <= 1 + 1e-10) || return false
        (0 - 1e-10 <= xyz[i, 2] <= 1 + 1e-10) || return false
    end
    return minimum_area(xyz, patch) > area_floor
end

"""
    relocate(xyz0, patch, metric_id, target; iterations=20, step=0.2)

Gradient descent on the free nodes. `step` is the largest coordinate change
admitted at `α = 1`, as a fraction of the initial element size `1/n`.

Returns `(xyz, history)` with `history` a vector of objective values, starting
from the initial mesh.
"""
function relocate(xyz0::AbstractMatrix, patch::QuadPatch, metric_id::Integer, target;
        iterations::Integer=20, step::Float64=0.2)
    nodes, comps = boundary_dofs(xyz0)
    h = 1 / patch.n
    area_floor = 1e-8 * h * h
    objective(z) = tmop_functional(displaced(xyz0, nodes, comps, z), patch, metric_id, target)
    z = zeros(length(nodes))
    history = Float64[]
    mesh_feasible(xyz0, patch; area_floor=area_floor) ||
        error("initial mesh is not feasible")
    F = objective(z)
    push!(history, F)
    for _ in 1:iterations
        g = ForwardDiff.gradient(objective, z)
        gmax = maximum(abs, g)
        gmax < 1e-10 && break
        direction = -g .* (step * h / gmax)
        slope = dot(g, direction)
        slope >= 0 && break
        α = 1.0
        accepted = false
        for _attempt in 1:16
            ztrial = z .+ α .* direction
            xyz = displaced(xyz0, nodes, comps, ztrial)
            if mesh_feasible(xyz, patch; area_floor=area_floor)
                Ftrial = objective(ztrial)
                if isfinite(Ftrial) && Ftrial <= F + 1e-4 * α * slope
                    z = ztrial
                    F = Ftrial
                    accepted = true
                    break
                end
            end
            α *= 0.5
        end
        accepted || break
        push!(history, F)
    end
    return displaced(xyz0, nodes, comps, z), history
end
