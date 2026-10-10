# Integral of a TMOP metric on a quad mesh.
#
# F = Σ_q μ w_ξ det(W), which is the quadrature of ∫_{E_t} μ(T) dx_t
# on the unit reference square. See Dobrev et al., hr-adaptivity, equation (5).

"""
    tmop_functional(xyz, patch, metric_id, target)

`target(e, q, x, y)` returns the 2×2 matrix `W` at element `e`, quadrature
point `q`. The line search rejects steps with a non-positive Jacobian before
this integral is used as a gradient.
"""
function tmop_functional(xyz::AbstractMatrix, patch::QuadPatch, metric_id::Integer, target)
    acc = zero(promote_type(eltype(xyz), Float64))
    for e in 1:nelems(patch)
        for q in 1:nquad(patch)
            A, x = element_map(xyz, patch, e, q)
            W = target(e, q, x[1], x[2])
            acc += metric_value(metric_id, A, W) * det2(W) * patch.wξ[q]
        end
    end
    return acc
end

"""Wrap a physical field `build(x, y) -> W` as `target(e, q, x, y)`."""
field_target(build) = (e, q, x, y) -> build(x, y)

"""
    element_metrics(xyz, patch, metric_id, target)

Quadrature-weighted mean of μ on each element of the current mesh.
"""
function element_metrics(xyz::AbstractMatrix, patch::QuadPatch, metric_id::Integer, target)
    μ = Vector{Float64}(undef, nelems(patch))
    for e in 1:nelems(patch)
        acc = 0.0
        wsum = 0.0
        for q in 1:nquad(patch)
            A, x = element_map(xyz, patch, e, q)
            W = target(e, q, x[1], x[2])
            w = patch.wξ[q]
            acc += metric_value(metric_id, A, W) * w
            wsum += w
        end
        μ[e] = acc / wsum
    end
    return μ
end
