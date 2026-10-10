# Cubic polynomial recovery on nearest-neighbour patches, and the Gram budget
# that separates its error sources.
#
# A patch is fitted in coordinates scaled by its own radius, which keeps the
# Vandermonde well conditioned without claiming superconvergent patch recovery.
# Exact nodal values and finite element values are fitted on precisely the same
# patch geometry, so their difference isolates the error carried by the
# computed nodal data.

"""Result of a patch fit: scaled coefficients, patch radius, conditioning and size."""
struct Recovery
    centers::Matrix{Float64}
    scales::Vector{Float64}
    coefficients::Matrix{Float64}
    conditions::Vector{Float64}
    counts::Vector{Int}
end

"""Coefficients referred to physical coordinates, undoing the patch scaling."""
function physical_coefficients(r::Recovery)
    degree = [sum(e) for e in EXP]
    return r.coefficients ./ (r.scales .^ permutedims(degree))
end

"""Indices and distances of the `k` nodes of `points` nearest to `c`."""
function nearest(points::AbstractMatrix, c::AbstractVector, k::Integer)
    n = size(points, 1)
    d = Vector{Float64}(undef, n)
    for i in 1:n
        d[i] = hypot(points[i, 1] - c[1], points[i, 2] - c[2])
    end
    ids = partialsortperm(d, 1:min(k, n))
    return collect(ids), d[ids]
end

"""
    fit(points, values, centers; minimum=24)

Weighted least squares fit of the cubic Taylor basis on a nearest-neighbour
patch about each centre. The patch grows until the Vandermonde has full rank.
"""
function fit(points::AbstractMatrix, values::AbstractVector, centers::AbstractMatrix;
        minimum_nodes::Integer=24)
    nc, np = size(centers, 1), size(points, 1)
    nb = length(EXP)
    k0 = min(max(minimum_nodes, 2 * nb), np)
    coeff = Matrix{Float64}(undef, nc, nb)
    scales = Vector{Float64}(undef, nc)
    conditions = Vector{Float64}(undef, nc)
    counts = Vector{Int}(undef, nc)
    for c in 1:nc
        cc = [centers[c, 1], centers[c, 2]]
        k = k0
        while true
            ids, dist = nearest(points, cc, k)
            s = max(dist[end], 1e-14)
            q = Matrix{Float64}(undef, length(ids), 2)
            for (r, i) in enumerate(ids)
                q[r, 1] = (points[i, 1] - cc[1]) / s
                q[r, 2] = (points[i, 2] - cc[2]) / s
            end
            V = monomials(q)
            wgt = @. 1 / sqrt(1 + (dist / s)^2)
            A = V .* wgt
            b = values[ids] .* wgt
            F = svd(A)
            tol = 1e-12 * F.S[1]
            rank = count(>(tol), F.S)
            if rank == nb
                y = F.U' * b
                coeff[c, :] = F.V * (y ./ F.S)
                scales[c] = s
                conditions[c] = F.S[1] / F.S[end]
                counts[c] = k
                break
            end
            k == np && error("polynomial recovery patch is rank deficient")
            k = min(np, k + 10)
        end
    end
    return Recovery(Matrix{Float64}(centers), scales, coeff, conditions, counts)
end

"""Element centroids in physical coordinates."""
function element_centers(m::Mesh)
    N = shapes(m.kind, zeros(1, 2))[1]
    C = Matrix{Float64}(undef, ncells(m), 2)
    for e in 1:ncells(m)
        c1 = c2 = 0.0
        for a in axes(m.cells, 2)
            c1 += N[1, a] * m.points[m.cells[e, a], 1]
            c2 += N[1, a] * m.points[m.cells[e, a], 2]
        end
        C[e, 1], C[e, 2] = c1, c2
    end
    return C
end

"""Patch fit centred on every element centroid."""
element_recovery(m::Mesh, values::AbstractVector) =
    fit(m.points, values, element_centers(m))

"""
    gram(m, centers; order=7)

Local Gram operator of the degree-2 and degree-3 monomial residuals about each
element centre, with the residual gradients `R[e, q, k, i]` and the weights.
"""
function gram(m::Mesh, centers::AbstractMatrix; order::Integer=7)
    _, B, x, w, _ = mesh_geometry(m, order)
    ne, na, nq, nh = ncells(m), size(m.cells, 2), size(w, 2), length(HIGH)
    R = Array{Float64,4}(undef, ne, nq, nh, 2)
    G = Array{Float64,3}(undef, ne, nh, nh)
    nodes = Matrix{Float64}(undef, na, 2)
    delta = Matrix{Float64}(undef, nq, 2)
    for e in 1:ne
        for a in 1:na
            nodes[a, 1] = m.points[m.cells[e, a], 1] - centers[e, 1]
            nodes[a, 2] = m.points[m.cells[e, a], 2] - centers[e, 2]
        end
        for q in 1:nq
            delta[q, 1] = x[e, q, 1] - centers[e, 1]
            delta[q, 2] = x[e, q, 2] - centers[e, 2]
        end
        M = monomials(nodes)[:, HIGH]
        Dx = monomials(delta; derivative=(1, 0))[:, HIGH]
        Dy = monomials(delta; derivative=(0, 1))[:, HIGH]
        for q in 1:nq, k in 1:nh
            r1 = Dx[q, k]
            r2 = Dy[q, k]
            for a in 1:na
                r1 -= B[e, q, a, 1] * M[a, k]
                r2 -= B[e, q, a, 2] * M[a, k]
            end
            R[e, q, k, 1] = r1
            R[e, q, k, 2] = r2
        end
        for k in 1:nh, l in 1:nh
            acc = 0.0
            for q in 1:nq
                acc += w[e, q] * (R[e, q, k, 1] * R[e, q, l, 1] + R[e, q, k, 2] * R[e, q, l, 2])
            end
            G[e, k, l] = acc
        end
    end
    return G, R, w
end

"""Element-wise Gram energy of the coefficient set `C` (one row per element, ten columns)."""
function contraction(C::AbstractMatrix, R::Array{Float64,4}, w::AbstractMatrix)
    ne, nq, nh = size(R, 1), size(R, 2), size(R, 3)
    out = zeros(ne)
    for e in 1:ne
        acc = 0.0
        for q in 1:nq
            r1 = r2 = 0.0
            for k in 1:nh
                c = C[e, HIGH[k]]
                r1 += R[e, q, k, 1] * c
                r2 += R[e, q, k, 2] * c
            end
            acc += w[e, q] * (r1 * r1 + r2 * r2)
        end
        out[e] = acc
    end
    return out
end

"""
    diagnostics(m, uh, name; order=7)

Recovery budget in the Gram seminorm. Returns a summary dictionary, the
element-wise squared contributions, and the patch fit.
"""
function diagnostics(m::Mesh, uh::AbstractVector, name::AbstractString="paper1";
        order::Integer=7)
    u = manufactured(name)
    rec = element_recovery(m, uh)
    exact_nodal_values = [u(m.points[i, 1], m.points[i, 2]) for i in 1:npoints(m)]
    exactnodal = fit(m.points, exact_nodal_values, rec.centers)
    truth = coefficients_at(u, rec.centers)
    a = physical_coefficients(rec)
    ai = physical_coefficients(exactnodal)
    _, R, w = gram(m, rec.centers; order=order)
    locals = Dict(
        "oracle2" => contraction(truth, R, w),
        "recovered2" => contraction(a, R, w),
        "recovery_budget2" => contraction(a .- truth, R, w),
        "sampling_budget2" => contraction(ai .- truth, R, w),
        "fem_data_budget2" => contraction(a .- ai, R, w),
    )
    summary = Dict{String,Float64}()
    for (k, v) in locals
        summary[k[1:end - 1]] = sqrt(sum(v))
    end
    volumes = vec(sum(w; dims=2))
    h = sqrt.(volumes)
    for degree in (2, 3)
        ids = [i for (i, e) in enumerate(EXP) if sum(e) == degree]
        weights = volumes .* h .^ (2 * degree)
        denom = sum(weights .* vec(sum(truth[:, ids] .^ 2; dims=2)))
        for (tag, coef) in (("recovery", a), ("exact_nodal", ai))
            num = sum(weights .* vec(sum((coef[:, ids] .- truth[:, ids]) .^ 2; dims=2)))
            summary["$(tag)_degree$(degree)"] = sqrt(num / max(denom, 1e-30))
        end
    end
    summary["patch_condition_max"] = maximum(rec.conditions)
    summary["patch_nodes_max"] = Float64(maximum(rec.counts))
    summary["budget_inequality_slack"] =
        summary["recovery_budget"] - abs(summary["recovered"] - summary["oracle"])
    return summary, locals, rec
end
