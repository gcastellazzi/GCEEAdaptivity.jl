# Two reference indicators for the relocation studies: a ZZ-type recovered
# gradient and an uncalibrated explicit residual estimator.
#
# Neither is calibrated against the true error. They are kept as independent
# witnesses: the ZZ indicator rebuilds a continuous gradient from the same
# nodal values, and the residual indicator uses the interior residual and the
# normal flux jumps, with no boundary terms because the test traces vanish on
# a Dirichlet boundary.

"""
    zz_indicator(m, uh; order=7)

Element-wise squared difference between the patch-recovered continuous
gradient and the finite element gradient.
"""
function zz_indicator(m::Mesh, uh::AbstractVector; order::Integer=7)
    rec = fit(m.points, uh, m.points)
    # EXP starts 1, x, y: columns 2 and 3 divided by the patch scale are the gradient.
    nodal = rec.coefficients[:, 2:3] ./ rec.scales
    N, B, _, w, _ = mesh_geometry(m, order)
    ne, na, nq = ncells(m), size(m.cells, 2), size(w, 2)
    out = zeros(ne)
    for e in 1:ne
        ids = @view m.cells[e, :]
        acc = 0.0
        for q in 1:nq
            g1 = g2 = h1 = h2 = 0.0
            for a in 1:na
                g1 += N[q, a] * nodal[ids[a], 1]
                g2 += N[q, a] * nodal[ids[a], 2]
                h1 += B[e, q, a, 1] * uh[ids[a]]
                h2 += B[e, q, a, 2] * uh[ids[a]]
            end
            acc += w[e, q] * ((g1 - h1)^2 + (g2 - h2)^2)
        end
        out[e] = acc
    end
    return out
end

"""Reference points and unit tangent of side `s` of the reference square."""
function _edge_points(s::Integer, p::AbstractVector)
    n = length(p)
    pts = Matrix{Float64}(undef, n, 2)
    if s == 1
        pts[:, 1] .= p; pts[:, 2] .= -1.0
        t = (1.0, 0.0)
    elseif s == 2
        pts[:, 1] .= 1.0; pts[:, 2] .= p
        t = (0.0, 1.0)
    elseif s == 3
        pts[:, 1] .= -p; pts[:, 2] .= 1.0
        t = (-1.0, 0.0)
    else
        pts[:, 1] .= -1.0; pts[:, 2] .= -p
        t = (0.0, -1.0)
    end
    return pts, t
end

"""Outward-ish normal flux of `uh` along side `s` of element `e`, and the edge metric."""
function _edge_flux(m::Mesh, uh::AbstractVector, e::Integer, s::Integer,
        p::AbstractVector, flip::Bool)
    pts, t = _edge_points(s, p)
    if flip
        pts = pts[end:-1:1, :]
    end
    D = shapes(m.kind, pts)[2]
    na = size(m.cells, 2)
    nq = size(pts, 1)
    flux = Vector{Float64}(undef, nq)
    len = Vector{Float64}(undef, nq)
    for q in 1:nq
        J11 = J12 = J21 = J22 = 0.0
        for a in 1:na
            px, py = m.points[m.cells[e, a], 1], m.points[m.cells[e, a], 2]
            J11 += D[q, a, 1] * px
            J12 += D[q, a, 2] * px
            J21 += D[q, a, 1] * py
            J22 += D[q, a, 2] * py
        end
        dJ = J11 * J22 - J12 * J21
        i11, i12, i21, i22 = J22 / dJ, -J12 / dJ, -J21 / dJ, J11 / dJ
        g1 = g2 = 0.0
        for a in 1:na
            b1 = D[q, a, 1] * i11 + D[q, a, 2] * i21
            b2 = D[q, a, 1] * i12 + D[q, a, 2] * i22
            g1 += b1 * uh[m.cells[e, a]]
            g2 += b2 * uh[m.cells[e, a]]
        end
        tx = J11 * t[1] + J12 * t[2]
        ty = J21 * t[1] + J22 * t[2]
        L = hypot(tx, ty)
        flux[q] = (g1 * ty - g2 * tx) / L
        len[q] = L
    end
    return flux, len
end

"""
    residual_indicator(m, uh, name; order=7)

Element-wise explicit residual estimator: `h_K² ‖f + Δu_h‖²` plus the normal
flux jumps on the interior edges, split evenly between the two elements.
"""
function residual_indicator(m::Mesh, uh::AbstractVector, name::AbstractString="paper1";
        order::Integer=7)
    u = manufactured(name)
    N, D, H, weights, q = reference(m.kind, order)
    ne, na, nq = ncells(m), size(m.cells, 2), length(weights)
    out = zeros(ne)
    for e in 1:ne
        ids = @view m.cells[e, :]
        diameter = 0.0
        for a in 1:na, b in 1:na
            diameter = max(diameter, hypot(m.points[ids[a], 1] - m.points[ids[b], 1],
                                           m.points[ids[a], 2] - m.points[ids[b], 2]))
        end
        acc = 0.0
        for k in 1:nq
            J11 = J12 = J21 = J22 = 0.0
            x1 = x2 = 0.0
            for a in 1:na
                px, py = m.points[ids[a], 1], m.points[ids[a], 2]
                J11 += D[k, a, 1] * px
                J12 += D[k, a, 2] * px
                J21 += D[k, a, 1] * py
                J22 += D[k, a, 2] * py
                x1 += N[k, a] * px
                x2 += N[k, a] * py
            end
            dJ = J11 * J22 - J12 * J21
            iv = (J22 / dJ, -J12 / dJ, -J21 / dJ, J11 / dJ)   # (i11, i12, i21, i22)
            g1 = g2 = 0.0
            for a in 1:na
                g1 += (D[k, a, 1] * iv[1] + D[k, a, 2] * iv[3]) * uh[ids[a]]
                g2 += (D[k, a, 1] * iv[2] + D[k, a, 2] * iv[4]) * uh[ids[a]]
            end
            # Second derivatives in reference coordinates, corrected for the
            # curvature of the map, then pulled back to physical coordinates.
            eff = zeros(2, 2)
            for b in 1:2, c in 1:2
                hr = 0.0
                m1 = m2 = 0.0
                for a in 1:na
                    hr += H[k, a, b, c] * uh[ids[a]]
                    m1 += H[k, a, b, c] * m.points[ids[a], 1]
                    m2 += H[k, a, b, c] * m.points[ids[a], 2]
                end
                eff[b, c] = hr - (g1 * m1 + g2 * m2)
            end
            inv2x2 = [iv[1] iv[2]; iv[3] iv[4]]
            hphys = inv2x2' * eff * inv2x2
            lap = hphys[1, 1] + hphys[2, 2]
            f = field_force(u, x1, x2)
            acc += weights[k] * dJ * (f + lap)^2
        end
        out[e] = diameter^2 * acc
    end

    p, wp = gauss_legendre(order)
    book = Dict{Tuple{Int,Int},NTuple{4,Int}}()
    for e in 1:ne, k in 1:4
        a, b = m.cells[e, k], m.cells[e, mod1(k + 1, 4)]
        key = (min(a, b), max(a, b))
        if !haskey(book, key)
            book[key] = (e, k, a, b)
            continue
        end
        other, side, oa, _ = pop!(book, key)
        f1, l1 = _edge_flux(m, uh, e, k, p, false)
        f2, _ = _edge_flux(m, uh, other, side, p, a != oa)
        jump = f1 .+ f2
        edge_length = dot(wp, l1)
        value = edge_length * dot(wp .* l1, jump .^ 2)
        out[e] += value / 2
        out[other] += value / 2
    end
    return out
end
