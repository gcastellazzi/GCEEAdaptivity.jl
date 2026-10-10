# Structured quadrilateral meshes with the distortion families of the first
# paper, plus the quarter-annulus curved beam used by the relocation study.
#
# A `Mesh` carries its own connectivity, boundary flags and logical
# coordinates, so it is independent of the `QuadPatch` used by the TMOP
# shortlist. Node keys are integer pairs on the 2n grid: even entries are
# vertices, odd entries midsides, so Q4, Q8 and Q9 share one numbering.

"""
    jitter(seed, i, j, d)

Deterministic value in `(-1, 1)`, from a SplitMix64 hash of the node indices
and the seed. This replaces a stateful generator so that a mesh depends only
on its arguments, not on the order in which its nodes were visited, and stays
identical across machines and versions.
"""
function jitter(seed::Integer, i::Integer, j::Integer, d::Integer)
    z = (UInt64(seed) * 0x9E3779B97F4A7C15) ⊻ (UInt64(i) * 0xBF58476D1CE4E5B9) ⊻
        (UInt64(j) * 0x94D049BB133111EB) ⊻ (UInt64(d) * 0xD6E8FEB86659FD93)
    z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
    z = z ⊻ (z >> 31)
    return 2.0 * (Float64(z >> 11) / Float64(UInt64(1) << 53)) - 1.0
end

"""
    Mesh

Structured quadrilateral mesh. `points` are physical coordinates, `cells` the
connectivity in the node order of [`basis`](@ref), `boundary` marks nodes on
`∂Ω`, and `logical` holds the parameter coordinates the node came from.
"""
struct Mesh
    points::Matrix{Float64}
    cells::Matrix{Int}
    boundary::Vector{Bool}
    logical::Matrix{Float64}
    corner::Vector{Bool}
    n::Int
    kind::String
    family::String
    amplitude::Float64
    seed::Int
end

"""Same mesh with new node coordinates."""
with_points(m::Mesh, points::AbstractMatrix) =
    Mesh(Matrix{Float64}(points), m.cells, m.boundary, m.logical, m.corner,
         m.n, m.kind, m.family, m.amplitude, m.seed)

"""Number of elements."""
ncells(m::Mesh) = size(m.cells, 1)

"""Number of nodes."""
npoints(m::Mesh) = size(m.points, 1)

"""
    node_layout(n, kind)

Keys, connectivity, boundary and corner flags of an `n`×`n` patch of `kind`
elements. Key `(a, b)` runs over `1:2n+1` in each direction before scaling.
"""
function node_layout(n::Integer, kind::AbstractString)
    nodes = basis(kind)[1]
    keys = Dict{Tuple{Int,Int},Int}()
    order = Tuple{Int,Int}[]
    cells = Matrix{Int}(undef, n * n, size(nodes, 1))
    e = 0
    for i in 0:(n - 1), j in 0:(n - 1)
        e += 1
        for a in axes(nodes, 1)
            key = (2i + Int(nodes[a, 1]) + 1, 2j + Int(nodes[a, 2]) + 1)
            id = get(keys, key, 0)
            if id == 0
                push!(order, key)
                id = length(order)
                keys[key] = id
            end
            cells[e, a] = id
        end
    end
    logical = Matrix{Float64}(undef, length(order), 2)
    boundary = Vector{Bool}(undef, length(order))
    corner = Vector{Bool}(undef, length(order))
    for (i, k) in enumerate(order)
        logical[i, 1] = k[1] / (2n)
        logical[i, 2] = k[2] / (2n)
        boundary[i] = k[1] == 0 || k[1] == 2n || k[2] == 0 || k[2] == 2n
        corner[i] = iseven(k[1]) && iseven(k[2])
    end
    return order, keys, cells, logical, boundary, corner
end

"""
    structured_mesh(n; kind="Q4", family="random", amplitude=0.2, seed=7)

Mesh of the unit square in one of the first paper's families: `"regular"`,
`"random"`, `"checkerboard"`, or the curved families `"curved_A"`,
`"curved_B"`, `"curved_C"`, which need `Q8` or `Q9`.
"""
function structured_mesh(n::Integer; kind::AbstractString="Q4",
        family::AbstractString="random", amplitude::Real=0.2, seed::Integer=7)
    kind in ("Q4", "Q8", "Q9") || error("unknown element kind $kind")
    (startswith(family, "curved") && kind == "Q4") && error("curved sides need Q8 or Q9")
    family in ("regular", "random", "checkerboard", "curved_A", "curved_B", "curved_C") ||
        error("unknown family $family")

    vertices = Array{Float64,3}(undef, n + 1, n + 1, 2)
    for i in 0:n, j in 0:n
        vertices[i + 1, j + 1, 1] = i / n
        vertices[i + 1, j + 1, 2] = j / n
    end
    if family == "random"
        a = amplitude / n
        for i in 2:n, j in 2:n, d in 1:2
            vertices[i, j, d] += jitter(seed, i, j, d) * a
        end
    elseif family == "checkerboard"
        for i in 1:(n - 1), j in 1:(n - 1)
            vertices[i + 1, j + 1, 1] += (-1)^(i + j) * amplitude / (2n)
        end
    end

    order, keys, cells, logical, boundary, corner = node_layout(n, kind)
    nodes = basis(kind)[1]
    Nq4 = shapes("Q4", nodes)[1]
    points = Matrix{Float64}(undef, length(order), 2)
    filled = falses(length(order))
    e = 0
    for i in 0:(n - 1), j in 0:(n - 1)
        e += 1
        X = [vertices[i + 1, j + 1, 1] vertices[i + 1, j + 1, 2];
             vertices[i + 2, j + 1, 1] vertices[i + 2, j + 1, 2];
             vertices[i + 2, j + 2, 1] vertices[i + 2, j + 2, 2];
             vertices[i + 1, j + 2, 1] vertices[i + 1, j + 2, 2]]
        xy = Nq4 * X
        if startswith(family, "curved")
            c = family == "curved_C" ? amplitude / n : amplitude
            # Local nodes 5 and 7 are the midsides at t = j and t = j+1.
            for (local_node, row) in ((5, j), (7, j + 1))
                if 0 < row < n
                    s = (family == "curved_A" && isodd(row)) ? -1.0 : 1.0
                    xy[local_node, 2] += s * c / (2n)
                end
            end
        end
        for a in axes(nodes, 1)
            id = cells[e, a]
            if !filled[id]
                points[id, :] = xy[a, :]
                filled[id] = true
            end
        end
    end
    return Mesh(points, cells, boundary, logical, corner, n, kind, family,
                Float64(amplitude), Int(seed))
end

"""
    mesh_geometry(m, order)

Shape values `N`, physical gradients `B[e, q, a, i]`, quadrature points
`x[e, q, i]`, weighted areas `w[e, q]` and Jacobian determinants.
"""
function mesh_geometry(m::Mesh, order::Integer=7)
    N, D, _, weights, _ = reference(m.kind, order)
    ne, na, nq = ncells(m), size(N, 2), length(weights)
    B = Array{Float64,4}(undef, ne, nq, na, 2)
    x = Array{Float64,3}(undef, ne, nq, 2)
    w = Matrix{Float64}(undef, ne, nq)
    det = Matrix{Float64}(undef, ne, nq)
    for e in 1:ne
        ids = @view m.cells[e, :]
        for q in 1:nq
            J11 = J12 = J21 = J22 = 0.0
            x1 = x2 = 0.0
            for a in 1:na
                px, py = m.points[ids[a], 1], m.points[ids[a], 2]
                J11 += D[q, a, 1] * px
                J12 += D[q, a, 2] * px
                J21 += D[q, a, 1] * py
                J22 += D[q, a, 2] * py
                x1 += N[q, a] * px
                x2 += N[q, a] * py
            end
            dJ = J11 * J22 - J12 * J21
            dJ <= 0 && error("nonpositive quadrature Jacobian in element $e")
            i11, i12, i21, i22 = J22 / dJ, -J12 / dJ, -J21 / dJ, J11 / dJ
            for a in 1:na
                B[e, q, a, 1] = D[q, a, 1] * i11 + D[q, a, 2] * i21
                B[e, q, a, 2] = D[q, a, 1] * i12 + D[q, a, 2] * i22
            end
            x[e, q, 1] = x1
            x[e, q, 2] = x2
            det[e, q] = dJ
            w[e, q] = weights[q] * dJ
        end
    end
    return N, B, x, w, det
end

"""Determinant of the geometric Jacobian at the four corners of every element."""
function corner_determinants(m::Mesh)
    D = shapes(m.kind, CORN)[2]
    out = Matrix{Float64}(undef, ncells(m), 4)
    for e in 1:ncells(m), q in 1:4
        J11 = J12 = J21 = J22 = 0.0
        for a in axes(m.cells, 2)
            px, py = m.points[m.cells[e, a], 1], m.points[m.cells[e, a], 2]
            J11 += D[q, a, 1] * px
            J12 += D[q, a, 2] * px
            J21 += D[q, a, 1] * py
            J22 += D[q, a, 2] * py
        end
        out[e, q] = J11 * J22 - J12 * J21
    end
    return out
end

"""
    certify(m)

Lower bound on the geometric Jacobian. Bilinear geometry is certified by its
corner values; a curved family needs the Bernstein coefficients of the
degree-3 determinant, which bound it from below on the whole element.
"""
function certify(m::Mesh)
    if m.kind == "Q4" || !startswith(m.family, "curved")
        return minimum(corner_determinants(m)), "bilinear-corners"
    end
    t = range(0.0, 1.0; length=4)
    B = [binomial(3, j) * s^j * (1 - s)^(3 - j) for s in t, j in 0:3]
    inv_B = inv(B)
    q = Matrix{Float64}(undef, 16, 2)
    k = 0
    for s in t, r in t
        k += 1
        q[k, 1] = 2s - 1
        q[k, 2] = 2r - 1
    end
    D = shapes(m.kind, q)[2]
    lower = Inf
    for e in 1:ncells(m)
        v = Matrix{Float64}(undef, 4, 4)
        for idx in 1:16
            J11 = J12 = J21 = J22 = 0.0
            for a in axes(m.cells, 2)
                px, py = m.points[m.cells[e, a], 1], m.points[m.cells[e, a], 2]
                J11 += D[idx, a, 1] * px
                J12 += D[idx, a, 2] * px
                J21 += D[idx, a, 1] * py
                J22 += D[idx, a, 2] * py
            end
            v[cld(idx, 4), mod1(idx, 4)] = J11 * J22 - J12 * J21
        end
        lower = min(lower, minimum(inv_B * v * inv_B'))
    end
    lower <= 0 && error("Bernstein certificate inconclusive or invalid")
    return lower, "Bernstein-degree3"
end
