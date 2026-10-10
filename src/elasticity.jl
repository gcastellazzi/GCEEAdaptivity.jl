# Displacement-only plane-strain controls and the generalized Gram spectrum.
#
# No body force and no prescribed traction enter the solve: every non-affine
# probe is the gradient of a harmonic polynomial, so it satisfies the
# homogeneous Navier equations identically and the only data are the boundary
# displacements. The spectrum at the end measures the geometric response of an
# element on the whole equilibrated cubic class, independently of any probe.

"""Plane-strain material matrix for `E = 1`, `ν = 0.3`, engineering shear."""
function material()
    mu = 5 / 13
    lam = 15 / 26
    return [lam+2mu lam 0.0; lam lam+2mu 0.0; 0.0 0.0 mu]
end

"""Engineering strain `(ε_xx, ε_yy, γ_xy)` from a displacement gradient `g[i, j] = ∂u_i/∂x_j`."""
strain(g11, g12, g21, g22) = (g11, g22, g12 + g21)

"""
    elastic_probe(name)

One of the five displacement probes. `"affine"` is a rigid-plus-uniform-strain
field; the other four are gradients of the real and imaginary parts of
`(X + iY)^d` about the centre of the unit square, which are harmonic, so the
displacement is divergence free and Navier equilibrium holds with no load.
"""
function elastic_probe(name::AbstractString)
    if name == "affine"
        return (x, y) -> (x + y / 5, -x / 10 + 2y / 5)
    end
    X(x) = x - 0.5
    Y(y) = y - 0.5
    if name == "quadratic_cos"      # φ = Re(X+iY)³ / 3
        return (x, y) -> (X(x)^2 - Y(y)^2, -2 * X(x) * Y(y))
    elseif name == "quadratic_sin"  # φ = Im(X+iY)³ / 3
        return (x, y) -> (2 * X(x) * Y(y), X(x)^2 - Y(y)^2)
    elseif name == "cubic_cos"      # φ = Re(X+iY)⁴ / 4
        return (x, y) -> (X(x)^3 - 3 * X(x) * Y(y)^2, -3 * X(x)^2 * Y(y) + Y(y)^3)
    elseif name == "cubic_sin"      # φ = Im(X+iY)⁴ / 4
        return (x, y) -> (3 * X(x)^2 * Y(y) - Y(y)^3, X(x)^3 - 3 * X(x) * Y(y)^2)
    end
    error("unknown probe $name")
end

"""Residual of the Navier operator at `(x, y)`; zero for every probe."""
function navier_residual(name::AbstractString, x::Real, y::Real)
    p = elastic_probe(name)
    mu, lam = 5 / 13, 15 / 26
    u1(z) = p(z[1], z[2])[1]
    u2(z) = p(z[1], z[2])[2]
    H1 = ForwardDiff.hessian(u1, [x, y])
    H2 = ForwardDiff.hessian(u2, [x, y])
    # ∂(div u) = (u1_xx + u2_xy, u1_xy + u2_yy)
    d1 = H1[1, 1] + H2[1, 2]
    d2 = H1[1, 2] + H2[2, 2]
    return (mu * (H1[1, 1] + H1[2, 2]) + (lam + mu) * d1,
            mu * (H2[1, 1] + H2[2, 2]) + (lam + mu) * d2)
end

"""Taylor coefficients of both components of a probe, as a `ncells`×2×10 array."""
function probe_coefficients(name::AbstractString, points::AbstractMatrix)
    p = elastic_probe(name)
    C = Array{Float64,3}(undef, size(points, 1), 2, length(EXP))
    for c in 1:2
        f = (x, y) -> p(x, y)[c]
        for i in axes(points, 1)
            C[i, c, :] = field_coefficients(f, points[i, 1], points[i, 2])
        end
    end
    return C
end

"""Nodal values of a probe, as a `npoints`×2 matrix."""
function probe_values(name::AbstractString, points::AbstractMatrix)
    p = elastic_probe(name)
    V = Matrix{Float64}(undef, size(points, 1), 2)
    for i in axes(points, 1)
        V[i, 1], V[i, 2] = p(points[i, 1], points[i, 2])
    end
    return V
end

"""Assembled plane-strain stiffness, reused across probes on one mesh."""
struct PatchSystem
    mesh::Mesh
    K::Matrix{Float64}
    free::Vector{Int}
    fixed::Vector{Int}
end

"""Assemble the displacement-only stiffness of `m`."""
function PatchSystem(m::Mesh; order::Integer=11)
    _, B, _, w, _ = mesh_geometry(m, order)
    D = material()
    ne, na, nq = ncells(m), size(m.cells, 2), size(w, 2)
    nd = 2 * npoints(m)
    K = zeros(nd, nd)
    S = zeros(3, 2na)
    for e in 1:ne
        Ke = zeros(2na, 2na)
        for q in 1:nq
            fill!(S, 0.0)
            for a in 1:na
                S[1, 2a - 1] = B[e, q, a, 1]
                S[2, 2a] = B[e, q, a, 2]
                S[3, 2a - 1] = B[e, q, a, 2]
                S[3, 2a] = B[e, q, a, 1]
            end
            Ke .+= w[e, q] .* (S' * D * S)
        end
        dofs = Vector{Int}(undef, 2na)
        for a in 1:na
            dofs[2a - 1] = 2 * m.cells[e, a] - 1
            dofs[2a] = 2 * m.cells[e, a]
        end
        K[dofs, dofs] .+= Ke
    end
    fixedmask = repeat(m.boundary; inner=2)
    return PatchSystem(m, K, findall(.!fixedmask), findall(fixedmask))
end

"""Solve the patch for one probe. Returns nodal displacements and solver diagnostics."""
function patch_solve(sys::PatchSystem, name::AbstractString)
    v = vec(permutedims(probe_values(name, sys.mesh.points)))
    rhs = -sys.K[sys.free, sys.fixed] * v[sys.fixed]
    v[sys.free] = sys.K[sys.free, sys.free] \ rhs
    forces = sys.K * v
    reactions = zeros(length(v))
    reactions[sys.fixed] = forces[sys.fixed]
    R = reshape(reactions, 2, :)
    info = Dict{String,Float64}(
        "algebraic_residual" => norm(forces[sys.free]) / max(norm(rhs), 1e-30),
        "reaction_norm" => norm(reactions),
        "reaction_resultant" => norm(vec(sum(R; dims=2))))
    return permutedims(reshape(v, 2, :)), info
end

"""Element-wise strain energy of a strain field `S[e, q, :]` under the material matrix."""
function strain_energy(w::AbstractMatrix, S::Array{Float64,3})
    D = material()
    out = zeros(size(S, 1))
    for e in axes(S, 1), q in axes(S, 2)
        s = @view S[e, q, :]
        out[e] += w[e, q] * dot(s, D * s)
    end
    return out
end

"""Discrete strain of a nodal vector field on `m`."""
function discrete_strain(m::Mesh, B::Array{Float64,4}, v::AbstractMatrix)
    ne, na, nq = ncells(m), size(m.cells, 2), size(B, 2)
    S = zeros(ne, nq, 3)
    for e in 1:ne, q in 1:nq
        g11 = g12 = g21 = g22 = 0.0
        for a in 1:na
            id = m.cells[e, a]
            g11 += B[e, q, a, 1] * v[id, 1]
            g12 += B[e, q, a, 2] * v[id, 1]
            g21 += B[e, q, a, 1] * v[id, 2]
            g22 += B[e, q, a, 2] * v[id, 2]
        end
        S[e, q, 1], S[e, q, 2], S[e, q, 3] = strain(g11, g12, g21, g22)
    end
    return S
end

"""
    assess_patch(m, uh, name; order=13)

Interpolation and Galerkin errors in the energy norm, the oracle and recovered
Gram values, the equilibrium defect of the recovered polynomial, and the
recovery budget.
"""
function assess_patch(m::Mesh, uh::AbstractMatrix, name::AbstractString; order::Integer=13)
    p = elastic_probe(name)
    _, B, x, w, _ = mesh_geometry(m, order)
    ne, na, nq = ncells(m), size(m.cells, 2), size(w, 2)
    exact = zeros(ne, nq, 3)
    for e in 1:ne, q in 1:nq
        g = ForwardDiff.jacobian(z -> collect(p(z[1], z[2])), [x[e, q, 1], x[e, q, 2]])
        exact[e, q, 1], exact[e, q, 2], exact[e, q, 3] = strain(g[1, 1], g[1, 2], g[2, 1], g[2, 2])
    end
    nodal = probe_values(name, m.points)
    interpolated = discrete_strain(m, B, nodal)
    computed = discrete_strain(m, B, uh)
    EI2 = strain_energy(w, exact .- interpolated)
    Eh2 = strain_energy(w, exact .- computed)
    distance = strain_energy(w, computed .- interpolated)

    centers = element_centers(m)
    rec = [fit(m.points, uh[:, c], centers) for c in 1:2]
    nod = [fit(m.points, nodal[:, c], centers) for c in 1:2]
    coeff = cat(physical_coefficients(rec[1]), physical_coefficients(rec[2]); dims=3)
    exactfit = cat(physical_coefficients(nod[1]), physical_coefficients(nod[2]); dims=3)
    truth = probe_coefficients(name, centers)

    nb = length(EXP)
    R = Array{Float64,4}(undef, ne, nq, nb, 2)
    dn = Matrix{Float64}(undef, na, 2)
    dq = Matrix{Float64}(undef, nq, 2)
    for e in 1:ne
        for a in 1:na
            dn[a, 1] = m.points[m.cells[e, a], 1] - centers[e, 1]
            dn[a, 2] = m.points[m.cells[e, a], 2] - centers[e, 2]
        end
        for q in 1:nq
            dq[q, 1] = x[e, q, 1] - centers[e, 1]
            dq[q, 2] = x[e, q, 2] - centers[e, 2]
        end
        M = monomials(dn)
        Dx = monomials(dq; derivative=(1, 0))
        Dy = monomials(dq; derivative=(0, 1))
        for q in 1:nq, k in 1:nb
            r1, r2 = Dx[q, k], Dy[q, k]
            for a in 1:na
                r1 -= B[e, q, a, 1] * M[a, k]
                r2 -= B[e, q, a, 2] * M[a, k]
            end
            R[e, q, k, 1], R[e, q, k, 2] = r1, r2
        end
    end

    function residual_energy(A)   # A[e, k, c]
        S = zeros(ne, nq, 3)
        for e in 1:ne, q in 1:nq
            g11 = g12 = g21 = g22 = 0.0
            for k in 1:nb
                g11 += A[e, k, 1] * R[e, q, k, 1]
                g12 += A[e, k, 1] * R[e, q, k, 2]
                g21 += A[e, k, 2] * R[e, q, k, 1]
                g22 += A[e, k, 2] * R[e, q, k, 2]
            end
            S[e, q, 1], S[e, q, 2], S[e, q, 3] = strain(g11, g12, g21, g22)
        end
        return strain_energy(w, S)
    end

    truthA = permutedims(truth, (1, 3, 2))
    oracle2 = residual_energy(truthA)
    recovered2 = residual_energy(coeff)
    budget2 = residual_energy(coeff .- truthA)

    # A fitted polynomial need not satisfy equilibrium even when the exact
    # field has no body force. This is recorded separately from interpolation.
    mu, lam = 5 / 13, 15 / 26
    eqd = zeros(ne)
    for e in 1:ne
        dxx = monomials(dq; derivative=(2, 0))
        acc = 0.0
        for q in 1:nq
            dq[q, 1] = x[e, q, 1] - centers[e, 1]
            dq[q, 2] = x[e, q, 2] - centers[e, 2]
        end
        Dxx = monomials(dq; derivative=(2, 0))
        Dyy = monomials(dq; derivative=(0, 2))
        Dxy = monomials(dq; derivative=(1, 1))
        for q in 1:nq
            f1 = f2 = 0.0
            for k in 1:nb
                f1 += coeff[e, k, 1] * ((lam + 2mu) * Dxx[q, k] + mu * Dyy[q, k]) +
                      coeff[e, k, 2] * (lam + mu) * Dxy[q, k]
                f2 += coeff[e, k, 2] * (mu * Dxx[q, k] + (lam + 2mu) * Dyy[q, k]) +
                      coeff[e, k, 1] * (lam + mu) * Dxy[q, k]
            end
            acc += w[e, q] * (f1 * f1 + f2 * f2)
        end
        eqd[e] = sum(@view w[e, :]) * acc
    end

    EI = sqrt(sum(EI2))
    Eh = sqrt(sum(Eh2))
    s = Dict{String,Any}(
        "E_I" => EI, "E_h" => Eh,
        "exact_energy_norm" => sqrt(sum(strain_energy(w, exact))),
        "oracle" => sqrt(sum(oracle2)), "recovered" => sqrt(sum(recovered2)),
        "recovered_equilibrium_defect" => sqrt(sum(eqd)),
        "recovery_budget" => sqrt(sum(budget2)),
        "exact_nodal_coefficient_error" =>
            norm(exactfit .- truthA) / max(norm(truthA), 1e-30),
        "pythagoras_absolute_defect" => abs(sum(EI2) - sum(Eh2) - sum(distance)),
        "patch_condition_max" => maximum(maximum(r.conditions) for r in rec))
    s["oracle_relative_identity_error"] = abs(s["oracle"] - EI) / max(EI, 1e-12)
    s["rec_I"] = EI > 1e-10 ? s["recovered"] / EI : NaN
    s["rec_h"] = Eh > 1e-10 ? s["recovered"] / Eh : NaN
    return s, Dict("E_I2" => EI2, "E_h2" => Eh2, "oracle2" => oracle2,
                   "recovered2" => recovered2, "recovery_budget2" => budget2)
end

"""Laplacian and first derivatives of the monomial basis, as index maps on [`EXP`](@ref)."""
function _monomial_index(i::Integer, j::Integer)
    (i < 0 || j < 0) && return 0
    for (k, e) in enumerate(EXP)
        e == (i, j) && return k
    end
    return 0
end

"""
    navier_basis()

The eight homogeneous degree-2 and degree-3 displacement fields satisfying the
Navier equations, as `C[component, monomial, mode]`. Affine modes are excluded,
so every equilibrated polynomial field of degree at most three is spanned
modulo an affine field. The basis is the orthonormal null space of the six
equilibrium coefficient equations in the fourteen unknowns.
"""
function navier_basis()
    # Unknowns: a[(c-1)*7 + k] multiplies HIGH monomial k in component c.
    A = zeros(6, 14)
    rows = [(0, 0), (1, 0), (0, 1)]
    for c in 1:2, k in 1:7
        col = (c - 1) * 7 + k
        i, j = EXP[HIGH[k]]
        # Contribution to equation `eq_d = ∇²u_d + (5/2) ∂_d(div u)`.
        for d in 1:2
            poly = Dict{Tuple{Int,Int},Float64}()
            add!(p, a, b, v) = (v != 0 && (p[(a, b)] = get(p, (a, b), 0.0) + v))
            if d == c                      # ∇² u_c
                add!(poly, i - 2, j, i * (i - 1))
                add!(poly, i, j - 2, j * (j - 1))
            end
            # (5/2) ∂_d ∂_c (the c-th term of div u is ∂_c u_c)
            di = (d == 1 ? 1 : 0) + (c == 1 ? 1 : 0)
            dj = (d == 2 ? 1 : 0) + (c == 2 ? 1 : 0)
            fac = _falling(i, di) * _falling(j, dj)
            (i >= di && j >= dj) && add!(poly, i - di, j - dj, 2.5 * fac)
            for (r, key) in enumerate(rows)
                A[(d - 1) * 3 + r, col] += get(poly, key, 0.0)
            end
        end
    end
    F = svd(A; full=true)
    null = F.V[:, (rank(A) + 1):end]
    size(null, 2) == 8 || error("expected an 8-dimensional Navier null space")
    C = zeros(2, length(EXP), 8)
    for mode in 1:8, c in 1:2, k in 1:7
        C[c, HIGH[k], mode] = null[(c - 1) * 7 + k, mode]
    end
    return C
end

"""
    geometric_response(m; order=13)

Generalized Gram spectrum on the equilibrated cubic class. The denominator is
the exact strain energy after removing the mean strain, which quotients out the
arbitrary affine representative of each class. No load, no solve, no recovery.
"""
function geometric_response(m::Mesh; order::Integer=13)
    _, B, x, w, _ = mesh_geometry(m, order)
    C = navier_basis()
    centers = element_centers(m)
    ne, na, nq, nm = ncells(m), size(m.cells, 2), size(w, 2), 8
    D = material()
    G = zeros(ne, nm, nm)
    H = zeros(ne, nm, nm)
    spectra = zeros(ne, nm)
    for e in 1:ne
        scale = sqrt(sum(@view w[e, :]))
        dn = Matrix{Float64}(undef, na, 2)
        dq = Matrix{Float64}(undef, nq, 2)
        for a in 1:na
            dn[a, 1] = (m.points[m.cells[e, a], 1] - centers[e, 1]) / scale
            dn[a, 2] = (m.points[m.cells[e, a], 2] - centers[e, 2]) / scale
        end
        for q in 1:nq
            dq[q, 1] = (x[e, q, 1] - centers[e, 1]) / scale
            dq[q, 2] = (x[e, q, 2] - centers[e, 2]) / scale
        end
        M = monomials(dn)
        Dx = monomials(dq; derivative=(1, 0)) ./ scale
        Dy = monomials(dq; derivative=(0, 1)) ./ scale
        S = zeros(nq, nm, 3)
        Rr = zeros(nq, nm, 3)
        for q in 1:nq, mode in 1:nm
            g11 = g12 = g21 = g22 = 0.0
            h11 = h12 = h21 = h22 = 0.0
            for k in 1:length(EXP)
                c1 = C[1, k, mode]
                c2 = C[2, k, mode]
                g11 += c1 * Dx[q, k]; g12 += c1 * Dy[q, k]
                g21 += c2 * Dx[q, k]; g22 += c2 * Dy[q, k]
            end
            for a in 1:na
                v1 = v2 = 0.0
                for k in 1:length(EXP)
                    v1 += C[1, k, mode] * M[a, k]
                    v2 += C[2, k, mode] * M[a, k]
                end
                h11 += B[e, q, a, 1] * v1; h12 += B[e, q, a, 2] * v1
                h21 += B[e, q, a, 1] * v2; h22 += B[e, q, a, 2] * v2
            end
            S[q, mode, 1], S[q, mode, 2], S[q, mode, 3] = strain(g11, g12, g21, g22)
            i1, i2, i3 = strain(h11, h12, h21, h22)
            Rr[q, mode, 1] = S[q, mode, 1] - i1
            Rr[q, mode, 2] = S[q, mode, 2] - i2
            Rr[q, mode, 3] = S[q, mode, 3] - i3
        end
        total = sum(@view w[e, :])
        Sbar = zeros(nm, 3)
        for q in 1:nq, mode in 1:nm, i in 1:3
            Sbar[mode, i] += w[e, q] * S[q, mode, i] / total
        end
        for a in 1:nm, b in 1:nm
            ga = hb = 0.0
            for q in 1:nq
                ra = @view Rr[q, a, :]
                rb = @view Rr[q, b, :]
                ga += w[e, q] * dot(ra, D * rb)
                s0a = [S[q, a, i] - Sbar[a, i] for i in 1:3]
                s0b = [S[q, b, i] - Sbar[b, i] for i in 1:3]
                hb += w[e, q] * dot(s0a, D * s0b)
            end
            G[e, a, b] = ga
            H[e, a, b] = hb
        end
        Ge = Symmetric((G[e, :, :] + G[e, :, :]') / 2)
        He = Symmetric((H[e, :, :] + H[e, :, :]') / 2)
        spectra[e, :] = max.(eigvals(Ge, He), 0.0)
    end
    return Dict("G" => G, "H" => H, "eigenvalues" => spectra,
                "qmax" => sqrt.(spectra[:, end]),
                "qrms" => sqrt.(vec(sum(spectra; dims=2)) ./ nm))
end

# ---------------------------------------------------------------------------
# Critical modes: the geometric worst case of one element, and its response in
# another interpolation space. The eigenproblem supplies the mode; the error is
# then integrated directly from the polynomial and its interpolant, never from
# a matrix contraction, so the two routes are independent checks of each other.

"""A worst-case equilibrated polynomial mode, in a canonical affine representative."""
struct CriticalMode
    center::Vector{Float64}
    scale::Float64
    coefficients::Matrix{Float64}      # component × monomial
    centroid::Vector{Float64}
    mean_value::Vector{Float64}
    mean_gradient::Matrix{Float64}     # component × direction
end

"""Displacement of the mode at `points`, as a `size(points,1)`×2 matrix."""
function mode_value(m::CriticalMode, points::AbstractMatrix)
    q = (points .- m.center') ./ m.scale
    raw = monomials(q) * m.coefficients'
    return raw .- m.mean_value' .- (points .- m.centroid') * m.mean_gradient'
end

"""Displacement gradient of the mode at `points`, as `g[p, i, d] = ∂u_i/∂x_d`."""
function mode_gradient(m::CriticalMode, points::AbstractMatrix)
    q = (points .- m.center') ./ m.scale
    gx = monomials(q; derivative=(1, 0)) * m.coefficients' ./ m.scale
    gy = monomials(q; derivative=(0, 1)) * m.coefficients' ./ m.scale
    g = Array{Float64,3}(undef, size(points, 1), 2, 2)
    g[:, :, 1] = gx
    g[:, :, 2] = gy
    return g .- reshape(m.mean_gradient, 1, 2, 2)
end

"""Strain of a gradient array `g[p, i, d]`, as a `size(g,1)`×3 matrix."""
function strain_of(g::Array{Float64,3})
    S = Matrix{Float64}(undef, size(g, 1), 3)
    for p in axes(g, 1)
        S[p, 1], S[p, 2], S[p, 3] = strain(g[p, 1, 1], g[p, 1, 2], g[p, 2, 1], g[p, 2, 2])
    end
    return S
end

"""Gradient of the nodal interpolant of a vector field, at the quadrature points of `cell`."""
function _interpolant_gradient(m::Mesh, B::Array{Float64,4}, cell::Integer, v::AbstractMatrix)
    nq, na = size(B, 2), size(m.cells, 2)
    g = zeros(nq, 2, 2)
    for q in 1:nq, i in 1:2, d in 1:2
        acc = 0.0
        for a in 1:na
            acc += B[cell, q, a, d] * v[a, i]
        end
        g[q, i, d] = acc
    end
    return g
end

"""
    critical_mode(m, cell, G, H; order=21)

Largest generalized eigenmode of `(G, H)` on `cell`, put into its canonical
affine representative, with statistics that compare the eigenvalue against a
direct integration of the same mode.
"""
function critical_mode(m::Mesh, cell::Integer, G::AbstractMatrix, H::AbstractMatrix;
        order::Integer=21)
    F = eigen(Symmetric((G + G') / 2), Symmetric((H + H') / 2))
    values, vectors = F.values, F.vectors
    a = vectors[:, end]
    basis3 = navier_basis()
    C = zeros(2, length(EXP))
    for mode in 1:8, c in 1:2, k in 1:length(EXP)
        C[c, k] += basis3[c, k, mode] * a[mode]
    end
    # Sign changes do not alter the response; fix one for plotting.
    if C[argmax(abs.(C))] < 0
        a = -a
        C = -C
    end
    _, B, x, w, _ = mesh_geometry(m, order)
    nq, na = size(w, 2), size(m.cells, 2)
    centers = element_centers(m)
    center = [centers[cell, 1], centers[cell, 2]]
    total = sum(@view w[cell, :])
    scale = sqrt(total)
    xq = Matrix{Float64}(undef, nq, 2)
    for q in 1:nq
        xq[q, 1], xq[q, 2] = x[cell, q, 1], x[cell, q, 2]
    end
    q = (xq .- center') ./ scale
    raw = monomials(q) * C'
    gx = monomials(q; derivative=(1, 0)) * C' ./ scale
    gy = monomials(q; derivative=(0, 1)) * C' ./ scale
    wq = @view w[cell, :]
    centroid = vec(wq' * xq) ./ total
    mean_value = vec(wq' * raw) ./ total
    mean_gradient = hcat(vec(wq' * gx) ./ total, vec(wq' * gy) ./ total)
    mode = CriticalMode(center, scale, C, centroid, mean_value, mean_gradient)

    nodes = m.points[m.cells[cell, :], :]
    S = strain_of(mode_gradient(mode, xq))
    Rg = S .- strain_of(_interpolant_gradient(m, B, cell, mode_value(mode, nodes)))
    D = material()
    norm2 = sum(wq[k] * dot(S[k, :], D * S[k, :]) for k in 1:nq)
    error2 = sum(wq[k] * dot(Rg[k, :], D * Rg[k, :]) for k in 1:nq)

    pieces = Matrix{Float64}[]
    for degree in (2, 3)
        Cd = copy(C)
        for k in 1:length(EXP)
            sum(EXP[k]) == degree || (Cd[:, k] .= 0.0)
        end
        pn = monomials((nodes .- center') ./ scale) * Cd'
        gpx = monomials(q; derivative=(1, 0)) * Cd' ./ scale
        gpy = monomials(q; derivative=(0, 1)) * Cd' ./ scale
        gp = Array{Float64,3}(undef, nq, 2, 2)
        gp[:, :, 1] = gpx
        gp[:, :, 2] = gpy
        push!(pieces, strain_of(gp) .- strain_of(_interpolant_gradient(m, B, cell, pn)))
    end
    e2 = sum(wq[k] * dot(pieces[1][k, :], D * pieces[1][k, :]) for k in 1:nq)
    e3 = sum(wq[k] * dot(pieces[2][k, :], D * pieces[2][k, :]) for k in 1:nq)
    cross = 2 * sum(wq[k] * dot(pieces[1][k, :], D * pieces[2][k, :]) for k in 1:nq)

    stats = Dict{String,Float64}(
        "q_eigen" => sqrt(values[end]),
        "q_direct" => sqrt(error2 / norm2),
        "normalization" => norm2,
        "eigen_residual" => norm(G * a - values[end] * H * a) / max(norm(G * a), 1e-30),
        "relative_top_gap" => (values[end] - values[end - 1]) / values[end],
        "error2" => error2, "quadratic_error2" => e2, "cubic_error2" => e3,
        "cross_error2" => cross,
        "diagonal_degree_ratio" => sqrt((e2 + e3) / error2))
    return mode, stats, values, a
end

"""
    direct_response(m, cell, mode; order=21)

Evaluate one fixed physical polynomial in the interpolation space of another
mesh, by direct integration. Used to check that a mode extracted on one space
cannot exceed the worst case of the space it is evaluated in.
"""
function direct_response(m::Mesh, cell::Integer, mode::CriticalMode; order::Integer=21)
    _, B, x, w, _ = mesh_geometry(m, order)
    nq = size(w, 2)
    xq = Matrix{Float64}(undef, nq, 2)
    for q in 1:nq
        xq[q, 1], xq[q, 2] = x[cell, q, 1], x[cell, q, 2]
    end
    wq = @view w[cell, :]
    nodes = m.points[m.cells[cell, :], :]
    S = strain_of(mode_gradient(mode, xq))
    Rg = S .- strain_of(_interpolant_gradient(m, B, cell, mode_value(mode, nodes)))
    D = material()
    norm2 = sum(wq[k] * dot(S[k, :], D * S[k, :]) for k in 1:nq)
    err2 = sum(wq[k] * dot(Rg[k, :], D * Rg[k, :]) for k in 1:nq)
    return sqrt(err2 / norm2)
end
