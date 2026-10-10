# Poisson solve and error measures on a `Mesh`.
#
# The load is `f = -Δu` of the manufactured solution and the Dirichlet data are
# the exact nodal values, so the discrete problem is the Galerkin projection of
# a known solution and the true error needs no reference solve.

"""
    solve_field(m, name; order=7)

Nodal values of the finite element solution, and the relative algebraic
residual on the free degrees of freedom.
"""
function solve_field(m::Mesh, name::AbstractString="paper1"; order::Integer=7)
    u = manufactured(name)
    N, B, x, w, _ = mesh_geometry(m, order)
    nn, ne, na, nq = npoints(m), ncells(m), size(m.cells, 2), size(w, 2)
    K = zeros(nn, nn)
    F = zeros(nn)
    for e in 1:ne
        ids = @view m.cells[e, :]
        Ke = zeros(na, na)
        Fe = zeros(na)
        for q in 1:nq
            fq = field_force(u, x[e, q, 1], x[e, q, 2])
            for a in 1:na
                Fe[a] += N[q, a] * w[e, q] * fq
                for b in 1:na
                    Ke[a, b] += w[e, q] * (B[e, q, a, 1] * B[e, q, b, 1] +
                                           B[e, q, a, 2] * B[e, q, b, 2])
                end
            end
        end
        for a in 1:na
            F[ids[a]] += Fe[a]
            for b in 1:na
                K[ids[a], ids[b]] += Ke[a, b]
            end
        end
    end
    uh = [u(m.points[i, 1], m.points[i, 2]) for i in 1:nn]
    free = findall(.!m.boundary)
    fixed = findall(m.boundary)
    rhs = F[free] - K[free, fixed] * uh[fixed]
    uh[free] = K[free, free] \ rhs
    residual = norm((K * uh - F)[free]) / max(norm(rhs), 1e-30)
    return uh, residual
end

"""
    field_errors(m, uh, name; order=9)

Element-wise squared interpolation error, squared Galerkin error, and the
squared energy norm of the exact solution.
"""
function field_errors(m::Mesh, uh::AbstractVector, name::AbstractString="paper1";
        order::Integer=9)
    u = manufactured(name)
    _, B, x, w, _ = mesh_geometry(m, order)
    ne, na, nq = ncells(m), size(m.cells, 2), size(w, 2)
    ui = [u(m.points[i, 1], m.points[i, 2]) for i in 1:npoints(m)]
    EI = zeros(ne)
    Eh = zeros(ne)
    norm2 = 0.0
    for e in 1:ne
        ids = @view m.cells[e, :]
        for q in 1:nq
            g = field_gradient(u, x[e, q, 1], x[e, q, 2])
            gi1 = gi2 = gh1 = gh2 = 0.0
            for a in 1:na
                gi1 += B[e, q, a, 1] * ui[ids[a]]
                gi2 += B[e, q, a, 2] * ui[ids[a]]
                gh1 += B[e, q, a, 1] * uh[ids[a]]
                gh2 += B[e, q, a, 2] * uh[ids[a]]
            end
            EI[e] += w[e, q] * ((g[1] - gi1)^2 + (g[2] - gi2)^2)
            Eh[e] += w[e, q] * ((g[1] - gh1)^2 + (g[2] - gh2)^2)
            norm2 += w[e, q] * (g[1]^2 + g[2]^2)
        end
    end
    return EI, Eh, norm2
end
