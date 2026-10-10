# Q4 Poisson solve on a QuadPatch, Dirichlet data taken from a manufactured solution.
#
# The load is the exact gradient, ∫ ∇u · ∇φ. For test functions that vanish on
# the boundary this is the weak form of −Δu = −∇²u. The hr paper writes ∇²u = f;
# with the exact boundary values the sign convention does not change the solution.

"""Nodal coefficients of the Q4 solution. Boundary nodes, including corners, are set to `u`."""
function solve_poisson(xyz::AbstractMatrix, patch::QuadPatch, u)
    nn = size(xyz, 1)
    fixed = falses(nn)
    ub = zeros(nn)
    for i in 1:nn
        x, y = xyz[i, 1], xyz[i, 2]
        if x <= 1e-12 || x >= 1 - 1e-12 || y <= 1e-12 || y >= 1 - 1e-12
            fixed[i] = true
            ub[i] = u(x, y)
        end
    end
    K = zeros(nn, nn)
    F = zeros(nn)
    for e in 1:nelems(patch)
        ids = patch.conn[e, :]
        X = xyz[ids, :]
        Ke = zeros(4, 4)
        Fe = zeros(4)
        for q in 1:nquad(patch)
            Gη = patch.dNη[q]
            A = X' * Gη
            dV = det2(A) * patch.wη[q]
            grads = Gη * inv2(A)
            Ke .+= (grads * grads') .* dV
            x = X' * patch.N[q]
            # Load from the exact gradient. For v vanishing on the boundary this
            # is the same as ∫ (−Δu) v, and the gradient is less peaked than the
            # Laplacian of a circular front, so the quadrature holds on a coarse mesh.
            g = ForwardDiff.gradient(z -> u(z[1], z[2]), [x[1], x[2]])
            Fe .+= grads * g .* dV
        end
        K[ids, ids] .+= Ke
        F[ids] .+= Fe
    end
    free = findall(.!fixed)
    uh = zeros(nn)
    uh[fixed] = ub[fixed]
    rhs = F[free] .- K[free, fixed] * ub[fixed]
    uh[free] = K[free, free] \ rhs
    return uh
end

"""
    poisson_errors(xyz, patch, u, uh)

`L²` norm and energy seminorm of `u − u_h`. The energy seminorm is `‖∇(u − u_h)‖`.
"""
function poisson_errors(xyz::AbstractMatrix, patch::QuadPatch, u, uh::AbstractVector)
    l2 = 0.0
    energy = 0.0
    for e in 1:nelems(patch)
        ids = patch.conn[e, :]
        X = xyz[ids, :]
        for q in 1:nquad(patch)
            Gη = patch.dNη[q]
            A = X' * Gη
            dV = det2(A) * patch.wη[q]
            N = patch.N[q]
            x = X' * N
            uhq = dot(N, uh[ids])
            diff = u(x[1], x[2]) - uhq
            l2 += diff^2 * dV
            grads = Gη * inv2(A)
            g_h = grads' * uh[ids]
            g = ForwardDiff.gradient(z -> u(z[1], z[2]), [x[1], x[2]])
            energy += sum((g .- g_h) .^ 2) * dV
        end
    end
    return sqrt(l2), sqrt(energy)
end
