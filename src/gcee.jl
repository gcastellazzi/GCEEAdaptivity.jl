# Scalar Q4 Gram estimator for the Poisson problem.
#
# The indicator is the square root of the H¹ residual energy of the quadratic
# Taylor modes (x², xy, y²). On a bilinear quad the geometry enters through the
# element size h_K and the dimensionless hourglass (δ_x, δ_y), the same
# parameters as the one-component Q4 a-priori proxy. The membrane block-Gram
# (plane-stress D, two displacement components) is a different operator and is
# not used here.
#
# The resulting element weight replaces ‖∇u‖ inside W. The metric stays μ₉.

"""Geometry of one bilinear quad: size `h_K`, hourglass `(δx, δy)`, centroid."""
function q4_geometry(X::AbstractMatrix)
    x1x, x1y = X[1, 1], X[1, 2]
    x2x, x2y = X[2, 1], X[2, 2]
    x3x, x3y = X[3, 1], X[3, 2]
    x4x, x4y = X[4, 1], X[4, 2]
    a1x = 0.25 * (-x1x + x2x + x3x - x4x)
    a1y = 0.25 * (-x1y + x2y + x3y - x4y)
    a2x = 0.25 * (-x1x - x2x + x3x + x4x)
    a2y = 0.25 * (-x1y - x2y + x3y + x4y)
    hgx = 0.25 * (x1x - x2x + x3x - x4x)
    hgy = 0.25 * (x1y - x2y + x3y - x4y)
    hK = 2 * sqrt(abs(a1x * a2y - a2x * a1y))
    hK = max(hK, 1e-14)
    return hK, 2 * hgx / hK, 2 * hgy / hK, (0.25 * (x1x + x2x + x3x + x4x), 0.25 * (x1y + x2y + x3y + x4y))
end

"""
    q4_gram_indicator(uxx, uxy, uyy, hK, δx, δy)

Square root of the Q4 Gram energy of the quadratic part of a scalar field.
`uxx`, `uxy`, `uyy` are second derivatives at the element centroid.
"""
function q4_gram_indicator(uxx, uxy, uyy, hK, δx, δy)
    val = hK^4 / 45 * (
        0.25 * uxx^2 * (15 + 60 * δx^2) +
        11 * uxy^2 * (δx^2 + δy^2) +
        0.25 * uyy^2 * (15 + 60 * δy^2)
    )
    return sqrt(max(val, 0.0))
end

"""Hessian of a scalar field at one point."""
function scalar_hessian(u, x, y)
    H = ForwardDiff.hessian(z -> u(z[1], z[2]), [x, y])
    return H[1, 1], H[1, 2], H[2, 2]
end

"""Eigenvalues of a symmetric 2×2 Hessian and a unit eigenvector of `λ₁`."""
function hessian_eigs(uxx, uxy, uyy)
    half = 0.5 * (uxx - uyy)
    disc = sqrt(max(half^2 + uxy^2, 0.0))
    λ1 = 0.5 * (uxx + uyy) + disc
    λ2 = 0.5 * (uxx + uyy) - disc
    if abs(uxy) >= abs(uxx - λ1)
        vx, vy = uxy, λ1 - uxx
    else
        vx, vy = λ1 - uyy, uxy
    end
    nrm = hypot(vx, vy)
    if nrm < 1e-14
        return λ1, λ2, (1.0, 0.0)
    end
    return λ1, λ2, (vx / nrm, vy / nrm)
end

"""
Target of area `ζ` whose short axis follows the stronger principal curvature.

`ρ = 1` when the curvature is negligible next to `κ_gate`. Otherwise `ρ` is
the square root of the eigenvalue ratio, clamped to `aspect_max`.
"""
function curvature_target(ζ, uxx, uxy, uyy, κ_gate; aspect_max=8.0)
    λ1, λ2, (nx, ny) = hessian_eigs(uxx, uxy, uyy)
    a1, a2 = abs(λ1), abs(λ2)
    if a2 > a1
        a1, a2 = a2, a1
        nx, ny = -ny, nx
    end
    if a1 < κ_gate || a1 == 0
        return isotropic_W(ζ)
    end
    ρ = clamp(sqrt(a1 / max(a2, 1e-14 * a1)), 1.0, aspect_max)
    s = sqrt(ζ / ρ)
    L = sqrt(ζ * ρ)
    tx, ty = ny, -nx
    return [L * tx s * nx; L * ty s * ny]
end

"""Signed area of element `e`."""
function element_area(xyz, patch::QuadPatch, e::Integer)
    acc = 0.0
    for q in 1:nquad(patch)
        A, _ = element_map(xyz, patch, e, q)
        acc += det2(A) * patch.wξ[q]
    end
    return acc
end

"""
    gcee_indicators(xyz, patch, u)

One scalar GCEE value per element, from the exact Hessian of `u` at the
centroid and the current quad geometry.
"""
function gcee_indicators(xyz::AbstractMatrix, patch::QuadPatch, u)
    η = Vector{Float64}(undef, nelems(patch))
    for e in 1:nelems(patch)
        X = xyz[patch.conn[e, :], :]
        hK, δx, δy, (cx, cy) = q4_geometry(X)
        uxx, uxy, uyy = scalar_hessian(u, cx, cy)
        η[e] = q4_gram_indicator(uxx, uxy, uyy, hK, δx, δy)
    end
    return η
end

"""
    gcee_target(xyz, patch, u; floor_fraction=0.05, aspect_max=8.0, curvature_fraction=0.25)

Frozen element target built from [`gcee_indicators`](@ref).

The size is equidistributed, `ζ = C / w` with `C = (∑ w_e |K_e|) / N`, and `w`
is the indicator floored at `floor_fraction` times its maximum. The aspect
follows the principal curvature where that curvature is at least
`curvature_fraction` of the mesh maximum; elsewhere `W` is isotropic.
The target is held fixed in reference space (`dW = 0`).
"""
function gcee_target(xyz::AbstractMatrix, patch::QuadPatch, u;
        floor_fraction::Real=0.05, aspect_max::Real=8.0, curvature_fraction::Real=0.25)
    n = nelems(patch)
    η = Vector{Float64}(undef, n)
    areas = Vector{Float64}(undef, n)
    hess = Vector{NTuple{3,Float64}}(undef, n)
    κ = Vector{Float64}(undef, n)
    for e in 1:n
        X = xyz[patch.conn[e, :], :]
        hK, δx, δy, (cx, cy) = q4_geometry(X)
        uxx, uxy, uyy = scalar_hessian(u, cx, cy)
        hess[e] = (uxx, uxy, uyy)
        λ1, λ2, _ = hessian_eigs(uxx, uxy, uyy)
        κ[e] = max(abs(λ1), abs(λ2))
        η[e] = q4_gram_indicator(uxx, uxy, uyy, hK, δx, δy)
        areas[e] = element_area(xyz, patch, e)
    end
    ηmax = maximum(η)
    κmax = maximum(κ)
    Ws = Vector{Matrix{Float64}}(undef, n)
    if ηmax <= 0
        ζ = sum(areas) / n
        W0 = isotropic_W(ζ)
        for e in 1:n
            Ws[e] = W0
        end
    else
        w = max.(η, Float64(floor_fraction) * ηmax)
        C = sum(w[e] * areas[e] for e in 1:n) / n
        gate = Float64(curvature_fraction) * κmax
        for e in 1:n
            uxx, uxy, uyy = hess[e]
            Ws[e] = curvature_target(C / w[e], uxx, uxy, uyy, gate; aspect_max=Float64(aspect_max))
        end
    end
    return function (e, q, x, y)
        return Ws[e]
    end
end
