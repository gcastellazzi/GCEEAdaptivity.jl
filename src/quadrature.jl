# Quadrilateral meshes on the unit square, with both FinEtools' [-1, 1]²
# reference and the unit-square reference used by MFEM / TMOP.
#
# If η ∈ [-1, 1]² and ξ = (η + 1) / 2 ∈ [0, 1]², then
# ∂/∂ξ = 2 ∂/∂η and dξ = dη / 4. Element area is det(∂x/∂ξ).

"""
    QuadPatch

Cartesian `n`×`n` mesh of bilinear quads on `[0, 1]²`, plus a tensor-product
Gauss rule. Basis values are stored at the quadrature points.
"""
struct QuadPatch
    n::Int
    conn::Matrix{Int}
    xyz0::Matrix{Float64}
    N::Vector{Vector{Float64}}
    dNη::Vector{Matrix{Float64}}
    dNξ::Vector{Matrix{Float64}}
    wη::Vector{Float64}
    wξ::Vector{Float64}
end

"""Number of quadrature points."""
nquad(patch::QuadPatch) = length(patch.wξ)

"""Number of elements."""
nelems(patch::QuadPatch) = size(patch.conn, 1)

"""
    unit_square_mesh(n; order=4)

`n`×`n` mesh of `FESetQ4` elements on `[0, 1]²`, built with [`Q4block`](@ref).

`order` is the number of Gauss points per direction (FinEtools `GaussRule`).
"""
function unit_square_mesh(n::Integer; order::Integer=4)
    n >= 2 || error("need at least a 2×2 mesh")
    fens, fes = Q4block(1.0, 1.0, n, n)
    rule = GaussRule(2, order)
    nq = rule.npts
    N = Vector{Vector{Float64}}(undef, nq)
    dNη = Vector{Matrix{Float64}}(undef, nq)
    dNξ = Vector{Matrix{Float64}}(undef, nq)
    wη = Vector{Float64}(undef, nq)
    wξ = Vector{Float64}(undef, nq)
    prototype = fes
    for q in 1:nq
        η = vec(rule.param_coords[q, :])
        Nq = vec(bfun(prototype, η))
        Gη = bfundpar(prototype, η)
        N[q] = Nq
        dNη[q] = Gη
        dNξ[q] = 2 .* Gη
        wη[q] = rule.weights[q]
        wξ[q] = rule.weights[q] / 4
    end
    return QuadPatch(n, connasarray(fes), copy(fens.xyz), N, dNη, dNξ, wη, wξ)
end

"""Physical image of quadrature point `(e, q)` and the unit-reference Jacobian `∂x/∂ξ`."""
function element_map(xyz::AbstractMatrix, patch::QuadPatch, e::Integer, q::Integer)
    X = xyz[patch.conn[e, :], :]
    A = X' * patch.dNξ[q]
    x = X' * patch.N[q]
    return A, x
end

"""
    integrate_scalar(xyz, patch, f)

Integrate `f(x, y)` over the mesh. The quadrature is the patch rule, so a
sharp `f` needs a fine mesh or a high `order`.
"""
function integrate_scalar(xyz::AbstractMatrix, patch::QuadPatch, f)
    acc = zero(promote_type(eltype(xyz), Float64))
    for e in 1:nelems(patch)
        for q in 1:nquad(patch)
            A, x = element_map(xyz, patch, e, q)
            acc += f(x[1], x[2]) * det2(A) * patch.wξ[q]
        end
    end
    return acc
end

"""Smallest `det(∂x/∂ξ)` on the patch. This is the signed element area at a quadrature point."""
function minimum_area(xyz::AbstractMatrix, patch::QuadPatch)
    m = typeof(zero(eltype(xyz)))(Inf)
    for e in 1:nelems(patch)
        for q in 1:nquad(patch)
            A, _ = element_map(xyz, patch, e, q)
            m = min(m, det2(A))
        end
    end
    return m
end
