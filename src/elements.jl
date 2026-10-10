# Q4 / Q8 / Q9 shape functions on the [-1, 1]² reference square, built by
# inverting the monomial Vandermonde at the nodes. This is the element family
# of the first paper, carried over so that the recovery, relocation and
# curved-beam studies run on the same geometry as the Gram framework.
#
# `EXP` is the cubic Taylor basis, ordered by total degree:
#   (0,0) (1,0) (0,1) (2,0) (1,1) (0,2) (3,0) (2,1) (1,2) (0,3)
# `HIGH` selects the entries of total degree 2 and 3, which are the modes the
# Gram operator acts on.

"""Corner nodes of the reference square, counter-clockwise from `(-1, -1)`."""
const CORN = [-1.0 -1.0; 1.0 -1.0; 1.0 1.0; -1.0 1.0]

"""Midside nodes, each between the corner of the same index and the next."""
const MIDS = [0.0 -1.0; 1.0 0.0; 0.0 1.0; -1.0 0.0]

"""Cubic Taylor exponents, ordered by total degree."""
const EXP = Tuple{Int,Int}[(d - j, j) for d in 0:3 for j in 0:d]

"""Indices into [`EXP`](@ref) of the modes of total degree 2 or 3."""
const HIGH = [i for (i, e) in enumerate(EXP) if sum(e) >= 2]

"""Falling factorial `a! / (a-b)!`, the constant left by differentiating `t^a` `b` times."""
function _falling(a::Integer, b::Integer)
    f = 1
    for k in 0:(b - 1)
        f *= (a - k)
    end
    return f
end

"""
    monomials(points, exponents; derivative=(0, 0))

Row `p`, column `k` is `∂^derivative` of the monomial `exponents[k]` at
`points[p, :]`. A monomial killed by the derivative gives an exact zero.
"""
function monomials(points::AbstractMatrix, exponents=EXP; derivative::Tuple{Int,Int}=(0, 0))
    np = size(points, 1)
    T = promote_type(eltype(points), Float64)
    M = zeros(T, np, length(exponents))
    for (k, e) in enumerate(exponents)
        (e[1] < derivative[1] || e[2] < derivative[2]) && continue
        c = _falling(e[1], derivative[1]) * _falling(e[2], derivative[2])
        a, b = e[1] - derivative[1], e[2] - derivative[2]
        for p in 1:np
            M[p, k] = c * points[p, 1]^a * points[p, 2]^b
        end
    end
    return M
end

"""Nodes and monomial exponents of element `kind`, one of `"Q4"`, `"Q8"`, `"Q9"`."""
function element_basis(kind::AbstractString)
    if kind == "Q4"
        nodes = CORN
        exps = Tuple{Int,Int}[(0, 0), (1, 0), (0, 1), (1, 1)]
    elseif kind == "Q8"
        nodes = vcat(CORN, MIDS)
        exps = Tuple{Int,Int}[(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 1), (1, 2)]
    elseif kind == "Q9"
        nodes = vcat(CORN, MIDS, [0.0 0.0])
        exps = Tuple{Int,Int}[(i, j) for i in 0:2 for j in 0:2]
    else
        error("unknown element kind $kind")
    end
    return nodes, exps, inv(monomials(nodes, exps))
end

const _BASIS = Dict{String,Tuple{Matrix{Float64},Vector{Tuple{Int,Int}},Matrix{Float64}}}()

"""Cached [`element_basis`](@ref)."""
function basis(kind::AbstractString)
    get!(_BASIS, kind) do
        element_basis(kind)
    end
end

"""Number of nodes of element `kind`."""
nnodes(kind::AbstractString) = size(basis(kind)[1], 1)

"""
    shapes(kind, points)

Values `N[p, a]`, first derivatives `D[p, a, d]` and second derivatives
`H[p, a, i, j]` of the shape functions at reference `points`.
"""
function shapes(kind::AbstractString, points::AbstractMatrix)
    _, exps, C = basis(kind)
    np, na = size(points, 1), length(exps)
    N = monomials(points, exps) * C
    D = Array{Float64,3}(undef, np, na, 2)
    D[:, :, 1] = monomials(points, exps; derivative=(1, 0)) * C
    D[:, :, 2] = monomials(points, exps; derivative=(0, 1)) * C
    H = Array{Float64,4}(undef, np, na, 2, 2)
    for i in 1:2, j in 1:2
        d = (i == 1 ? 1 : 0) + (j == 1 ? 1 : 0), (i == 2 ? 1 : 0) + (j == 2 ? 1 : 0)
        H[:, :, i, j] = monomials(points, exps; derivative=d) * C
    end
    return N, D, H
end

"""
    gauss_legendre(n)

Nodes and weights of the `n`-point Gauss-Legendre rule on `[-1, 1]`, from the
Newton iteration on the Legendre polynomial. Self-contained, so the element
quadrature does not depend on the FinEtools rule used by `QuadPatch`.
"""
function gauss_legendre(n::Integer)
    x = zeros(Float64, n)
    w = zeros(Float64, n)
    for i in 1:n
        # Chebyshev start, accurate enough for quadratic convergence below.
        z = cos(pi * (i - 0.25) / (n + 0.5))
        local pp = 0.0
        for _ in 1:100
            p0, p1 = 1.0, 0.0
            for k in 1:n
                p2 = p1
                p1 = p0
                p0 = ((2k - 1) * z * p1 - (k - 1) * p2) / k
            end
            pp = n * (z * p0 - p1) / (z * z - 1)
            dz = p0 / pp
            z -= dz
            abs(dz) < 1e-15 && break
        end
        x[i] = z
        w[i] = 2 / ((1 - z * z) * pp * pp)
    end
    perm = sortperm(x)
    return x[perm], w[perm]
end

"""Tensor-product Gauss-Legendre rule of `order` points per direction on `[-1, 1]²`."""
function tensor_gauss(order::Integer)
    p, w = gauss_legendre(order)
    nq = order * order
    q = Matrix{Float64}(undef, nq, 2)
    weights = Vector{Float64}(undef, nq)
    k = 0
    for i in 1:order, j in 1:order
        k += 1
        q[k, 1] = p[i]
        q[k, 2] = p[j]
        weights[k] = w[i] * w[j]
    end
    return q, weights
end

const _REFERENCE = Dict{Tuple{String,Int},
    Tuple{Matrix{Float64},Array{Float64,3},Array{Float64,4},Vector{Float64},Matrix{Float64}}}()

"""
    reference(kind, order)

Shape values and derivatives of element `kind` at the `order`-point tensor
Gauss rule, with the rule's weights and reference coordinates. Cached.
"""
function reference(kind::AbstractString, order::Integer=7)
    get!(_REFERENCE, (kind, order)) do
        q, w = tensor_gauss(order)
        N, D, H = shapes(kind, q)
        (N, D, H, w, q)
    end
end
