# Manufactured solutions and their exact data.
#
# Derivatives come from nested forward-mode AD rather than from a symbolic
# package, so a field is defined once, as a plain Julia function, and its
# gradient, load and cubic Taylor coefficients follow. `coefficients` returns
# the ten entries of [`EXP`](@ref), each already divided by `a! b!`, which is
# the oracle the recovery is measured against.

"""
    manufactured(name)

Named manufactured solution as a function `(x, y) -> u`. The families are
those of the first paper, plus `"beam"`, a localized feature at radius `3/2`
and `45°` inside the quarter annulus.
"""
function manufactured(name::AbstractString)
    if name == "paper1"
        return (x, y) -> sin(2x) * cosh(y)
    elseif name == "poisson"
        return (x, y) -> sin(3x) * cosh(y) + x * y * y
    elseif name == "layer"
        return (x, y) -> tanh(10 * (x + y - 1))
    elseif name == "bump"
        return (x, y) -> exp(-60 * (x - 0.65)^2 - 12 * (y - 0.35)^2)
    elseif name == "quadratic"
        return (x, y) -> x * x + 0.3 * x * y + 2 * y * y
    elseif name == "cubic"
        return (x, y) -> x^3 + 2 * x * x * y - 0.6 * x * y * y + y^3 + x - y
    elseif name == "affine"
        return (x, y) -> 1 + 2x - 3y
    elseif name == "beam"
        return (x, y) -> exp(-6 * ((x - 1.06)^2 + (y - 1.06)^2))
    else
        error("unknown field $name")
    end
end

"""Gradient of `u` at `(x, y)`."""
field_gradient(u, x, y) = ForwardDiff.gradient(z -> u(z[1], z[2]), [x, y])

"""Hessian of `u` at `(x, y)`."""
field_hessian(u, x, y) = ForwardDiff.hessian(z -> u(z[1], z[2]), [x, y])

"""Load `f = -Δu` of the manufactured solution."""
function field_force(u, x, y)
    H = field_hessian(u, x, y)
    return -(H[1, 1] + H[2, 2])
end

"""
    field_coefficients(u, x, y)

The ten cubic Taylor coefficients of `u` about `(x, y)`, in the order of
[`EXP`](@ref) and each divided by `a! b!`, so that
`u(x+s, y+t) ≈ Σ c_k s^a t^b`.
"""
function field_coefficients(u, x, y)
    p = [x, y]
    g = ForwardDiff.gradient(z -> u(z[1], z[2]), p)
    H = ForwardDiff.hessian(z -> u(z[1], z[2]), p)
    T = ForwardDiff.jacobian(z -> vec(ForwardDiff.hessian(ζ -> u(ζ[1], ζ[2]), z)), p)
    # vec(H) is column-major: H11, H21, H12, H22.
    uxxx, uxxy, uxyy, uyyy = T[1, 1], T[1, 2], T[2, 2], T[4, 2]
    return [u(x, y), g[1], g[2],
            H[1, 1] / 2, H[1, 2], H[2, 2] / 2,
            uxxx / 6, uxxy / 2, uxyy / 2, uyyy / 6]
end

"""Taylor coefficients at every row of `points`, as a `size(points,1)`×10 matrix."""
function coefficients_at(u, points::AbstractMatrix)
    C = Matrix{Float64}(undef, size(points, 1), length(EXP))
    for i in axes(points, 1)
        C[i, :] = field_coefficients(u, points[i, 1], points[i, 2])
    end
    return C
end
