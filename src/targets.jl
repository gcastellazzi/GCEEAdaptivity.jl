# Target matrices W for the shortlist in docs/tmop-benchmarks.md.
#
# W is the Jacobian from the unit reference square ξ ∈ [0, 1]² to the target
# element, so det(W) is the target area. Construction follows Dobrev et al.,
# IMR 2019, §6.1, the hr-adaptivity paper (arXiv:2010.02166) §4.1–4.2, and
# the MFEM miniapp `miniapps/meshing/mesh-optimizer.hpp` (master, October 2026).

"""Analytic band used by Dobrev et al., IMR 2019, equation (9). Peak value is `2 tanh(1)`."""
function sinusoidal_eta(x, y)
    a = 10 * (y - 0.5) + sin(4π * x)
    return tanh(a + 1) - tanh(a - 1)
end

"""
    band_size(η, s, ratio)

Local target area `η s + (1 − η) ratio s` of IMR 2019, equation (10) and `W₁`.
`s` is the small area, `ratio` the big-to-small ratio (`α` or `γ`).
"""
band_size(η, s, ratio) = η * s + (1 - η) * ratio * s

"""
    two_level_scale(patch, xyz, η, ratio)

Small area `s` from IMR 2019, equation (10): `Vη / s + (V − Vη) / (ratio s) = N`.

`η` is clamped to `[0, 1]` before integration. The raw sinusoidal band peaks at
`2 tanh(1) ≈ 1.52`; without the clamp the factor under the square root of `W₁`
is negative.
"""
function two_level_scale(patch::QuadPatch, xyz::AbstractMatrix, η, ratio)
    ηc(x, y) = clamp(η(x, y), 0.0, 1.0)
    Vη = integrate_scalar(xyz, patch, ηc)
    V = integrate_scalar(xyz, patch, (x, y) -> 1.0)
    N = nelems(patch)
    s = (Vη + (V - Vη) / ratio) / N
    return s, Vη, V
end

"""Isotropic target `W = √ζ I`."""
isotropic_W(ζ) = [sqrt(ζ) 0.0; 0.0 sqrt(ζ)]

"""
    sinusoidal_target(patch, xyz; variant=:size)

Frozen-in-reference-space targets for the three IMR 2019 §6.1 constructions.
`η` is sampled on the initial mesh and then held fixed, which is the discrete
field with `dW = 0`. It is not the advection remap used for the published
transient.

Variants:
- `:size` — `W₁`, `α = 7`, metric 7
- `:shape` — `W₂`, `β = 3`, metric 58
- `:both` — `W₃`, `γ = 3`, metric 9

`W₃` is the diagonal matrix that meets the prose statement in §6.1: at `η = 1`,
aspect ratio `γ` (more resolution in `y`) and area `s`; at `η = 0`, aspect 1 and
area `γ s`. The size `s` uses equation (10) with that same ratio.
"""
function sinusoidal_target(patch::QuadPatch, xyz::AbstractMatrix; variant::Symbol=:size)
    ηc(x, y) = clamp(sinusoidal_eta(x, y), 0.0, 1.0)
    if variant == :shape
        β = 3.0
        build = function (x, y)
            η = ηc(x, y)
            return [(η * β + (1 - η)) 0.0; 0.0 1.0]
        end
        return freeze_target(patch, xyz, build), 58
    elseif variant == :size
        α = 7.0
        s, _, _ = two_level_scale(patch, xyz, sinusoidal_eta, α)
        build = (x, y) -> isotropic_W(band_size(ηc(x, y), s, α))
        return freeze_target(patch, xyz, build), 7
    elseif variant == :both
        γ = 3.0
        s, _, _ = two_level_scale(patch, xyz, sinusoidal_eta, γ)
        build = function (x, y)
            η = ηc(x, y)
            a = sqrt(γ * s)
            b = a * (η / γ + (1 - η))
            return [a 0.0; 0.0 b]
        end
        return freeze_target(patch, xyz, build), 9
    else
        error("unknown sinusoidal variant $variant")
    end
end

"""Sample `build(x, y)` at the initial quadrature points and ignore later motion."""
function freeze_target(patch::QuadPatch, xyz::AbstractMatrix, build)
    nq = nquad(patch)
    W = Vector{Matrix{Float64}}(undef, nelems(patch) * nq)
    for e in 1:nelems(patch)
        for q in 1:nq
            _, x = element_map(xyz, patch, e, q)
            W[(e - 1) * nq + q] = build(x[1], x[2])
        end
    end
    return function (e, q, x, y)
        return W[(e - 1) * nq + q]
    end
end

"""Distance from the centre of the unit square."""
radius_center(x, y) = hypot(x - 0.5, y - 0.5)

"""Unclamped ring field. Inside the annulus the difference of `tanh` exceeds 1."""
function ring_feature(x, y; r1=0.2, r2=0.3, β=30.0)
    # The floor keeps the Hessian of the radius defined at the ring centre.
    r = max(radius_center(x, y), 1e-8)
    return tanh(β * (r - r1)) - tanh(β * (r - r2))
end

"""
Ring indicator clamped to `[0, 1]`, as in MFEM `HRHessianCoefficient`.
The paper text treats `η` as a value in `[0, 1]`.
"""
function ring_indicator(x, y; r1=0.2, r2=0.3, β=30.0)
    return clamp(ring_feature(x, y; r1=r1, r2=r2, β=β), 0.0, 1.0)
end

"""
    ring_size_target(; r1=0.2, r2=0.3, small=0.001, big=0.01, β=30.0)

`W = √ζ I` with `ζ = η small + (1 − η) big`.

Defaults are hr-adaptivity §4.2 (`r1 = 0.2`, `r2 = 0.3`, `δ = 0.001`, `ψ = 10`).
Pass `r1 = 0.15`, `r2 = 0.35` for MFEM `-tid 9`.
"""
function ring_size_target(; r1=0.2, r2=0.3, small=0.001, big=0.01, β=30.0)
    return function (x, y)
        η = ring_indicator(x, y; r1=r1, r2=r2, β=β)
        return isotropic_W(η * small + (1 - η) * big)
    end
end

"""
    ring_aspect_target()

MFEM `-tid 10` / `HRHessianCoefficient` with `hr_target_type == 1`:
size between radii 0.2 and 0.3, plus an angular aspect ratio.
"""
function ring_aspect_target()
    small = 0.001
    big = 0.01
    r1 = 0.2
    r2 = 0.3
    sf = 30.0
    szfac = 1.0
    asfac = 4.0
    eps2 = szfac / asfac
    eps1 = szfac
    return function (x, y)
        xc = x - 0.5
        yc = y - 0.5
        r = hypot(xc, yc)
        ind = clamp(tanh(sf * (r - r1)) - tanh(sf * (r - r2)), 0.0, 1.0)
        sz = ind * small + (1 - ind) * big
        wgt = 0.5 * (tanh(sf * (r - r1) + 1) - tanh(sf * (r - r2) - 1))
        wgt = clamp(wgt, 0.0, 1.0)
        maxval = eps2 + eps1 * (1 - wgt)^2
        minval = eps1
        avg = 0.5 * (maxval + minval)
        amp = 0.5 * (maxval - minval)
        θ = atan(yc, xc)
        val1 = avg + amp * cos(2θ)
        val2 = avg - amp * cos(2θ)
        s = sqrt(sz)
        return [val1 * s 0.0; 0.0 val2 * s]
    end
end

"""
Poisson circular front.

The hr-adaptivity HTML extract drops the `− r` inside the `arctan`. The radius
in the prose, the ICOSAHOM 2021 slides, and Kolev et al. 2026 all keep it.
`sharpness` is `α = 200` in the hr paper and `β = 20` in the 2026 paper.
"""
function circular_front(x, y; sharpness=20.0, radius=0.7, center=(-0.05, -0.05))
    ρ = hypot(x - center[1], y - center[2])
    return atan(sharpness * (ρ - radius))
end

"""Axis-aligned gradient of [`circular_front`](@ref). The centre lies outside `[0, 1]²`."""
function circular_front_grad(x, y; sharpness=20.0, radius=0.7, center=(-0.05, -0.05))
    dx = x - center[1]
    dy = y - center[2]
    ρ = hypot(dx, dy)
    φ = ρ - radius
    s = sharpness / (1 + (sharpness * φ)^2)
    return s * dx / ρ, s * dy / ρ
end

"""Oblique shock. `slope` is the factor in front of the linear argument; §8.2 does not fix it."""
function oblique_shock(x, y; slope=20.0)
    ψ = x - 0.5 - 0.2 * (y - 0.5)
    return atan(slope * ψ)
end

"""Gradient of [`oblique_shock`](@ref)."""
function oblique_shock_grad(x, y; slope=20.0)
    ψ = x - 0.5 - 0.2 * (y - 0.5)
    s = slope / (1 + (slope * ψ)^2)
    return s, -0.2 * s
end

"""
    equidistributed_isotropic(weight, patch, xyz; floor=1e-8)

`W = √ζ I` with `ζ = C / weight` and `C = (∫ weight dx) / N`, so `∫ dx/ζ = N`
on this mesh. A large weight asks for a small element.
"""
function equidistributed_isotropic(weight, patch::QuadPatch, xyz::AbstractMatrix; floor=1e-8)
    wsafe(x, y) = max(weight(x, y), floor)
    C = integrate_scalar(xyz, patch, wsafe) / nelems(patch)
    return function (x, y)
        return isotropic_W(C / wsafe(x, y))
    end
end

"""
    gradient_shape_size(grad, patch, xyz; aspect_max=8.0, floor=1e-6)

Diagonal target (15) of the hr paper: `W = diag(√(ζ/ρ), √(ζ ρ))`, with `ζ`
equidistributed from `‖∇u‖` and `ρ = |u_x| / |u_y|` clamped to a positive range.
The paper writes `ζ ∝ ‖∇u‖`; the figure and the ring example both put the
*small* elements on the feature, so `ζ` here is the target area and scales as
`1 / ‖∇u‖`. `ρ` uses absolute values so the aspect ratio stays positive.
"""
function gradient_shape_size(grad, patch::QuadPatch, xyz::AbstractMatrix; aspect_max=8.0, floor=1e-6)
    magnitude(x, y) = max(hypot(grad(x, y)...), floor)
    C = integrate_scalar(xyz, patch, magnitude) / nelems(patch)
    return function (x, y)
        gx, gy = grad(x, y)
        g = max(hypot(gx, gy), floor)
        ζ = C / g
        ρ = clamp(max(abs(gx), floor) / max(abs(gy), floor), inv(aspect_max), aspect_max)
        return [sqrt(ζ / ρ) 0.0; 0.0 sqrt(ζ * ρ)]
    end
end

"""
Orientation target of Kolev et al. 2026, §8.2, and of MFEM `-tid 8`:
`θ = π y (1 − y) cos(2π x)`, aspect 1, skew `π/2`.

The target area is the current uniform element area `1/n²`, so the metric sees
orientation rather than a global size change. MFEM's `-tid 8` instead uses the
constant area `0.01`.
"""
function orientation_target(n::Integer)
    area = 1 / (n * n)
    s = sqrt(area)
    return function (x, y)
        θ = π * y * (1 - y) * cos(2π * x)
        c, sn = cos(θ), sin(θ)
        return s .* [c -sn; sn c]
    end
end

"""
Unclamped MFEM `size_indicator`: semicircular annulus centred at `(0, 0.5)`,
radii 0.45 and 0.55, sharpness 30.
"""
function size_feature(x, y)
    r = hypot(x, y - 0.5)
    return 0.5 * (1 + tanh(30 * (r - 0.45))) - 0.5 * (1 + tanh(30 * (r - 0.55)))
end

"""`size_feature` clamped to `[0, 1]`."""
function size_indicator(x, y)
    return clamp(size_feature(x, y), 0.0, 1.0)
end

"""
Target area used by MFEM `ConstructSizeGF` in 2D (`size_ratio = 9`).

`volume_ind` and `volume` are integrals of the indicator and of 1. `NE` is the
element count. Returns `(small, big)`.
"""
function discrete_size_levels(volume_ind, volume, NE; ratio=9.0)
    small = volume_ind / NE + (volume - volume_ind) / (ratio * NE)
    return small, ratio * small
end

"""Ideal-shape target for `-tid 5`: `W = √(size) I` with the `ConstructSizeGF` size."""
function discrete_size_target(patch::QuadPatch, xyz::AbstractMatrix)
    Vη = integrate_scalar(xyz, patch, size_indicator)
    V = integrate_scalar(xyz, patch, (x, y) -> 1.0)
    small, big = discrete_size_levels(Vη, V, nelems(patch))
    return function (x, y)
        η = size_indicator(x, y)
        return isotropic_W(η * small + (1 - η) * big)
    end
end

"""Unclamped MFEM `material_indicator_2d`: sinusoidal interface rotated by 22.5°."""
function material_feature(x, y)
    xc = x - 0.5
    yc = y - 0.5
    th = 22.5 * π / 180
    xn = cos(th) * xc + sin(th) * yc
    yn = -sin(th) * xc + cos(th) * yc
    stretch = 1 / cos(th)
    xc2 = xn / stretch
    yc2 = yn / stretch
    arg = 20 * yc2 + 3 * sin(3π * xc2) + 1
    return tanh(arg)
end

"""`material_feature` clamped to `[0, 1]`, as in the miniapp."""
function material_indicator(x, y)
    return clamp(material_feature(x, y), 0.0, 1.0)
end

"""
    material_interface_target(patch, xyz; ratio=9.0)

Size and aspect built from `‖∇η‖` of [`material_indicator`](@ref).

The blend weight is the gradient magnitude divided by its maximum on the
initial mesh. Target areas follow the IMR two-level rule. The aspect ratio
grows from 1 to 3 across the interface, with the short axis along the
gradient. Columns of `W` are ordered so that `det(W) > 0`.
"""
function material_interface_target(patch::QuadPatch, xyz::AbstractMatrix; ratio::Float64=9.0)
    weight(x, y) = hypot(material_indicator_grad(x, y)...)
    gmax = 0.0
    for e in 1:nelems(patch)
        for q in 1:nquad(patch)
            _, x = element_map(xyz, patch, e, q)
            gmax = max(gmax, weight(x[1], x[2]))
        end
    end
    gmax = max(gmax, 1e-12)
    η(x, y) = clamp(weight(x, y) / gmax, 0.0, 1.0)
    s, _, _ = two_level_scale(patch, xyz, η, ratio)
    return function (x, y)
        gx, gy = material_indicator_grad(x, y)
        gn = hypot(gx, gy)
        blend = clamp(gn / gmax, 0.0, 1.0)
        ζ = band_size(blend, s, ratio)
        ρ = 1 + 2 * blend
        scale = sqrt(ζ)
        if gn < 1e-10
            return isotropic_W(ζ)
        end
        normal = [gx, gy] ./ gn
        # Right-handed frame: tangent × normal = +1, short axis along the gradient.
        tangent = [normal[2], -normal[1]]
        λt = sqrt(ρ)
        λn = 1 / sqrt(ρ)
        col1 = scale * λt * tangent
        col2 = scale * λn * normal
        return [col1[1] col2[1]; col1[2] col2[2]]
    end
end

"""Gradient of the unclamped `material_indicator` argument, times `sech²`."""
function material_indicator_grad(x, y)
    xc = x - 0.5
    yc = y - 0.5
    th = 22.5 * π / 180
    c, s = cos(th), sin(th)
    xn = c * xc + s * yc
    yn = -s * xc + c * yc
    # MFEM divides both rotated coordinates by `stretch = 1/cos(θ)`.
    xc2 = xn * c
    yc2 = yn * c
    arg = 20 * yc2 + 3 * sin(3π * xc2) + 1
    # Slope of the unclamped tanh. Inside (0, 1) this is also the slope of the
    # clamped indicator used by the miniapp.
    dtda = 1 - tanh(arg)^2
    darg_dxc2 = 3 * 3π * cos(3π * xc2)
    darg_dyc2 = 20.0
    dxc2_dx = c * c
    dxc2_dy = c * s
    dyc2_dx = c * (-s)
    dyc2_dy = c * c
    darg_dx = darg_dxc2 * dxc2_dx + darg_dyc2 * dyc2_dx
    darg_dy = darg_dxc2 * dxc2_dy + darg_dyc2 * dyc2_dy
    return dtda * darg_dx, dtda * darg_dy
end
