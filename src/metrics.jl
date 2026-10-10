# TMOP quality metrics, 2D.
#
# `T = A W^{-1}` is the transformation from the target element to the current
# element. Formulas follow the implementations in MFEM `fem/tmop.cpp`
# (EvalW / EvalWMatrixForm), not the shorter header comments where those two
# disagree. `|·|` is the Frobenius norm.

"""Frobenius norm squared of a 2×2 matrix."""
frobenius2(T) = T[1, 1]^2 + T[1, 2]^2 + T[2, 1]^2 + T[2, 2]^2

"""Determinant of a 2×2 matrix."""
det2(T) = T[1, 1] * T[2, 2] - T[1, 2] * T[2, 1]

"""Inverse of a 2×2 matrix."""
function inv2(T)
    d = det2(T)
    return [T[2, 2] -T[1, 2]; -T[2, 1] T[1, 1]] ./ d
end

"""μ₁ = |T|². Non-barrier, no normalization."""
mu1(T) = frobenius2(T)

"""μ₂ = ½ |T|² / det(T) − 1. Barrier shape metric. Zero on scaled rotations."""
mu2(T) = 0.5 * frobenius2(T) / det2(T) - 1

"""μ₇ = |T − T^{-t}|² = |T|² + |T^{-1}|² − 4. Barrier shape+size."""
function mu7(T)
    τ = det2(T)
    I1 = frobenius2(T)
    return I1 * (1 + 1 / (τ * τ)) - 4
end

"""μ₉ = det(T) |T − T^{-t}|². Barrier shape+size."""
function mu9(T)
    τ = det2(T)
    I1 = frobenius2(T)
    I1b = I1 / τ
    return (I1 - 4) * τ + I1b
end

"""μ₁₄ = |T − I|². Shape+size+orientation."""
function mu14(T)
    D = [T[1, 1] - 1 T[1, 2]; T[2, 1] T[2, 2] - 1]
    return frobenius2(D)
end

"""μ₅₅ = (det(T) − 1)². Size."""
mu55(T) = (det2(T) - 1)^2

"""μ₅₆ = ½ (det(T) + 1/det(T)) − 1. Barrier size."""
function mu56(T)
    τ = det2(T)
    return 0.5 * (τ + 1 / τ) - 1
end

"""μ₅₈ = I₁b (I₁b − 2). Barrier shape. Scale invariant."""
function mu58(T)
    I1b = frobenius2(T) / det2(T)
    return I1b * (I1b - 2)
end

"""μ₇₇ = ½ (det(T)² + 1/det(T)²) − 1. Barrier size."""
function mu77(T)
    τ = det2(T)
    return 0.5 * (τ * τ + 1 / (τ * τ)) - 1
end

"""μ₉₄ = μ₂ + 1.5 μ₅₆. Balanced shape+size combo used by MFEM `-mid 94`."""
mu94(T) = mu2(T) + 1.5 * mu56(T)

"""
    mu107(A, W)

Shape+orientation metric from the MFEM `mesh-optimizer` help string,
`(1/2) |A − (|A|/|W|) W|²`. Zero when `A` is a positive multiple of `W`.
"""
function mu107(A, W)
    nA = sqrt(frobenius2(A))
    nW = sqrt(frobenius2(W))
    D = A - (nA / nW) * W
    return 0.5 * frobenius2(D)
end

"""Evaluate metric `id` on the pair `(A, W)`. Metric 107 uses `A` and `W` directly."""
function metric_value(id::Integer, A, W)
    if id == 107
        return mu107(A, W)
    end
    T = A * inv2(W)
    if id == 1
        return mu1(T)
    elseif id == 2
        return mu2(T)
    elseif id == 7
        return mu7(T)
    elseif id == 9
        return mu9(T)
    elseif id == 14
        return mu14(T)
    elseif id == 55
        return mu55(T)
    elseif id == 56
        return mu56(T)
    elseif id == 58
        return mu58(T)
    elseif id == 77
        return mu77(T)
    elseif id == 94
        return mu94(T)
    else
        error("TMOP metric id $id is not implemented")
    end
end
