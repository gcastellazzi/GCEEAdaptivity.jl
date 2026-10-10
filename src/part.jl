# A structured quadrilateral mesh of a real part: the quarter annulus.
#
# The curved beam of the MacNeal-Harder set is a genuine mechanical shape, and
# both of its long boundaries are curved, so the curved-edge channels of the
# first paper are active rather than incidental. The mesh is structured, which
# keeps the logical coordinates the relocation objective needs, and its
# boundary nodes lie on the exact arcs.
#
# Interior vertices are perturbed to build a deliberately poor initial mesh.
# Interior midside nodes follow their edge and the Q9 centre the four vertices,
# so the distortion is angular and the element edges stay straight inside the
# domain; only the two boundary arcs are curved, and they are exact.

"""Unit square to the quarter annulus; `s` is radial and `t` angular."""
function annulus_map(s::Real, t::Real, inner::Real, outer::Real)
    r = inner + s * (outer - inner)
    a = t * (pi / 2)
    return r * cos(a), r * sin(a)
end

"""
    curved_beam(n; kind="Q8", inner=1.0, outer=2.0, amplitude=0.0, seed=11)

Quarter-annulus mesh with `n` by `n` cells. `amplitude` perturbs the interior
vertices by that fraction of the radial cell size; boundary nodes are
untouched and stay on the exact geometry.
"""
function curved_beam(n::Integer; kind::AbstractString="Q8", inner::Real=1.0,
        outer::Real=2.0, amplitude::Real=0.0, seed::Integer=11)
    kind in ("Q4", "Q8", "Q9") || error("unknown element kind $kind")
    order, keys, cells, logical, boundary, corner = node_layout(n, kind)
    np = length(order)
    points = Matrix{Float64}(undef, np, 2)
    for i in 1:np
        x, y = annulus_map(logical[i, 1], logical[i, 2], inner, outer)
        points[i, 1], points[i, 2] = x, y
    end

    if amplitude != 0
        h = (outer - inner) / n
        for i in 1:np
            (corner[i] && !boundary[i]) || continue
            k = order[i]
            points[i, 1] += jitter(seed, k[1], k[2], 1) * amplitude * h
            points[i, 2] += jitter(seed, k[1], k[2], 2) * amplitude * h
        end
    end

    # Interior midside nodes follow their edge; the centre follows the vertices.
    for i in 1:np
        (corner[i] || boundary[i]) && continue
        k = order[i]
        ox, oy = isodd(k[1]), isodd(k[2])
        ends = if ox && oy                      # Q9 centre
            ((k[1] - 1, k[2] - 1), (k[1] + 1, k[2] - 1),
             (k[1] + 1, k[2] + 1), (k[1] - 1, k[2] + 1))
        elseif ox                               # midside along s
            ((k[1] - 1, k[2]), (k[1] + 1, k[2]))
        else                                    # midside along t
            ((k[1], k[2] - 1), (k[1], k[2] + 1))
        end
        sx = sy = 0.0
        for e in ends
            j = keys[e]
            sx += points[j, 1]
            sy += points[j, 2]
        end
        points[i, 1] = sx / length(ends)
        points[i, 2] = sy / length(ends)
    end

    return Mesh(points, cells, boundary, logical, corner, n, kind,
                "curved_beam", Float64(amplitude), Int(seed))
end
