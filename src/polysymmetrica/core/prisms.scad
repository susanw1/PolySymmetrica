// ---------------------------------------------------------------------------
// PolySymmetrica - Prism generators
// Version: 0.1.0
// Copyright 2026 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>

function _ps_prism_height(edge, height, height_scale) =
    ((is_undef(height) ? edge : height) * height_scale);

function _ps_antiprism_height(n, p, edge, angle, eps=1e-12) =
    let(
        r = _ps_polygram_radius(n, p, edge),
        theta = 180 * p / n + angle,
        d2 = 2 * r * r * (1 - cos(theta)),
        h2 = edge * edge - d2,
        _ok = assert(h2 >= -eps, str("poly_antiprism: invalid angle/edge (no real height), n=", n, " p=", p, " edge=", edge, " angle=", angle))
    )
    (h2 <= 0) ? 0 : sqrt(h2);

// Regular/star prism with {n,p} caps.
//
// - n: number of sides (n >= 3)
// - p: polygon step (1 <= p < n/2, gcd(n,p)=1). p=1 => regular prism.
// - edge: target edge length
// - height: explicit prism height (undef => regular default height=edge)
// - height_scale: multiplier applied to the chosen base height
function poly_prism(n=3, p=1, edge=1, height=undef, height_scale=1) =
    let(
        np = _ps_validate_np(n, p, "poly_prism"),
        n_eff = np[0],
        p_eff = np[1],
        _e_ok = assert(edge > 0, "poly_prism: edge must be > 0"),
        _hs_ok = assert(height_scale > 0, "poly_prism: height_scale must be > 0"),
        h = _ps_prism_height(edge, height, height_scale),
        _h_ok = assert(h > 0, "poly_prism: height must be > 0"),
        r = _ps_polygram_radius(n_eff, p_eff, edge),
        z0 = -h / 2,
        z1 = h / 2,
        bottom = _ps_ngon_ring(n_eff, r, z0, 0),
        top = _ps_ngon_ring(n_eff, r, z1, 0),
        cyc = _ps_polygram_cycle(n_eff, p_eff),
        verts = concat(bottom, top),
        faces_raw = concat(
            [[for (k = [0:1:n_eff-1]) cyc[k]]], // bottom cap ({n,p})
            [[for (k = [0:1:n_eff-1]) n_eff + cyc[k]]], // top cap ({n,p})
            [for (k = [0:1:n_eff-1]) // side quads follow cap-cycle adjacency
                let(
                    a = cyc[k],
                    b = cyc[(k+1) % n_eff]
                )
                [a, b, n_eff + b, n_eff + a]
            ]
        ),
        faces = ps_orient_all_faces_outward(verts, faces_raw),
        ir = _ps_poly_ir(verts, faces),
        e_over_ir = edge / ir
    )
    make_poly(verts, faces, e_over_ir);

// Regular/star antiprism with {n,p} caps and 2n side triangles.
//
// - n: number of sides (n >= 3)
// - p: polygon step (1 <= p < n/2, gcd(n,p)=1). p=1 => regular antiprism.
// - edge: target edge length
// - angle: additive twist offset in degrees relative to exact {n,p} antiprism twist (180*p/n)
//          angle=0 => exact regular/star antiprism
// - height: explicit antiprism height (undef => solved from edge + twist)
// - height_scale: multiplier applied to chosen base height
function poly_antiprism(n=3, p=1, edge=1, angle=0, height=undef, height_scale=1) =
    let(
        np = _ps_validate_np(n, p, "poly_antiprism"),
        n_eff = np[0],
        p_eff = np[1],
        _e_ok = assert(edge > 0, "poly_antiprism: edge must be > 0"),
        _hs_ok = assert(height_scale > 0, "poly_antiprism: height_scale must be > 0"),
        h_base = is_undef(height) ? _ps_antiprism_height(n_eff, p_eff, edge, angle) : height,
        h = h_base * height_scale,
        _h_ok = assert(h > 0, "poly_antiprism: height must be > 0"),
        r = _ps_polygram_radius(n_eff, p_eff, edge),
        z0 = -h / 2,
        z1 = h / 2,
        theta = 180 * p_eff / n_eff + angle,
        bottom = _ps_ngon_ring(n_eff, r, z0, 0),
        top = _ps_ngon_ring(n_eff, r, z1, theta),
        cyc = _ps_polygram_cycle(n_eff, p_eff),
        verts = concat(bottom, top),
        faces_raw = concat(
            [[for (k = [0:1:n_eff-1]) cyc[k]]], // bottom cap ({n,p})
            [[for (k = [0:1:n_eff-1]) n_eff + cyc[k]]], // top cap ({n,p})
            concat(
                // Side strip: two triangles per cap-cycle edge.
                [for (k = [0:1:n_eff-1])
                    let(
                        a = cyc[k],
                        b = cyc[(k+1) % n_eff]
                    )
                    [a, b, n_eff + a]
                ],
                [for (k = [0:1:n_eff-1])
                    let(
                        a = cyc[k],
                        b = cyc[(k+1) % n_eff]
                    )
                    [b, n_eff + b, n_eff + a]
                ]
            )
        ),
        faces = ps_orient_all_faces_outward(verts, faces_raw),
        ir = _ps_poly_ir(verts, faces),
        e_over_ir = edge / ir
    )
    make_poly(verts, faces, e_over_ir);
