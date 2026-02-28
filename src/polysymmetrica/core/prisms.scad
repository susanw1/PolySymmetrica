// ---------------------------------------------------------------------------
// PolySymmetrica - Prism generators
// Version: 0.1.0
// Copyright 2026 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>

// Regular n-gon radius for given edge length.
function _ps_ngon_radius(n, edge) =
    edge / (2 * sin(180 / n));

function _ps_ngon_ring(n, radius, z, phase=0) =
    [ for (k = [0:1:n-1]) [radius * cos(360*k/n + phase), radius * sin(360*k/n + phase), z] ];

function _ps_poly_ir(verts, faces) =
    let(
        edges = _ps_edges_from_faces(faces),
        mids = [for (e = edges) norm((verts[e[0]] + verts[e[1]]) / 2)],
        ir = min(mids),
        _ok = assert(ir > 0, "prism/antiprism: inter-radius must be > 0")
    ) ir;

function _ps_prism_height(edge, height, height_scale) =
    ((is_undef(height) ? edge : height) * height_scale);

function _ps_antiprism_height(n, edge, angle, eps=1e-12) =
    let(
        r = _ps_ngon_radius(n, edge),
        theta = 180 / n + angle,
        d2 = 2 * r * r * (1 - cos(theta)),
        h2 = edge * edge - d2,
        _ok = assert(h2 >= -eps, str("poly_antiprism: invalid angle/edge (no real height), n=", n, " edge=", edge, " angle=", angle))
    )
    (h2 <= 0) ? 0 : sqrt(h2);

// Regular prism with n-gon caps.
//
// - n: number of sides (n >= 3)
// - edge: target edge length
// - height: explicit prism height (undef => regular default height=edge)
// - height_scale: multiplier applied to the chosen base height
function poly_prism(n=3, edge=1, height=undef, height_scale=1) =
    let(
        _n_ok = assert(n >= 3, "poly_prism: n must be >= 3"),
        _e_ok = assert(edge > 0, "poly_prism: edge must be > 0"),
        _hs_ok = assert(height_scale > 0, "poly_prism: height_scale must be > 0"),
        h = _ps_prism_height(edge, height, height_scale),
        _h_ok = assert(h > 0, "poly_prism: height must be > 0"),
        r = _ps_ngon_radius(n, edge),
        z0 = -h / 2,
        z1 = h / 2,
        bottom = _ps_ngon_ring(n, r, z0, 0),
        top = _ps_ngon_ring(n, r, z1, 0),
        verts = concat(bottom, top),
        faces_raw = concat(
            [[for (k = [0:1:n-1]) k]],                             // bottom
            [[for (k = [0:1:n-1]) n + k]],                         // top
            [for (k = [0:1:n-1]) [k, (k+1)%n, n+((k+1)%n), n+k]] // sides
        ),
        faces = ps_orient_all_faces_outward(verts, faces_raw),
        ir = _ps_poly_ir(verts, faces),
        e_over_ir = edge / ir
    )
    make_poly(verts, faces, e_over_ir);

// Regular antiprism with n-gon caps and 2n side triangles.
//
// - n: number of sides (n >= 3)
// - edge: target edge length
// - angle: additive twist offset in degrees relative to regular antiprism twist (180/n)
//          angle=0 => exact regular antiprism
// - height: explicit antiprism height (undef => solved from edge + twist)
// - height_scale: multiplier applied to chosen base height
function poly_antiprism(n=3, edge=1, angle=0, height=undef, height_scale=1) =
    let(
        _n_ok = assert(n >= 3, "poly_antiprism: n must be >= 3"),
        _e_ok = assert(edge > 0, "poly_antiprism: edge must be > 0"),
        _hs_ok = assert(height_scale > 0, "poly_antiprism: height_scale must be > 0"),
        h_base = is_undef(height) ? _ps_antiprism_height(n, edge, angle) : height,
        h = h_base * height_scale,
        _h_ok = assert(h >= 0, "poly_antiprism: height must be >= 0"),
        r = _ps_ngon_radius(n, edge),
        z0 = -h / 2,
        z1 = h / 2,
        theta = 180 / n + angle,
        bottom = _ps_ngon_ring(n, r, z0, 0),
        top = _ps_ngon_ring(n, r, z1, theta),
        verts = concat(bottom, top),
        faces_raw = concat(
            [[for (k = [0:1:n-1]) k]], // bottom
            [[for (k = [0:1:n-1]) n + k]],
            concat(
                // Standard antiprism side strip: two triangles per bottom edge.
                [for (k = [0:1:n-1]) [k, (k+1)%n, n+k]],
                [for (k = [0:1:n-1]) [(k+1)%n, n+((k+1)%n), n+k]]
            )
        ),
        faces = ps_orient_all_faces_outward(verts, faces_raw),
        ir = _ps_poly_ir(verts, faces),
        e_over_ir = edge / ir
    )
    make_poly(verts, faces, e_over_ir);
