// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral transform helpers
// Version: 0.1.0
// Copyright 2026 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>

// Generic site/cycle poly transform kernel.
// Keep this file minimal and reusable across future operators (snub, cantitruncate, etc.).

// index of point p in list (or -1)
function _ps_find_point(list, p, eps, i=0) =
    (i >= len(list)) ? -1 :
    (ps_point_eq(list[i], p, eps) ? i : _ps_find_point(list, p, eps, i+1));

// Build unique vertex list from a flat list of points
function _ps_unique_points(points, eps, acc=[], i=0) =
    (i >= len(points)) ? acc :
    let(p = points[i])
    (_ps_find_point(acc, p, eps) >= 0)
        ? _ps_unique_points(points, eps, acc, i+1)
        : _ps_unique_points(points, eps, concat(acc, [p]), i+1);

// Remap a face described by points -> indices in uniq[]
function _ps_face_points_to_indices(uniq, face_pts, eps) =
    [ for (p = face_pts) _ps_find_point(uniq, p, eps) ];

// Build a poly descriptor from face point lists (dedup + orient + unit-edge scale).
function _ps_poly_from_face_points(faces_pts_all, eps, len_eps=undef) =
    let(
        len_eps_eff = is_undef(len_eps) ? eps : len_eps,
        all_pts = [ for (fp = faces_pts_all) for (p = fp) p ],
        uniq_verts = _ps_unique_points(all_pts, len_eps_eff),
        faces_idx = [ for (fp = faces_pts_all) _ps_face_points_to_indices(uniq_verts, fp, len_eps_eff) ],
        faces_out = ps_orient_all_faces_outward(uniq_verts, faces_idx),
        edges_new = _ps_edges_from_faces(faces_out),
        e0 = edges_new[0],
        vA = uniq_verts[e0[0]],
        vB = uniq_verts[e0[1]],
        unit_e = norm(vB - vA),
        mid = (vA + vB) / 2,
        ir  = norm(mid),
        e_over_ir = unit_e / ir
    )
    make_poly(uniq_verts / unit_e, faces_out, e_over_ir);

// Build a poly from site-based face cycles.
// cycle entries: [0, v_idx] for original vertex, [1, site_idx] for site point.
function ps_poly_transform_from_sites(verts0, sites, site_points, face_cycles, eps=1e-8, len_eps=1e-6) =
    let(
        faces_pts_all = [
            for (cy = face_cycles)
                [ for (c = cy) (c[0] == 0) ? verts0[c[1]] : site_points[c[1]] ]
        ]
    )
    _ps_poly_from_face_points(faces_pts_all, eps, len_eps);
