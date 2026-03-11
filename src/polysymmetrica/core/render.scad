// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2026 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <placement.scad>
use <segments.scad>
use <validate.scad>

function _ps_render_cross2(a, b) = a[0]*b[1] - a[1]*b[0];

function _ps_render_face_pts2d(verts, f, eps=1e-9) =
    let(
        n = ps_face_frame_normal(verts, f, eps),
        c = ps_face_centroid(verts, f),
        a0 = verts[f[0]] - c,
        proj = a0 - n * v_dot(a0, n),
        proj_len = norm(proj),
        ex = (proj_len <= eps) ? [1,0,0] : proj / proj_len,
        ey = v_cross(n, ex)
    )
    [ for (vid = f)
        let(vdir = verts[vid] - c)
        [v_dot(vdir, ex), v_dot(vdir, ey)]
    ];

function _ps_render_face_is_concave(verts, f, eps=1e-9) =
    let(
        pts = _ps_render_face_pts2d(verts, f, eps),
        n = len(pts),
        turns = [
            for (i = [0:1:n-1])
                let(
                    p0 = pts[(i - 1 + n) % n],
                    p1 = pts[i],
                    p2 = pts[(i + 1) % n],
                    e0 = [p1[0] - p0[0], p1[1] - p0[1]],
                    e1 = [p2[0] - p1[0], p2[1] - p1[1]]
                )
                _ps_render_cross2(e0, e1)
        ],
        has_pos = max([for (t = turns) (t > eps) ? 1 : 0]) == 1,
        has_neg = max([for (t = turns) (t < -eps) ? 1 : 0]) == 1
    )
    has_pos && has_neg;

function _ps_render_face_needs_triangulation(verts, f, eps=1e-8) =
    (_ps_face_planarity_err(verts, f, eps) > eps) ||
    (_ps_face_self_intersections(verts, f, eps) > 0) ||
    _ps_render_face_is_concave(verts, f, eps);

function _ps_render_face_centroid_fan_tris3(verts, f) =
    let(
        n = len(f),
        c = v_scale(v_sum([for (vi = f) verts[vi]]), 1 / n)
    )
    (n < 3) ? [] :
    [ for (k = [0:1:n-1]) [verts[f[k]], verts[f[(k+1)%n]], c] ];

function _ps_render_face_tris3(verts, f, eps=1e-8) =
    // Self-intersecting faces use centroid fan to preserve original boundary
    // edges (keeps shell closed against adjacent faces).
    (_ps_face_self_intersections(verts, f, eps) > 0)
        ? _ps_render_face_centroid_fan_tris3(verts, f)
        : _ps_seg_face_tris3(f, verts, eps, "nonzero");

function _ps_render_find_point(list, p, eps, i=0) =
    (i >= len(list)) ? -1 :
    (ps_point_eq(list[i], p, eps) ? i : _ps_render_find_point(list, p, eps, i+1));

function _ps_render_unique_points(points, eps, acc=[], i=0) =
    (i >= len(points)) ? acc :
    let(p = points[i])
    (_ps_render_find_point(acc, p, eps) >= 0)
        ? _ps_render_unique_points(points, eps, acc, i+1)
        : _ps_render_unique_points(points, eps, concat(acc, [p]), i+1);

function _ps_render_face_indices(face_pts, uniq_pts, eps) =
    [ for (p = face_pts) _ps_render_find_point(uniq_pts, p, eps) ];

function _ps_render_tri_nondegenerate(verts, f, eps=1e-12) =
    let(
        a = verts[f[0]],
        b = verts[f[1]],
        c = verts[f[2]],
        n = v_cross(b - a, c - a)
    )
    norm(n) > eps;

function _ps_render_tri_mesh(verts, faces, eps=1e-8) =
    let(
        tris3 = [
            for (f = faces)
                for (t = _ps_render_face_tris3(verts, f, eps))
                    t
        ],
        raw_pts = [for (t = tris3) for (p = t) p],
        uniq_pts = _ps_render_unique_points(raw_pts, eps),
        tri_faces_raw = [for (t = tris3) _ps_render_face_indices(t, uniq_pts, eps)],
        tri_faces_valid = [
            for (f = tri_faces_raw)
                if (len(f) == 3 && f[0] >= 0 && f[1] >= 0 && f[2] >= 0 &&
                    f[0] != f[1] && f[1] != f[2] && f[0] != f[2] &&
                    _ps_render_tri_nondegenerate(uniq_pts, f, eps))
                    f
        ],
        tri_faces = ps_orient_all_faces_outward(uniq_pts, tri_faces_valid)
    )
    [uniq_pts, tri_faces];

// Return [points, faces] suitable for OpenSCAD polyhedron().
// triangulate:
// - "auto" (default): triangulate only when a face is non-planar/concave/self-intersecting
// - "always"/true: always triangulate
// - "raw"/false: pass original faces through unchanged
function ps_render_mesh(poly, inter_radius = 1, triangulate = "auto", eps=1e-8) =
    let(
        k = inter_radius * poly_e_over_ir(poly),
        verts = poly_verts(poly) * k,
        faces = poly_faces(poly),
        force_tri = (triangulate == true) || (triangulate == "always") || (triangulate == "triangulate"),
        force_raw = (triangulate == false) || (triangulate == "never") || (triangulate == "raw"),
        need_tri = force_tri ? true :
            force_raw ? false :
            (max([for (f = faces) _ps_render_face_needs_triangulation(verts, f, eps) ? 1 : 0]) == 1)
    )
    need_tri ? _ps_render_tri_mesh(verts, faces, eps) : [verts, faces];

module poly_render(poly, inter_radius = 1, triangulate = "auto", eps=1e-8) {
    let(
        mesh = ps_render_mesh(poly, inter_radius, triangulate, eps)
    )
    polyhedron(
        points = mesh[0],
        faces  = mesh[1]
    );
}


module poly_describe(poly, name = undef, detail = 1) {
    echo ("==========", is_undef(name)? "Polyhedron" : name, "============");
    verts = poly_verts(poly);
    faces = poly_faces(poly);

    echo ("#vertices: ", len(verts));
    echo ("#faces: ", len(faces));
    echo ("#edges: ", len(_ps_edges_from_faces(faces)));
    echo ("e_over_ir: ", poly_e_over_ir(poly));
    
    if (detail > 0) {
        place_on_faces(poly, 1) {
            let(
                face_verts = faces[$ps_face_idx],
                basic1 = str("face#", $ps_face_idx, ": (", $ps_vertex_count, " verts) vert_idx: ", face_verts),
                more2 = str(" poly_rad: ", $ps_face_midradius, " verts: ", [for (f = face_verts) verts[f]],
                    " lengths: ", [for (fi = [0:1:len(face_verts)-1]) norm(verts[face_verts[fi]]-verts[face_verts[(fi+1)%len(face_verts)]]) ]),
                max_plane_err = _ps_face_planarity_err(verts, face_verts),
                more3 = str(" max_plane_err: ", max_plane_err)
            )
            echo (str("  ", basic1, detail > 1 ? more2 : "", detail > 2 ? more3 : ""));
        }
        place_on_vertices(poly, 1) {
            let(
                basic1 = str("vertex#", $ps_vertex_idx, ": (", $ps_vertex_valence, " edges) neighbours_idx: ", $ps_vertex_neighbors_idx),
                more2 = str("vert_idx: ", verts[$ps_vertex_idx])
            )
            echo (str("  ", basic1, detail > 1? more2 : ""));
        }
        place_on_edges(poly, 1) {
            let(
                basic1 = str("edge#", $ps_edge_idx, ": len: ", $ps_edge_len, " verts_idx: ", $ps_edge_verts_idx, " faces_idx: ", $ps_edge_adj_faces_idx),
                more2 = str(" ")
            )
            echo (str("  ", basic1, detail > 1? more2 : ""));
        }
    }
    echo ("======================");
}
