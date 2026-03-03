// ---------------------------------------------------------------------------
// PolySymmetrica - Face-to-face poly attachment
// Attach p2 onto p1 by aligning selected faces and removing seam faces.
//
// v1 contract:
// - selected faces must be planar and have equal arity
// - supports optional cyclic correspondence shift via rotate_step
// - supports optional p2 scaling by selected-face average edge length
//
// API:
// poly_attach(
//   p1, p2,
//   f1=0, f2=0,
//   rotate_step=0,
//   scale_mode="fit_edge",   // "fit_edge" | "none"
//   eps=1e-8,
//   cleanup=true,
//   cleanup_eps=1e-8
// )

use <funcs.scad>
use <cleanup.scad>
use <validate.scad>

function _ps_rotate_cycle(list, step) =
    let(
        n = len(list),
        k0 = (n == 0) ? 0 : (round(step) % n),
        k = (k0 < 0) ? (k0 + n) : k0
    )
    [for (i = [0:1:n-1]) list[(i + k) % n]];

function _ps_face_center_raw(verts, f) =
    v_scale(v_sum([for (vi = f) verts[vi]]), 1 / len(f));

function _ps_face_frame_raw(verts, f) =
    let(
        c = _ps_face_center_raw(verts, f),
        ez = ps_face_frame_normal(verts, f),
        p0 = verts[f[0]] - c,
        ex_proj = p0 - ez * v_dot(p0, ez),
        p1 = verts[f[1]] - verts[f[0]],
        ex_fallback = p1 - ez * v_dot(p1, ez),
        ex = (norm(ex_proj) > 1e-12) ? v_norm(ex_proj) : v_norm(ex_fallback),
        ey = v_norm(v_cross(ez, ex))
    )
    [c, ex, ey, ez];

function _ps_face_avg_edge_len(verts, f) =
    let(
        n = len(f),
        lens = [for (k = [0:1:n-1]) norm(verts[f[(k+1)%n]] - verts[f[k]])]
    )
    (n == 0) ? 0 : sum(lens) / n;

function _ps_attach_map_point(p, frame1, frame2) =
    let(
        c1 = frame1[0], ex1 = frame1[1], ey1 = frame1[2], ez1 = frame1[3],
        c2 = frame2[0], ex2 = frame2[1], ey2 = frame2[2], ez2 = frame2[3],
        d = p - c2,
        x = v_dot(d, ex2),
        y = v_dot(d, ey2),
        z = v_dot(d, ez2)
    )
    c1 + x * ex1 + y * ey1 - z * ez1;

function _ps_drop_face_by_idx(faces, fi_drop) =
    [for (fi = [0:1:len(faces)-1]) if (fi != fi_drop) faces[fi]];

function poly_attach(
    p1, p2,
    f1=0, f2=0,
    rotate_step=0,
    scale_mode="fit_edge",
    eps=1e-8,
    cleanup=true,
    cleanup_eps=1e-8
) =
    let(
        _0 = assert(poly_valid(p1, "closed"), "attach: p1 must be closed-valid"),
        _1 = assert(poly_valid(p2, "closed"), "attach: p2 must be closed-valid"),
        _2 = assert(scale_mode == "fit_edge" || scale_mode == "none", "attach: scale_mode must be 'fit_edge' or 'none'"),

        v1 = poly_verts(p1),
        f1_all = ps_orient_all_faces_outward(v1, poly_faces(p1)),
        _3 = assert(f1 >= 0 && f1 < len(f1_all), "attach: f1 out of range"),
        face1 = f1_all[f1],

        v2_raw = poly_verts(p2),
        f2_all = ps_orient_all_faces_outward(v2_raw, poly_faces(p2)),
        _4 = assert(f2 >= 0 && f2 < len(f2_all), "attach: f2 out of range"),
        face2_base = f2_all[f2],
        face2 = _ps_rotate_cycle(face2_base, rotate_step),

        _5 = assert(len(face1) >= 3 && len(face2) >= 3, "attach: faces must have arity >= 3"),
        _6 = assert(len(face1) == len(face2), "attach: selected faces must have equal arity"),
        _7 = assert(_ps_face_planarity_err(v1, face1) <= eps, "attach: p1 selected face must be planar"),
        _8 = assert(_ps_face_planarity_err(v2_raw, face2) <= eps, "attach: p2 selected face must be planar"),

        len1 = _ps_face_avg_edge_len(v1, face1),
        len2 = _ps_face_avg_edge_len(v2_raw, face2),
        s2 = (scale_mode == "fit_edge") ? (len1 / len2) : 1,
        v2 = [for (v = v2_raw) v * s2],

        frame1 = _ps_face_frame_raw(v1, face1),
        frame2 = _ps_face_frame_raw(v2, face2),
        v2_t = [for (p = v2) _ps_attach_map_point(p, frame1, frame2)],

        faces1_keep = _ps_drop_face_by_idx(f1_all, f1),
        faces2_keep = _ps_drop_face_by_idx(f2_all, f2),
        off = len(v1),
        faces2_off = [for (f = faces2_keep) [for (vi = f) vi + off]],

        verts_join = concat(v1, v2_t),
        faces_join = concat(faces1_keep, faces2_off),
        p_join = [verts_join, faces_join, poly_e_over_ir(p1)],

        // Seam-merge always required for attachment.
        p_clean = poly_cleanup(
            p_join,
            eps=cleanup_eps,
            fix_winding=true,
            drop_degenerate=cleanup,
            triangulate_nonplanar=false,
            merge_vertices=true,
            remove_unreferenced=true
        ),
        _9 = assert(poly_valid(p_clean, "closed"), "attach: attached poly must be closed-valid")
    )
    p_clean;
