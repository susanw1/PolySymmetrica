// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Structural mesh cleanup helpers.
//
// poly_cleanup(poly, eps=1e-8, fix_winding=true, drop_degenerate=true,
//              triangulate_nonplanar=false, merge_vertices=true,
//              remove_unreferenced=true, return_report=false)
//
// Notes:
// - This is a structural normalizer, not a full topology repair tool.
// - It does not resolve self-intersections or non-manifold surgery.
// - It is useful after parameterized transforms that can produce coincident
//   or near-coincident indices/faces.

use <funcs.scad>

function _ps_identity_map(n) =
    [for (i = [0:1:n-1]) i];

function _ps_face_strip_adjacent_dups(f) =
    let(n = len(f))
    (n == 0) ? [] :
    [for (i = [0:1:n-1]) if (f[i] != f[(i-1+n)%n]) f[i]];

function _ps_face_trim_closing_dup(f) =
    (len(f) >= 2 && f[0] == f[len(f)-1]) ? [for (i = [0:1:len(f)-2]) f[i]] : f;

function _ps_face_clean_cycle(f) =
    _ps_face_trim_closing_dup(_ps_face_strip_adjacent_dups(f));

function _ps_faces_clean_cycles(faces) =
    [for (f = faces) _ps_face_clean_cycle(f)];

function _ps_distinct_count(list) =
    len([for (i = [0:1:len(list)-1]) if (_ps_index_of(list, list[i]) == i) 1]);

function _ps_face_area_mag(verts, f) =
    (len(f) < 3) ? 0 :
    sum([
        for (i = [1:1:len(f)-2])
            norm(v_cross(verts[f[i]] - verts[f[0]], verts[f[i+1]] - verts[f[0]])) / 2
    ]);

function _ps_face_planarity_err(verts, f, eps=1e-12) =
    (len(f) < 3) ? 0 :
    let(
        n_raw = ps_face_normal(verts, f),
        n_len = norm(n_raw),
        n = (n_len <= eps) ? [0,0,1] : (n_raw / n_len),
        d = v_dot(n, verts[f[0]]),
        errs = [for (vi = f) abs(v_dot(n, verts[vi]) - d)]
    )
    (len(errs) == 0) ? 0 : max(errs);

function _ps_faces_max_planarity_err(verts, faces, eps=1e-12) =
    (len(faces) == 0) ? 0 : max([for (f = faces) _ps_face_planarity_err(verts, f, eps)]);

function _ps_face_is_degenerate(verts, f, eps) =
    (len(f) < 3) ||
    (_ps_distinct_count(f) < 3) ||
    (_ps_distinct_count(f) != len(f)) ||
    (_ps_face_area_mag(verts, f) <= eps);

function _ps_faces_drop_degenerate(verts, faces, eps) =
    [for (f = faces) if (!_ps_face_is_degenerate(verts, f, eps)) f];

function _ps_face_triangulate_fan(f) =
    (len(f) < 3) ? [] : [for (i = [1:1:len(f)-2]) [f[0], f[i], f[i+1]]];

function _ps_face_is_planar(verts, f, eps) =
    _ps_face_planarity_err(verts, f, eps) <= eps;

function _ps_faces_triangulate_nonplanar(verts, faces, eps) =
    [
        for (f = faces)
            for (tri = (_ps_face_is_planar(verts, f, eps) ? [f] : _ps_face_triangulate_fan(f)))
                tri
    ];

function _ps_old_to_rep_by_eps(verts, eps) =
    [
        for (i = [0:1:len(verts)-1])
            let(idxs = [for (j = [0:1:i-1]) if (ps_point_eq(verts[j], verts[i], eps)) j])
            (len(idxs) == 0) ? i : idxs[0]
    ];

function _ps_merge_vertices_eps(verts, eps) =
    let(
        old_to_rep = _ps_old_to_rep_by_eps(verts, eps),
        rep_idxs = [for (i = [0:1:len(verts)-1]) if (old_to_rep[i] == i) i],
        old_to_new = [for (i = [0:1:len(verts)-1]) _ps_index_of(rep_idxs, old_to_rep[i])],
        verts_new = [for (ri = rep_idxs) verts[ri]]
    )
    [verts_new, old_to_new];

function _ps_faces_remap(faces, old_to_new) =
    [for (f = faces) [for (vi = f) old_to_new[vi]]];

function _ps_compact_unreferenced(verts, faces) =
    let(
        used = [for (f = faces) for (vi = f) vi],
        used_unique = [for (i = [0:1:len(used)-1]) if (_ps_index_of(used, used[i]) == i) used[i]],
        old_to_new = [for (i = [0:1:len(verts)-1]) _ps_index_of(used_unique, i)],
        verts_new = [for (oi = used_unique) verts[oi]],
        faces_new = _ps_faces_remap(faces, old_to_new)
    )
    [verts_new, faces_new];

// Returns cleaned poly, or [cleaned_poly, report] when return_report=true.
// report tuple:
// [faces_in, faces_out, faces_dropped, faces_triangulated,
//  verts_in, verts_out, verts_merged, verts_unreferenced_removed,
//  max_planarity_before, max_planarity_after]
function poly_cleanup(
    poly,
    eps=1e-8,
    fix_winding=true,
    drop_degenerate=true,
    triangulate_nonplanar=false,
    merge_vertices=true,
    remove_unreferenced=true,
    return_report=false
) =
    let(
        verts0 = poly_verts(poly),
        faces0 = poly_faces(poly),

        faces1 = _ps_faces_clean_cycles(faces0),

        merged = merge_vertices ? _ps_merge_vertices_eps(verts0, eps) : [verts0, _ps_identity_map(len(verts0))],
        verts2 = merged[0],
        map2 = merged[1],
        faces2 = _ps_faces_clean_cycles(_ps_faces_remap(faces1, map2)),

        faces3 = drop_degenerate ? _ps_faces_drop_degenerate(verts2, faces2, eps) : faces2,
        max_planarity_before = _ps_faces_max_planarity_err(verts2, faces3, eps),

        faces4 = triangulate_nonplanar ? _ps_faces_triangulate_nonplanar(verts2, faces3, eps) : faces3,
        faces5 = drop_degenerate ? _ps_faces_drop_degenerate(verts2, faces4, eps) : faces4,

        compacted = remove_unreferenced ? _ps_compact_unreferenced(verts2, faces5) : [verts2, faces5],
        verts6 = compacted[0],
        faces6 = compacted[1],
        max_planarity_after = _ps_faces_max_planarity_err(verts6, faces6, eps),

        // make_poly recomputes e_over_ir from current geometry.
        p0 = make_poly(verts6, faces6),
        p = fix_winding ? poly_fix_winding(p0) : p0,

        report = [
            len(faces0),
            len(faces6),
            len(faces0) - len(faces6),
            len(faces4) - len(faces3),
            len(verts0),
            len(verts6),
            len(verts0) - len(verts2),
            len(verts2) - len(verts6),
            max_planarity_before,
            max_planarity_after
        ]
    )
    return_report ? [p, report] : p;
