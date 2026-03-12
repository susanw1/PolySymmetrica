// ---------------------------------------------------------------------------
// PolySymmetrica - Construction helpers
// Topological construction primitives for open/capped mesh workflows.
//
// Initial API:
// - poly_delete_faces(poly, fids, cap=false, cleanup=true, cleanup_eps=1e-8)
// - poly_boundary_loops(poly)
// - poly_cap_loops(poly, loops=undef, cleanup=true, cleanup_eps=1e-8)
// - poly_slice(poly, plane_pt, plane_n, keep="above", cap=true, cleanup=true, cleanup_eps=1e-8)
//
// Notes:
// - These are explicit topology tools, intended as a substrate for later
//   slice/delete-vertex/delete-edge convenience wrappers.
// - `poly_delete_faces(..., cap=true)` deletes the selected faces, recovers
//   boundary loops, then adds caps over those loops.
// - Cap winding is normalized through `poly_cleanup(..., fix_winding=true)`,
//   so callers do not need to determine the boundary-loop orientation.

use <funcs.scad>
use <cleanup.scad>
use <transform.scad>
use <validate.scad>

function _ps_face_drop_by_idx_list(faces, idxs) =
    [for (fi = [0:1:len(faces)-1]) if (!_ps_list_contains(idxs, fi)) faces[fi]];

function _ps_faces_directed_edges(faces) =
    [
        for (f = faces)
            for (i = [0:1:len(f)-1])
                [f[i], f[(i+1)%len(f)]]
    ];

// Directed boundary edges: edges that are not matched by an opposite-directed
// edge in any other face.
function _ps_boundary_edges_directed(faces) =
    let(dirs = _ps_faces_directed_edges(faces))
    [
        for (i = [0:1:len(dirs)-1])
            let(
                e = dirs[i],
                hits = [
                    for (fi = [0:1:len(faces)-1])
                        if (ps_face_has_edge(faces[fi], e[0], e[1])) fi
                ]
            )
            if (len(hits) == 1) e
    ];

function _ps_find_boundary_edge_from(edges, v) =
    let(idxs = [for (i = [0:1:len(edges)-1]) if (edges[i][0] == v) i])
    (len(idxs) == 0) ? -1 : idxs[0];

function _ps_remove_edge_idx(edges, idx) =
    [for (i = [0:1:len(edges)-1]) if (i != idx) edges[i]];

function _ps_boundary_loop_from_edges(edges, start_edge, acc=[]) =
    let(
        e = edges[start_edge],
        rest = _ps_remove_edge_idx(edges, start_edge),
        acc2 = concat(acc, [e[0]]),
        next_v = e[1],
        done = (len(acc2) > 1 && next_v == acc2[0]),
        next_idx = done ? -1 : _ps_find_boundary_edge_from(rest, next_v)
    )
    done ? [acc2, rest]
         : (next_idx < 0)
            ? [acc2, rest]
            : _ps_boundary_loop_from_edges(rest, next_idx, acc2);

function _ps_boundary_loops_from_edges(edges, acc=[]) =
    (len(edges) == 0)
        ? acc
        : let(
            step = _ps_boundary_loop_from_edges(edges, 0, []),
            loop = step[0],
            rest = step[1]
        )
        _ps_boundary_loops_from_edges(rest, concat(acc, [loop]));

// Boundary vertex loops for an open mesh.
// Returns loops as vertex-index cycles, without repeating the first vertex.
function poly_boundary_loops(poly) =
    let(
        faces = poly_faces(poly),
        edges = _ps_boundary_edges_directed(faces),
        loops = _ps_boundary_loops_from_edges(edges)
    )
    [for (loop = loops) if (len(loop) >= 3) loop];

// Add caps for the provided boundary loops (or all boundary loops if omitted).
function poly_cap_loops(poly, loops=undef, cleanup=true, cleanup_eps=1e-8) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        loops_eff = is_undef(loops) ? poly_boundary_loops(poly) : loops,
        faces_new = concat(faces, loops_eff),
        p0 = [verts, faces_new, poly_e_over_ir(poly)]
    )
    cleanup
        ? poly_cleanup(
            p0,
            eps=cleanup_eps,
            fix_winding=true,
            drop_degenerate=true,
            triangulate_nonplanar=false,
            merge_vertices=true,
            remove_unreferenced=true
        )
        : p0;

// Delete faces by index. If cap=true, recover and cap the resulting openings.
function poly_delete_faces(poly, fids, cap=false, cleanup=true, cleanup_eps=1e-8) =
    let(
        fids_raw = is_list(fids) ? fids : [fids],
        fids_int = [for (fi = fids_raw) round(fi)],
        faces0 = poly_faces(poly),
        bad = [for (fi = fids_int) if (fi < 0 || fi >= len(faces0)) fi],
        _0 = assert(len(bad) == 0, "delete_faces: face index out of range"),
        faces1 = _ps_face_drop_by_idx_list(faces0, fids_int),
        p1 = [poly_verts(poly), faces1, poly_e_over_ir(poly)]
    )
    cap
        ? poly_cap_loops(p1, cleanup=cleanup, cleanup_eps=cleanup_eps)
        : (cleanup
            ? poly_cleanup(
                p1,
                eps=cleanup_eps,
                fix_winding=false,
                drop_degenerate=true,
                triangulate_nonplanar=false,
                merge_vertices=true,
                remove_unreferenced=true
            )
            : p1);

function _ps_plane_signed_dist(pt, plane_pt, plane_n) =
    v_dot(pt - plane_pt, plane_n);

function _ps_clip_edge_point(a, b, da, db, eps=1e-12) =
    let(den = da - db)
    (abs(den) <= eps) ? a : (a + (da / den) * (b - a));

function _ps_points_eq3(a, b, eps) =
    ps_point_eq(a, b, eps);

function _ps_clean_point_cycle(pts, eps=1e-8) =
    let(
        n = len(pts),
        no_adj = [
            for (i = [0:1:n-1])
                if (i == 0 || !_ps_points_eq3(pts[i], pts[i-1], eps))
                    pts[i]
        ],
        m = len(no_adj)
    )
    (m >= 2 && _ps_points_eq3(no_adj[0], no_adj[m-1], eps))
        ? [for (i = [0:1:m-2]) no_adj[i]]
        : no_adj;

function _ps_clip_face_to_halfspace(face_pts, plane_pt, plane_n, keep="above", eps=1e-8) =
    let(
        sign = (keep == "below") ? -1 : 1,
        n = len(face_pts),
        raw = [
            for (i = [0:1:n-1])
                let(
                    s = face_pts[i],
                    e = face_pts[(i+1)%n],
                    ds0 = sign * _ps_plane_signed_dist(s, plane_pt, plane_n),
                    ds1 = sign * _ps_plane_signed_dist(e, plane_pt, plane_n),
                    in0 = ds0 >= -eps,
                    in1 = ds1 >= -eps,
                    x = _ps_clip_edge_point(s, e, ds0, ds1, eps)
                )
                each (
                    in0 && in1 ? [e] :
                    in0 && !in1 ? [x] :
                    !in0 && in1 ? [x, e] :
                    []
                )
        ],
        cleaned = _ps_clean_point_cycle(raw, eps)
    )
    (len(cleaned) >= 3) ? cleaned : [];

function _ps_poly_from_face_points_preserve_scale(faces_pts_all, eps=1e-8) =
    let(
        all_pts = [for (fp = faces_pts_all) for (p = fp) p],
        uniq_verts = _ps_unique_points(all_pts, eps),
        faces_idx = [for (fp = faces_pts_all) _ps_face_points_to_indices(uniq_verts, fp, eps)],
        faces_simpl = [for (f = faces_idx) _ps_face_clean_cycle(f)],
        faces_keep = [for (f = faces_simpl) if (len(f) >= 3 && _ps_distinct_count(f) >= 3) f],
        _0 = assert(len(faces_keep) > 0, "slice: no faces remain after clipping"),
        faces_out = ps_orient_all_faces_outward(uniq_verts, faces_keep)
    )
    make_poly(uniq_verts, faces_out);

// Slice a closed polyhedron by a plane and keep one side.
// plane is given as point+normal; `keep` is "above" or "below" relative to the normal.
// With cap=true, the cut opening is recovered via boundary loops and capped.
function poly_slice(
    poly,
    plane_pt,
    plane_n,
    keep="above",
    cap=true,
    cleanup=true,
    cleanup_eps=1e-8
) =
    let(
        _0 = assert(poly_valid(poly, "closed"), "slice: source poly must be closed-valid"),
        _1 = assert(keep == "above" || keep == "below", "slice: keep must be 'above' or 'below'"),
        n_hat = v_norm(plane_n),
        _2 = assert(norm(n_hat) > 0, "slice: plane normal must be non-zero"),
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        clipped_faces = [
            for (f = faces)
                let(face_pts = [for (vi = f) verts[vi]])
                _ps_clip_face_to_halfspace(face_pts, plane_pt, n_hat, keep, cleanup_eps)
        ],
        keep_faces = [for (fp = clipped_faces) if (len(fp) >= 3) fp],
        p_open = _ps_poly_from_face_points_preserve_scale(keep_faces, cleanup_eps),
        p_norm = cleanup
            ? poly_cleanup(
                p_open,
                eps=cleanup_eps,
                fix_winding=false,
                drop_degenerate=true,
                triangulate_nonplanar=false,
                merge_vertices=true,
                remove_unreferenced=true
            )
            : p_open
    )
    cap
        ? poly_cap_loops(p_norm, cleanup=cleanup, cleanup_eps=cleanup_eps)
        : p_norm;
