// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <duals.scad>  // for faces_around_vertex helpers
use <transform.scad>
use <transform_util.scad>
use <solvers.scad>
use <validate.scad>

// --- internal helpers ---

// Site index for edge-point near a vertex.
function _ps_edge_site_index(edges, a, b, near_v) =
    let(
        ei = ps_find_edge_index(edges, a, b),
        e  = edges[ei]
    )
    2 * ei + ((near_v == e[0]) ? 0 : 1);

// True if face f contains directed edge v0->v1.
function _ps_face_has_dir(f, v0, v1) =
    let(pos = _ps_index_of(f, v0))
    (pos >= 0) ? (f[(pos+1)%len(f)] == v1) : false;

function _ps_all_near(vals, eps) =
    (len(vals) == 0) ? true : (max(vals) - min(vals) <= eps);

// True if poly has a single face size and single edge length (within eps).
function _ps_is_regular_base(poly, eps=1e-6) =
    let(
        verts = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        face_sizes = [for (f = faces0) len(f)],
        edge_lens = [for (e = edges) norm(verts[e[1]] - verts[e[0]])]
    )
    _ps_all_near(face_sizes, eps) && _ps_all_near(edge_lens, eps);

// Project edge points to face plane and order along face edge direction.
function _ps_project_edge_pts_for_face_edge(verts0, edges, edge_pts, n_f, p0, v0, v1) =
    let(
        ei = ps_find_edge_index(edges, v0, v1),
        e = edges[ei],
        p_a = edge_pts[ei][0],
        p_b = edge_pts[ei][1],
        p_near_v0 = (e[0] == v0) ? p_a : p_b,
        p_near_v1 = (e[0] == v1) ? p_a : p_b,
        proj0 = p_near_v0 + n_f * v_dot(n_f, (p0 - p_near_v0)),
        proj1 = p_near_v1 + n_f * v_dot(n_f, (p0 - p_near_v1)),
        e_vec = verts0[v1] - verts0[v0],
        e_dir = v_norm(e_vec - n_f * v_dot(n_f, e_vec)),
        flip = v_dot((proj1 - proj0), e_dir) < 0
    )
    [flip ? proj1 : proj0, flip ? proj0 : proj1];

// Intersect 2D lines n0·x=d0 and n1·x=d1.
function _ps_line2_intersect(n0, d0, n1, d1, eps=1e-12) =
    let(det = n0[0]*n1[1] - n0[1]*n1[0])
    (abs(det) < eps) ? undef
  : [
        (d0*n1[1] - n0[1]*d1) / det,
        (n0[0]*d1 - d0*n1[0]) / det
    ];

function _ps_rot2d(p, ang) =
    let(c = cos(ang), s = sin(ang))
    [p[0]*c - p[1]*s, p[0]*s + p[1]*c];

// Inset polygon vertices using bisector-plane intersection lines per edge.
// d_f shifts the face plane; d_e shifts the edge-bisector planes along their normals.
// Inset polygon vertices using bisector-plane intersection lines per edge.
// d_e may be a single value or a per-edge list (length n). Per-edge values
// allow mixed-family cantitruncation to bias individual edge planes.
function _ps_face_inset_bisector_2d(f, fi, d_f, d_e, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0) =
    let(
        n = len(f),
        p0 = center - n_f * d_f,
        d_e_list = (is_list(d_e) ? d_e : [for (_ = [0:1:n-1]) d_e]),
        lines = [
            for (k = [0:1:n-1])
                let(
                    v0 = f[k],
                    v1 = f[(k+1)%n],
                    ei = ps_find_edge_index(edges, v0, v1),
                    adj = edge_faces[ei],
                    f_adj = (len(adj) < 2) ? fi : ((adj[0] == fi) ? adj[1] : adj[0]),
                    n_adj = face_n[f_adj],
                    n_edge_raw = n_f + n_adj,
                    n_edge = (v_len(n_edge_raw) < 1e-8) ? n_f : v_norm(n_edge_raw),
                    n2d = [v_dot(n_edge, ex), v_dot(n_edge, ey)],
                    n2d_len = v_len(n2d),
                    p1_2d = pts2d[(k+1)%n],
                    e2d = [p1_2d[0]-pts2d[k][0], p1_2d[1]-pts2d[k][1]],
                    n2d_fb = v_norm([e2d[1], -e2d[0]]),
                    n2d_use = (n2d_len < 1e-8) ? n2d_fb : (n2d / n2d_len),
                    d2 = v_dot(n_edge, verts0[v0]) - v_dot(n_edge, p0) + d_e_list[k],
                    d2_use = (n2d_len < 1e-8) ? (v_dot(n2d_use, pts2d[k]) + d_f) : (d2 / n2d_len)
                )
                [n2d_use, d2_use]
        ],
        inset2d_raw = [
            for (k = [0:1:n-1])
                _ps_line2_intersect(
                    lines[(k-1+n)%n][0], lines[(k-1+n)%n][1],
                    lines[k][0], lines[k][1]
                )
        ],
        inset2d = [
            for (k = [0:1:n-1])
                is_undef(inset2d_raw[k]) ? pts2d[k] : inset2d_raw[k]
        ]
    )
    inset2d;

// 2D polygon area (absolute), pts may be undef.
function _ps_poly_area_abs_2d(pts2d) =
    (is_undef(pts2d) || len(pts2d) < 3) ? 0
  : let(
        n = len(pts2d),
        a = sum([
            for (i = [0:1:n-1])
                let(
                    p0 = pts2d[i],
                    p1 = pts2d[(i+1)%n]
                )
                (is_undef(p0) || is_undef(p1)) ? 0 : (p0[0]*p1[1] - p1[0]*p0[1])
        ])
    )
    abs(a) / 2;

function _ps_poly_area_abs_2d_valid(pts2d) =
    let(
        n = is_undef(pts2d) ? 0 : len(pts2d),
        ok = (n >= 3) && (sum([for (p = pts2d) is_undef(p) ? 1 : 0]) == 0)
    )
    ok ? _ps_poly_area_abs_2d(pts2d) : undef;

function _ps_face_inset_area_2d(f, fi, d_f, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0) =
    let(poly2d = _ps_face_inset_bisector_2d(f, fi, d_f, 0, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0))
    _ps_poly_area_abs_2d_valid(poly2d);

function _ps_expand_collapse_d(f, fi, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0, prev_d, prev_area, d, iter, max_iter) =
    let(area = _ps_face_inset_area_2d(f, fi, d, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0))
    (iter >= max_iter) ? [d, area] :
    (is_undef(area) || area <= 0) ? [d, area] :
    (area > prev_area) ? [prev_d, prev_area] :
    _ps_expand_collapse_d(f, fi, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0, d, area, d * 2, iter + 1, max_iter);

// Find d_f at which bisector-based inset polygon collapses (min area).
function _ps_face_bisector_collapse_d(f, fi, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0, steps=40) =
    let(
        area0 = _ps_poly_area_abs_2d(pts2d),
        r_max = max([for (p = pts2d) norm(p)]),
        d0 = (r_max > 0) ? r_max : 1,
        area_p0 = _ps_face_inset_area_2d(f, fi, d0, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0),
        area_n0 = _ps_face_inset_area_2d(f, fi, -d0, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0),
        a_p0 = is_undef(area_p0) ? 1e30 : area_p0,
        a_n0 = is_undef(area_n0) ? 1e30 : area_n0,
        any_decrease = (a_n0 < area0) || (a_p0 < area0),
        dir = (a_n0 < area0) ? -1 : (a_p0 < area0) ? 1 : -1,
        d_start = dir * d0,
        dmax_pair = (area0 <= 0) ? [0, area0]
                  : !any_decrease ? [d0, area0]
                  : _ps_expand_collapse_d(f, fi, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0, 0, area0, d_start, 0, 20),
        d_max = abs(dmax_pair[0]),
        ds_dir = [ for (i = [0:1:steps]) dir * d_max * i / steps ],
        areas_raw = [
            for (i = [0:1:len(ds_dir)-1])
                (ds_dir[i] == 0) ? area0 : _ps_face_inset_area_2d(f, fi, ds_dir[i], center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0)
        ],
        areas = [ for (a = areas_raw) is_undef(a) ? 1e30 : a ],
        min_i = _ps_index_of_min(areas),
        all_bad = (min(areas) >= 1e29)
    )
    all_bad ? d0 : abs(ds_dir[min_i]);

function _ps_index_of_min(list) =
    let(
        n = len(list),
        vals = [ for (i = [0:1:n-1]) [list[i], i] ],
        min_val = min([for (v = vals) v[0]]),
        idxs = [for (v = vals) if (v[0] == min_val) v[1]]
    )
    (len(idxs) == 0) ? 0 : idxs[0];

// --- main truncation ---

function _ps_truncate_norm_to_t(poly, c) =
    _ps_truncate_default_t(poly) * c;

function poly_truncate(poly, t=undef, c=undef, eps = 1e-8) =
    let(
        t_eff = !is_undef(t)
            ? t
            : (!is_undef(c) ? _ps_truncate_norm_to_t(poly, c) : _ps_truncate_default_t(poly))
    )
    // should this assert? Or should we allow, or silently call rectify?
    assert(t_eff != 0.5, "'t' cannot be 0.5 as this produces degenerate vertices - use poly_rectify() instead")
    (t_eff == 0) 
        ? poly
        :
    let(
        base = _ps_poly_base(poly),
        verts = base[0],
        faces = base[1],
        edges = base[2],

        // edge points: aligned with edges[ei]=[a,b]
        edge_pts = _ps_edge_points(verts, edges, t_eff),

        edge_faces = ps_edge_faces_table(faces, edges),

        // sites: two per edge (near each endpoint)
        sites = [
            for (ei = [0:1:len(edges)-1])
                each [[ei, edges[ei][0]], [ei, edges[ei][1]]]
        ],
        site_points = [
            for (ei = [0:1:len(edges)-1])
                each [edge_pts[ei][0], edge_pts[ei][1]]
        ],

        // truncated original faces -> 2n-gons
        face_cycles = [
            for (f = faces)
                let(n = len(f))
                [
                    for (k = [0:1:n-1])
                        let(
                            v      = f[k],
                            v_next = f[(k+1)%n],
                            v_prev = f[(k-1+n)%n],
                            idx_prev = _ps_edge_site_index(edges, v_prev, v, v),
                            idx_next = _ps_edge_site_index(edges, v, v_next, v)
                        )
                        each [[1, idx_prev], [1, idx_next]]
                ]
        ],

        // vertex faces: one per original vertex
        vert_cycles = [
            for (vi = [0:1:len(verts)-1])
                let(
                    fc = faces_around_vertex(poly, vi, edges, edge_faces),
                    neigh = [
                        for (idx = [0:1:len(fc)-1])
                            let(
                                f = faces[ fc[idx] ],
                                m = len(f),
                                pos = _ps_index_of(f, vi),
                                v_next = f[(pos+1)%m]
                            )
                            v_next
                    ]
                )
                [
                    for (vn = neigh)
                        [1, _ps_edge_site_index(edges, vi, vn, vi)]
                ]
        ],

        cycles_all = concat(face_cycles, vert_cycles)
    )
    ps_poly_transform_from_sites(verts, sites, site_points, cycles_all, eps, eps);

// Rectification: replace each vertex with the midpoint of each incident edge.
function poly_rectify(poly) =
    let(
        base = _ps_poly_base(poly),
        verts = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        poly0 = base[5],

        edge_mid = [
            for (e = edges)
                (verts[e[0]] + verts[e[1]]) / 2
        ],

        // Faces corresponding to original faces: walk edges in face order.
        face_faces = [
            for (f = faces0)
                let(n = len(f))
                [
                    for (k = [0:1:n-1])
                        let(
                            a = f[k],
                            b = f[(k+1)%n],
                            ei = ps_find_edge_index(edges, a, b)
                        )
                        ei
                ]
        ],

        // Faces corresponding to original vertices: cycle around the vertex.
        vert_faces = [
            for (vi = [0:1:len(verts)-1])
                let(
                    fc = faces_around_vertex(poly0, vi, edges, edge_faces),
                    neigh = [
                        for (idx = [0:1:len(fc)-1])
                            let(
                                f = faces0[fc[idx]],
                                m = len(f),
                                pos = _ps_index_of(f, vi),
                                v_next = f[(pos+1)%m]
                            )
                            v_next
                    ]
                )
                [
                    for (vn = neigh)
                        let(ei = ps_find_edge_index(edges, vi, vn))
                        ei
                ]
        ],

        faces_idx = concat(face_faces, vert_faces),
        cycles_all = [
            for (f = faces_idx)
                [ for (ei = f) [1, ei] ]
        ]
    )
    ps_poly_transform_from_sites(verts, [for (i = [0:1:len(edge_mid)-1]) [i]], edge_mid, cycles_all);

// Chamfer: face faces + edge faces (no vertex faces).
// t is a signed face-plane offset as a fraction of each face's collapse distance.
// Positive t = inward chamfer; negative t = anti-chamfer.
function poly_chamfer(poly, t=undef, c=undef, eps = 1e-8, len_eps = 1e-6) =
    let(
        t_eff = !is_undef(t)
            ? t
            : (!is_undef(c) ? _ps_truncate_norm_to_t(poly, c) : _ps_truncate_default_t(poly)),
        _t_ok = assert(abs(t_eff) != 1, "poly_chamfer: |t| must not be 1 (t=±1 collapses faces; use cleanup/hyper mode)"),
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],

        face_offsets = _ps_face_offsets(faces0),

        face_pts3d = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    n_f = face_n[fi],
                    center = poly_face_center(poly, fi, 1),
                    ex = poly_face_ex(poly, fi, 1),
                    ey = poly_face_ey(poly, fi, 1),
                    pts2d = [
                        for (k = [0:1:n-1])
                            let(p = verts0[f[k]] - center)
                                [v_dot(p, ex), v_dot(p, ey)]
                    ],
                    edge_lens = [
                        for (k = [0:1:n-1])
                            let(p0 = pts2d[k], p1 = pts2d[(k+1)%n])
                                norm([p1[0]-p0[0], p1[1]-p0[1]])
                    ],
                    face_collapse = abs(_ps_face_bisector_collapse_d(f, fi, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0)),
                    d_f0 = -t_eff * face_collapse,
                    d_f = d_f0,
                    inset2d = _ps_face_inset_bisector_2d(f, fi, d_f, 0, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0),
                    p0 = center - n_f * d_f
                )
                [
                    for (k = [0:1:n-1])
                        p0 + ex * inset2d[k][0] + ey * inset2d[k][1]
                ]
        ],

        sites = [
            for (fi = [0:1:len(faces0)-1])
                let(f = faces0[fi], n = len(f))
                for (k = [0:1:n-1])
                    [fi, f[k]]
        ],

        site_points = [
            for (fi = [0:1:len(faces0)-1])
                for (p = face_pts3d[fi])
                    p
        ],
        face_cycles = [
            for (fi = [0:1:len(faces0)-1])
                let(f = faces0[fi], n = len(f))
                [ for (k = [0:1:n-1]) [1, face_offsets[fi] + k] ]
        ],

        edge_cycles = [
            for (ei = [0:1:len(edges)-1])
                let(
                    e = edges[ei],
                    v0 = e[0],
                    v1 = e[1],
                    fpair = edge_faces[ei],
                    f0 = fpair[0],
                    f1 = fpair[1],
                    f_dir = _ps_face_has_dir(faces0[f0], v0, v1) ? f0 : f1,
                    f_opp = (f_dir == f0) ? f1 : f0,
                    k_dir_v0 = _ps_index_of(faces0[f_dir], v0),
                    k_dir_v1 = _ps_index_of(faces0[f_dir], v1),
                    k_opp_v1 = _ps_index_of(faces0[f_opp], v1),
                    k_opp_v0 = _ps_index_of(faces0[f_opp], v0),
                    s_dir_v0 = face_offsets[f_dir] + k_dir_v0,
                    s_dir_v1 = face_offsets[f_dir] + k_dir_v1,
                    s_opp_v1 = face_offsets[f_opp] + k_opp_v1,
                    s_opp_v0 = face_offsets[f_opp] + k_opp_v0
                )
                [
                    [0, v0],
                    [1, s_dir_v0],
                    [1, s_dir_v1],
                    [0, v1],
                    [1, s_opp_v1],
                    [1, s_opp_v0]
                ]
        ],

        cycles_all = concat(face_cycles, edge_cycles)
    )
    ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps);

// Cantellation/expansion:
// - face faces remain n-gons (one point per original vertex)
// - edge faces are quads (rectangles/squares)
// - vertex faces are valence-gons
// df offsets faces along their normals; edge/vertex faces are derived from those offsets.
function _ps_face_offset_pts(verts0, faces0, face_n, df) =
    [
        for (fi = [0:1:len(faces0)-1])
            let(
                f = faces0[fi],
                n = len(f),
                n_f = face_n[fi]
            )
            [
                for (k = [0:1:n-1])
                    let(v = f[k])
                    verts0[v] + df * n_f
            ]
    ];

function _ps_cantellate_df_from_c(poly, c, df_max=undef, steps=16, family_edge_idx=0) =
    let(
        verts = poly_verts(poly),
        edges = _ps_edges_from_faces(poly_faces(poly)),
        e0 = edges[0],
        edge_len = norm(verts[e0[1]] - verts[e0[0]]),
        ir = edge_len / poly_e_over_ir(poly),
        df_max_eff = is_undef(df_max) ? (2 * ir) : df_max,
        df_mid = cantellate_square_df(poly, 0, df_max_eff, steps, family_edge_idx),
        df = (c <= 0.5)
            ? (2 * c * df_mid)
            : (df_mid + (c - 0.5) * 2 * (df_max_eff - df_mid))
    )
    df;

function poly_cantellate(poly, df=undef, c=undef, df_max=undef, steps=16, family_edge_idx=0, eps = 1e-8, len_eps = 1e-6) =
    let(
        df_eff = !is_undef(df)
            ? df
            : (!is_undef(c) ? _ps_cantellate_df_from_c(poly, c, df_max, steps, family_edge_idx)
                            : _ps_cantellate_df_from_c(poly, 0.5, df_max, steps, family_edge_idx))
    )
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        // offset face corners: one point per (face, vertex) incidence
        face_pts = _ps_face_offset_pts(verts0, faces0, face_n, df_eff),
        face_offsets = _ps_face_offsets(faces0),
        sites = [
            for (fi = [0:1:len(faces0)-1])
                for (v = faces0[fi])
                    [fi, v]
        ],
        site_points = [
            for (fi = [0:1:len(faces0)-1])
                for (p = face_pts[fi])
                    p
        ],
        face_cycles = [
            for (fi = [0:1:len(faces0)-1])
                let(n = len(faces0[fi]))
                [ for (k = [0:1:n-1]) [1, face_offsets[fi] + k] ]
        ],
        edge_cycles = [
            for (ei = [0:1:len(edges)-1])
                let(
                    e = edges[ei],
                    fpair = edge_faces[ei],
                    f0 = fpair[0],
                    f1 = fpair[1],
                    v0 = e[0],
                    v1 = e[1],
                    pos0 = _ps_index_of(faces0[f0], v0),
                    pos1 = _ps_index_of(faces0[f0], v1),
                    pos2 = _ps_index_of(faces0[f1], v1),
                    pos3 = _ps_index_of(faces0[f1], v0)
                )
                [
                    [1, face_offsets[f0] + pos0],
                    [1, face_offsets[f0] + pos1],
                    [1, face_offsets[f1] + pos2],
                    [1, face_offsets[f1] + pos3]
                ]
        ],
        vert_cycles = [
            for (vi = [0:1:len(verts0)-1])
                let(fc = faces_around_vertex(poly0, vi, edges, edge_faces))
                [
                    for (fi = fc)
                        let(pos = _ps_index_of(faces0[fi], vi))
                            [1, face_offsets[fi] + pos]
                ]
        ],
        cycles_all = concat(face_cycles, edge_cycles, vert_cycles)
    )
    ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps);

// Measure how square an edge face is (edge length spread).
function _ps_face_edge_spread(verts, face) =
    let(
        n = len(face),
        ls = (n < 2) ? [] : [for (i = [0:1:n-1]) norm(verts[face[i]] - verts[face[(i+1)%n]])]
    )
    (n == 4) ? (max(ls) - min(ls)) : undef;

// Solve df so that the chosen edge-face family is as square as possible.
function cantellate_square_df(poly, df_min, df_max, steps=40, family_edge_idx=0, eps=1e-9) =
    let(
        faces0 = ps_orient_all_faces_outward(poly_verts(poly), poly_faces(poly)),
        edges0 = _ps_edges_from_faces(faces0),
        edge_faces0 = ps_edge_faces_table(faces0, edges0),
        fpair = edge_faces0[family_edge_idx],
        key = (len(fpair) == 2) ? ((len(faces0[fpair[0]]) <= len(faces0[fpair[1]]))
            ? [len(faces0[fpair[0]]), len(faces0[fpair[1]])]
            : [len(faces0[fpair[1]]), len(faces0[fpair[0]])]) : [0,0],
        family_edges = [
            for (ei = [0:1:len(edge_faces0)-1])
                let(
                    fpair_i = edge_faces0[ei],
                    key_i = (len(fpair_i) == 2) ? ((len(faces0[fpair_i[0]]) <= len(faces0[fpair_i[1]]))
                        ? [len(faces0[fpair_i[0]]), len(faces0[fpair_i[1]])]
                        : [len(faces0[fpair_i[1]]), len(faces0[fpair_i[0]])]) : [0,0]
                )
                if (key_i == key) ei
        ],
        df_vals = [for (i = [0:1:steps]) df_min + (df_max - df_min) * i / steps],
        cands = [
            for (df = df_vals)
                let(
                    q = poly_cantellate(poly, df),
                    faces_q = poly_faces(q),
                    face_count = len(faces0),
                    errs = [for (ei = family_edges) _ps_face_edge_spread(poly_verts(q), faces_q[face_count + ei])],
                    ok = len([for (e = errs) if (!is_undef(e)) 1]) == len(errs)
                )
                if (ok) [df, sum(errs)]
        ],
        _ = assert(len(cands) > 0, "cantellate_square_df: no valid candidates"),
        errs = [for (c = cands) c[1]],
        e_min = min(errs),
        idx = [for (i = [0:1:len(errs)-1]) if (abs(errs[i] - e_min) <= eps) i][0]
    )
    cands[idx][0];

// Normalized cantellation: map c in [0,1] to df in [0, df_max],
// with c=0.5 hitting the computed square-edge df.
function poly_cantellate_norm(poly, c, df_max=undef, steps=16, family_edge_idx=0, eps=1e-8, len_eps=1e-6) =
    let(df = _ps_cantellate_df_from_c(poly, c, df_max, steps, family_edge_idx))
    poly_cantellate(poly, df, undef, df_max, steps, family_edge_idx, eps, len_eps);

// Snub (chiral): twist cantellated face points within each face plane and triangulate edge cycles.
function _ps_snub_edge_spread(poly, c_eff, angle, handedness=1) =
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        df = _ps_cantellate_df_from_c(poly, c_eff),
        ei = 0,
        e = edges[ei],
        fpair = edge_faces[ei],
        f0 = fpair[0],
        f1 = fpair[1],
        v0 = e[0],
        // third face around v0
        fc = faces_around_vertex(poly0, v0, edges, edge_faces),
        f_other = [for (fi = fc) if (fi != f0 && fi != f1) fi][0],
        s_f0 = _ps_snub_face_point(verts0, faces0, face_n, df, f0, v0, handedness, angle, poly0),
        s_f1 = _ps_snub_face_point(verts0, faces0, face_n, df, f1, v0, handedness, angle, poly0),
        s_o = _ps_snub_face_point(verts0, faces0, face_n, df, f_other, v0, handedness, angle, poly0),
        tri = (handedness >= 0) ? [s_f0, s_f1, s_o] : [s_f1, s_f0, s_o],
        lens = [ norm(tri[0]-tri[1]), norm(tri[1]-tri[2]), norm(tri[2]-tri[0]) ]
    )
    max(lens) - min(lens);

function _ps_snub_default_angle(poly, c_eff, handedness=1, steps=60, eps=1e-9) =
    let(
        angs = [for (i = [0:1:steps]) 0 + 60 * i / steps],
        cands = [
            for (a = angs)
                let(sp = _ps_snub_edge_spread(poly, c_eff, a, handedness))
                if (!is_undef(sp)) [a, sp]
        ],
        _ = assert(len(cands) > 0, "snub: no valid angle candidates"),
        errs = [for (c = cands) c[1]],
        e_min = min(errs),
        idx = [for (i = [0:1:len(errs)-1]) if (abs(errs[i] - e_min) <= eps) i][0]
    )
    cands[idx][0];

function _ps_snub_face_point(verts0, faces0, face_n, df, fi, v0, handedness, angle, poly0) =
    let(
        f = faces0[fi],
        n_f = face_n[fi],
        center = poly_face_center(poly0, fi, 1),
        ex = poly_face_ex(poly0, fi, 1),
        ey = poly_face_ey(poly0, fi, 1),
        p_base = center + df * n_f,
        v_rel = verts0[v0] - center,
        p0_2d = _ps_rot2d([v_dot(v_rel, ex), v_dot(v_rel, ey)], handedness * angle)
    )
    p_base + ex * p0_2d[0] + ey * p0_2d[1];

function poly_snub(poly, angle=undef, t=undef, c=undef, handedness=1, family_k=undef, eps=1e-8, len_eps=1e-6) =
    let(
        _ = assert(poly_valid(poly, "struct"), "snub: requires structurally valid poly"),
        c_eff = is_undef(c) ? 0.5 : c,
        fam = is_undef(family_k) ? ps_face_family_mode(poly)[0] : family_k,
        df_eff = is_undef(family_k)
            ? _ps_cantellate_df_from_c(poly, c_eff)
            : (c_eff * ps_face_radius_stat(poly, fam)),
        // Default to no twist unless an explicit angle is provided.
        angle_eff = is_undef(angle) ? 0 : angle
    )
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        face_pts = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    n_f = face_n[fi],
                    center = poly_face_center(poly0, fi, 1),
                    ex = poly_face_ex(poly0, fi, 1),
                    ey = poly_face_ey(poly0, fi, 1)
                )
                [
                    for (k = [0:1:n-1])
                        let(
                            v = f[k],
                            p_base = center + df_eff * n_f,
                            v_rel = verts0[v] - center,
                            p0_2d = _ps_rot2d([v_dot(v_rel, ex), v_dot(v_rel, ey)], handedness * angle_eff)
                        )
                        p_base + ex * p0_2d[0] + ey * p0_2d[1]
                ]
        ],
        face_offsets = _ps_face_offsets(faces0),
        sites = [
            for (fi = [0:1:len(faces0)-1])
                for (v = faces0[fi])
                    [fi, v]
        ],
        site_points = [
            for (fi = [0:1:len(faces0)-1])
                for (p = face_pts[fi])
                    p
        ],
        face_cycles = [
            for (fi = [0:1:len(faces0)-1])
                let(n = len(faces0[fi]))
                [ for (k = [0:1:n-1]) [1, face_offsets[fi] + k] ]
        ],
        edge_cycles = [
            for (ei = [0:1:len(edges)-1])
                let(
                    e = edges[ei],
                    fpair = edge_faces[ei],
                    f0 = fpair[0],
                    f1 = fpair[1],
                    v0 = e[0],
                    v1 = e[1],
                    pos0 = _ps_index_of(faces0[f0], v0),
                    pos1 = _ps_index_of(faces0[f0], v1),
                    pos2 = _ps_index_of(faces0[f1], v1),
                    pos3 = _ps_index_of(faces0[f1], v0)
                )
                [
                    [1, face_offsets[f0] + pos0],
                    [1, face_offsets[f0] + pos1],
                    [1, face_offsets[f1] + pos2],
                    [1, face_offsets[f1] + pos3]
                ]
        ],
        vert_cycles = [
            for (vi = [0:1:len(verts0)-1])
                let(fc = faces_around_vertex(poly0, vi, edges, edge_faces))
                [
                    for (fi = fc)
                        let(pos = _ps_index_of(faces0[fi], vi))
                            [1, face_offsets[fi] + pos]
                ]
        ],
        cycles_all = concat(face_cycles, edge_cycles, vert_cycles),
        q = ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps),
        verts = poly_verts(q),
        faces_q = poly_faces(q),
        face_count = len(faces0),
        edge_count = len(edges),
        face_faces = [for (i = [0:1:face_count-1]) faces_q[i]],
        edge_faces_q = [for (i = [0:1:edge_count-1]) faces_q[face_count + i]],
        vert_faces = [for (i = [0:1:len(faces_q)-face_count-edge_count-1]) faces_q[face_count + edge_count + i]],
        edge_tris = [
            for (f = edge_faces_q)
                let(a = f[0], b = f[1], c = f[2], d = f[3])
                (handedness >= 0)
                    ? [ [a, b, c], [a, c, d] ]
                    : [ [a, b, d], [b, c, d] ]
        ],
        edge_tris_flat = [for (t = edge_tris) each t],
        faces_new = concat(face_faces, edge_tris_flat, vert_faces),
        faces_oriented = ps_orient_all_faces_outward(verts, faces_new)
    )
    make_poly(verts, faces_oriented, poly_e_over_ir(q));

// Cantitruncation: truncation + cantellation (two parameters).
// t controls face-plane shift (like chamfer), c controls edge/vertex expansion (like cantellate).
function poly_cantitruncate(poly, t=undef, c=undef, eps = 1e-8, len_eps = 1e-6) =
    let(
        sol = (is_undef(t) && is_undef(c) && _ps_is_regular_base(poly)) ? solve_cantitruncate_trig(poly) : [t, c],
        t_eff = is_undef(sol[0]) ? _ps_truncate_default_t(poly) : sol[0],
        c_eff = is_undef(sol[1]) ? 0 : sol[1]
    )
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        // scale by inter-radius for consistent parameterization
        ir = min([for (e = edges) norm((verts0[e[0]] + verts0[e[1]]) / 2)]),
        d_f = -c_eff * ir,
        d_e = c_eff * ir,
        face_offsets = _ps_face_offsets(faces0),
        // edge points (truncation-style): two per edge
        edge_pts = _ps_edge_points(verts0, edges, t_eff),
        // face-vertex sites: one per (face, vertex)
        face_sites = [
            for (fi = [0:1:len(faces0)-1])
                for (v = faces0[fi])
                    [fi, v]
        ],
        face_pts3d = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    n_f = face_n[fi],
                    center = poly_face_center(poly, fi, 1),
                    ex = poly_face_ex(poly, fi, 1),
                    ey = poly_face_ey(poly, fi, 1),
                    p0 = center - n_f * d_f,
                    pts2d = [
                        for (k = [0:1:n-1])
                            let(p = verts0[f[k]] - center)
                                [v_dot(p, ex), v_dot(p, ey)]
                    ],
                    inset2d = _ps_face_inset_bisector_2d(f, fi, d_f, d_e, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0)
                )
                [
                    for (k = [0:1:n-1])
                        p0 + ex * inset2d[k][0] + ey * inset2d[k][1]
                ]
        ],
        // face-edge points: project edge points onto shifted face plane
        face_edge_pts3d = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    n_f = face_n[fi],
                    center = poly_face_center(poly, fi, 1),
                    p0 = center - n_f * d_f
                )
                [
                    for (k = [0:1:n-1])
                        let(
                            v0 = f[k],
                            v1 = f[(k+1)%n],
                            p_pair = _ps_project_edge_pts_for_face_edge(verts0, edges, edge_pts, n_f, p0, v0, v1),
                            p0o = p_pair[0],
                            p1o = p_pair[1]
                        )
                        each [p0o, p1o]
                ]
        ],
        face_edge_offsets = _ps_face_edge_offsets(faces0),
        face_pts_flat = [ for (fi = [0:1:len(faces0)-1]) for (p = face_pts3d[fi]) p ],
        face_edge_pts_flat = [ for (fi = [0:1:len(faces0)-1]) for (p = face_edge_pts3d[fi]) p ],
        edge_site_offset = 0,
        face_site_offset = len(face_edge_pts_flat),
        sites = concat(
            [ for (i = [0:1:len(face_edge_pts_flat)-1]) [0, i] ],
            face_sites
        ),
        site_points = concat(face_edge_pts_flat, face_pts_flat),
        face_cycles = _ps_face_cycles_from_face_edge_sites(faces0, [for (x = face_edge_offsets) edge_site_offset + x]),
        edge_cycles = _ps_edge_cycles_from_face_edge_sites(faces0, edges, edge_faces, [for (x = face_edge_offsets) edge_site_offset + x]),
        vert_cycles = _ps_vert_cycles_from_face_edge_sites(verts0, faces0, edges, edge_faces, [for (x = face_edge_offsets) edge_site_offset + x], poly0),
        // Debug: echo one decagon face cycle points and their angles (first face with n>=5*2)
        cycles_all = concat(face_cycles, edge_cycles, vert_cycles)
    )
    ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps);

// Cantitruncation with per-face-family c values (indexed by face size).
// c_by_size: list of [face_size, c] pairs; default_c used if size not found.
function poly_cantitruncate_families(poly, t, c_by_size, default_c=0, c_edge_by_pair=undef, eps=1e-8, len_eps=1e-6) =
    let(
        t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t,
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        ir = min([for (e = edges) norm((verts0[e[0]] + verts0[e[1]]) / 2)]),
        face_offsets = _ps_face_offsets(faces0),
        // edge points (truncation-style): two per edge
        edge_pts = _ps_edge_points(verts0, edges, t_eff),
        face_sites = [
            for (fi = [0:1:len(faces0)-1])
                for (v = faces0[fi])
                    [fi, v]
        ],
        face_pts3d = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    c_face = _ps_map_face_c(n, c_by_size, default_c),
                    d_f = -c_face * ir,
                    n_f = face_n[fi],
                    center = poly_face_center(poly, fi, 1),
                    ex = poly_face_ex(poly, fi, 1),
                    ey = poly_face_ey(poly, fi, 1),
                    p0 = center - n_f * d_f,
                    pts2d = [
                        for (k = [0:1:n-1])
                            let(p = verts0[f[k]] - center)
                                [v_dot(p, ex), v_dot(p, ey)]
                    ],
                    d_e_edges = [
                        for (k = [0:1:n-1])
                            let(
                                v0 = f[k],
                                v1 = f[(k+1)%n],
                                ei = ps_find_edge_index(edges, v0, v1),
                                adj = edge_faces[ei],
                                f_adj = (len(adj) < 2) ? fi : ((adj[0] == fi) ? adj[1] : adj[0]),
                                c_edge = is_undef(c_edge_by_pair)
                                    ? c_face
                                    : _ps_map_edge_c(n, len(faces0[f_adj]), c_edge_by_pair, c_face)
                            )
                            c_edge * ir
                    ],
                    inset2d = _ps_face_inset_bisector_2d(f, fi, d_f, d_e_edges, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0)
                )
                [
                    for (k = [0:1:n-1])
                        p0 + ex * inset2d[k][0] + ey * inset2d[k][1]
                ]
        ],
        face_edge_pts3d = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    c_face = _ps_map_face_c(n, c_by_size, default_c),
                    d_f = -c_face * ir,
                    n_f = face_n[fi],
                    center = poly_face_center(poly, fi, 1),
                    p0 = center - n_f * d_f
                )
                [
                    for (k = [0:1:n-1])
                        let(
                            v0 = f[k],
                            v1 = f[(k+1)%n],
                            p_pair = _ps_project_edge_pts_for_face_edge(verts0, edges, edge_pts, n_f, p0, v0, v1),
                            p0o = p_pair[0],
                            p1o = p_pair[1]
                        )
                        each [p0o, p1o]
                ]
        ],
        face_edge_offsets = _ps_face_edge_offsets(faces0),
        face_pts_flat = [ for (fi = [0:1:len(faces0)-1]) for (p = face_pts3d[fi]) p ],
        face_edge_pts_flat = [ for (fi = [0:1:len(faces0)-1]) for (p = face_edge_pts3d[fi]) p ],
        edge_site_offset = 0,
        face_site_offset = len(face_edge_pts_flat),
        sites = concat(
            [ for (i = [0:1:len(face_edge_pts_flat)-1]) [0, i] ],
            face_sites
        ),
        site_points = concat(face_edge_pts_flat, face_pts_flat),
        face_cycles = _ps_face_cycles_from_face_edge_sites(faces0, [for (x = face_edge_offsets) edge_site_offset + x]),
        edge_cycles = _ps_edge_cycles_from_face_edge_sites(faces0, edges, edge_faces, [for (x = face_edge_offsets) edge_site_offset + x]),
        vert_cycles = _ps_vert_cycles_from_face_edge_sites(verts0, faces0, edges, edge_faces, [for (x = face_edge_offsets) edge_site_offset + x], poly0),
        cycles_all = concat(face_cycles, edge_cycles, vert_cycles)
    )
    ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps);



// Compute per-corner truncation estimate t = 1/(2+r),
// where r = |B-C| / mean(|V-B|,|V-C|) for face corner (...B,V,C...).
function _ps_corner_t(verts, face, k) =
    let(
        n  = len(face),
        vm = face[(k-1+n) % n],
        v  = face[k],
        vp = face[(k+1) % n],
        V  = verts[v],
        B  = verts[vm],
        C  = verts[vp],

        LB = norm(V - B),
        LC = norm(V - C),
        L  = (LB + LC) / 2,

        chord = norm(B - C),
        r = (L == 0) ? 0 : chord / L,

        t = (2 + r == 0) ? 0 : 1 / (2 + r)
    )
    t;

// Return a “best guess” truncation t if user passes t=undef.
// - If the poly is locally uniform, returns the mean corner t.
// - Otherwise returns fallback (default 0.2).
//
// tol is an absolute tolerance on (t_max - t_min).
function _ps_truncate_default_t(poly, tol = 1e-3, fallback = 0.2) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),

        // gather per-corner t estimates
        ts = [
            for (f = faces)
                for (k = [0 : len(f)-1])
                    _ps_corner_t(verts, f, k)
        ],

        _ = assert(len(ts) > 0, "ps_truncate_default_t: no face corners found"),

        tmin = min(ts),
        tmax = max(ts),
        tavg = sum(ts) / len(ts),

        // decide if "uniform enough"
        ok = (tmax - tmin) <= tol,

        t0 = ok ? tavg : fallback
    )
    t0;
