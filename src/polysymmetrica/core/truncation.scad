// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <duals.scad>  // for faces_around_vertex helpers
use <classify.scad>
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
// d_f shifts the face plane outward (along face normal) in caller space; note the
// internal p0 uses center - n_f*d_f so positive d_f still means outward overall.
// d_e shifts the edge-bisector planes along their normals.
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

function _ps_cantellate_df_map(poly, df_max=undef, steps=16, family_edge_idx=0) =
    let(
        verts = poly_verts(poly),
        edges = _ps_edges_from_faces(poly_faces(poly)),
        e0 = edges[0],
        edge_len = norm(verts[e0[1]] - verts[e0[0]]),
        ir = edge_len / poly_e_over_ir(poly),
        df_max_eff = is_undef(df_max) ? (2 * ir) : df_max,
        df_mid = cantellate_square_df(poly, 0, df_max_eff, steps, family_edge_idx)
    )
    [df_mid, df_max_eff, ir];

function _ps_cantellate_df_from_c_linear(c, df_mid, df_max_eff) =
    (c <= 0.5)
        ? (2 * c * df_mid)
        : (df_mid + (c - 0.5) * 2 * (df_max_eff - df_mid));

function _ps_cantellate_c_from_df_linear(df, df_mid, df_max_eff, eps=1e-12) =
    (df <= df_mid)
        ? (((2 * df_mid) <= eps) ? 0 : (df / (2 * df_mid)))
        : (((2 * (df_max_eff - df_mid)) <= eps) ? 0.5 : (0.5 + (df - df_mid) / (2 * (df_max_eff - df_mid))));

function _ps_cantellate_df_from_c(poly, c, df_max=undef, steps=16, family_edge_idx=0) =
    let(
        map = _ps_cantellate_df_map(poly, df_max, steps, family_edge_idx),
        df_mid = map[0],
        df_max_eff = map[1]
    )
    _ps_cantellate_df_from_c_linear(c, df_mid, df_max_eff);

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

function _ps_snub_face_cached_point(face_pts, faces0, fi, v) =
    let(pos = _ps_index_of(faces0[fi], v))
    face_pts[fi][pos];

function _ps_snub_required_faces(edge_faces, edge_list) =
    _ps_unique_ints([for (ei = edge_list) each edge_faces[ei]]);

// Uniformity objective for snub defaults:
// collect all edge types induced by edge-face quads and minimize global spread.
function _ps_snub_uniform_error_base(verts0, faces0, edges, edge_faces, face_n, poly0, d_f, d_e, angle, handedness=1, edge_reps=undef, eps=1e-12) =
    let(ev = _ps_snub_eval_errors_base(verts0, faces0, edges, edge_faces, face_n, poly0, d_f, d_e, angle, handedness, edge_reps, eps))
    ev[0];

// Combined snub objective terms in one pass:
// returns [uniform_spread, edge_triangle_error].
function _ps_snub_eval_errors_base(verts0, faces0, edges, edge_faces, face_n, poly0, d_f, d_e, angle, handedness=1, edge_reps=undef, eps=1e-12) =
    let(
        d_f_by_face = [for (_ = [0:1:len(faces0)-1]) d_f],
        edge_list = is_undef(edge_reps) ? [for (ei = [0:1:len(edges)-1]) ei] : edge_reps,
        req_faces = _ps_snub_required_faces(edge_faces, edge_list),
        face_pts = [
            for (fi = [0:1:len(faces0)-1])
                if (_ps_list_contains(req_faces, fi))
                    let(
                        f = faces0[fi],
                        center = poly_face_center(poly0, fi, 1),
                        ex = poly_face_ex(poly0, fi, 1),
                        ey = poly_face_ey(poly0, fi, 1),
                        n_f = face_n[fi],
                        d_fi = d_f_by_face[fi],
                        p_base = center + d_fi * n_f,
                        pts2d = [
                            for (k = [0:1:len(f)-1])
                                let(p = verts0[f[k]] - center)
                                    [v_dot(p, ex), v_dot(p, ey)]
                        ],
                        inset2d = _ps_face_inset_bisector_2d(f, fi, -d_fi, -d_e, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0),
                        pts2d_rot = [for (p2 = inset2d) _ps_rot2d(p2, handedness * angle)]
                    )
                    [for (p2 = pts2d_rot) p_base + ex * p2[0] + ey * p2[1]]
                else
                    undef
        ],
        per_edge = [
            for (ei = edge_list)
                let(
                    e = edges[ei],
                    fpair = edge_faces[ei],
                    f0 = fpair[0],
                    f1 = fpair[1],
                    v0 = e[0],
                    v1 = e[1],
                    p00 = _ps_snub_face_cached_point(face_pts, faces0, f0, v0),
                    p01 = _ps_snub_face_cached_point(face_pts, faces0, f0, v1),
                    p11 = _ps_snub_face_cached_point(face_pts, faces0, f1, v1),
                    p10 = _ps_snub_face_cached_point(face_pts, faces0, f1, v0),
                    l_face0 = norm(p00 - p01),
                    l_face1 = norm(p10 - p11),
                    l_vert0 = norm(p00 - p10),
                    l_vert1 = norm(p01 - p11),
                    l_diag = min(norm(p00 - p11), norm(p01 - p10)),
                    tris_a = [ [p00, p01, p11], [p00, p11, p10] ],
                    tris_b = [ [p00, p01, p10], [p01, p11, p10] ],
                    err_a = sum([
                        for (tri = tris_a)
                            let(
                                lens = [ norm(tri[0]-tri[1]), norm(tri[1]-tri[2]), norm(tri[2]-tri[0]) ],
                                avg = (lens[0] + lens[1] + lens[2]) / 3,
                                r0 = lens[0] / avg,
                                r1 = lens[1] / avg,
                                r2 = lens[2] / avg
                            )
                                (pow(r0 - 1, 2) + pow(r1 - 1, 2) + pow(r2 - 1, 2))
                    ]) / 2,
                    err_b = sum([
                        for (tri = tris_b)
                            let(
                                lens = [ norm(tri[0]-tri[1]), norm(tri[1]-tri[2]), norm(tri[2]-tri[0]) ],
                                avg = (lens[0] + lens[1] + lens[2]) / 3,
                                r0 = lens[0] / avg,
                                r1 = lens[1] / avg,
                                r2 = lens[2] / avg
                            )
                                (pow(r0 - 1, 2) + pow(r1 - 1, 2) + pow(r2 - 1, 2))
                    ]) / 2
                )
                [ [l_face0, l_face1, l_vert0, l_vert1, l_diag], min(err_a, err_b) ]
        ],
        lens_all = [for (x = per_edge) each x[0]],
        tri_errs = [for (x = per_edge) x[1]],
        good = (len(lens_all) > 0) && (len([for (l = lens_all) if (is_undef(l) || l <= eps) 1]) == 0),
        avg = good ? (sum(lens_all) / len(lens_all)) : undef,
        err_u = (!good || avg <= eps) ? undef : ((max(lens_all) - min(lens_all)) / avg),
        err_t = (len(tri_errs) == 0) ? undef : (sum(tri_errs) / len(tri_errs))
    )
    [err_u, err_t];

function _ps_snub_obj_from_errors(ev, bad=1e9) =
    (is_undef(ev) || len(ev) < 2 || is_undef(ev[0]) || is_undef(ev[1])) ? bad : max([ev[0], sqrt(ev[1])]);

function _ps_clamp(x, lo, hi) = min(max(x, lo), hi);

function _ps_snub_obj_regular(verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff, c, r, a, handedness=1, edge_reps=undef) =
    let(
        c0 = _ps_clamp(c, 0, 1),
        de = _ps_cantellate_df_from_c_linear(c0, df_mid, df_max_eff),
        rr = max(0.01, r),
        df = rr * de,
        ev = _ps_snub_eval_errors_base(verts0, faces0, edges, edge_faces, face_n, poly0, df, de, a, handedness, edge_reps)
    )
    _ps_snub_obj_from_errors(ev);

function _ps_snub_best3_regular(verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff, c, r, a, dc, dr, da, c_max, a_max, handedness=1, edge_reps=undef, bad=1e9) =
    let(
        cands = [
            for (ic = [-1:1:1])
                for (ir = [-1:1:1])
                    for (ia = [-1:1:1])
                        let(
                            c1 = _ps_clamp(c + ic * dc, 0, c_max),
                            r1 = _ps_clamp(r + ir * dr, 0.5, 1.5),
                            a1 = _ps_clamp(a + ia * da, 0, a_max),
                            e1 = _ps_snub_obj_regular(verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff, c1, r1, a1, handedness, edge_reps)
                        )
                        [c1, r1, a1, e1]
        ],
        vals = [for (x = cands) x[3]],
        e_min = min(vals),
        idx = [for (i = [0:1:len(vals)-1]) if (abs(vals[i] - e_min) <= 1e-12) i][0]
    )
    cands[idx];

function _ps_snub_local_search_regular(verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff, state, steps, c_max, a_max, handedness=1, edge_reps=undef, iter=0, max_iter=10, eps=1e-9) =
    let(
        c = state[0], r = state[1], a = state[2], e = state[3],
        dc = steps[0], dr = steps[1], da = steps[2],
        done = (iter >= max_iter) || (max([dc, dr, da]) <= 1e-4),
        best = _ps_snub_best3_regular(verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff, c, r, a, dc, dr, da, c_max, a_max, handedness, edge_reps),
        improved = best[3] + eps < e,
        next_state = improved ? best : state,
        next_steps = improved ? steps : [dc * 0.6, dr * 0.6, da * 0.6]
    )
    done ? state : _ps_snub_local_search_regular(verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff, next_state, next_steps, c_max, a_max, handedness, edge_reps, iter + 1, max_iter, eps);

function _ps_family_id_for_index(idx, fams) =
    let(hits = [for (fi = [0:1:len(fams)-1]) if (_ps_list_contains(fams[fi][1], idx)) fi])
    (len(hits) == 0) ? -1 : hits[0];

function _ps_family_ids_from_fams(n, fams) =
    [for (i = [0:1:n-1]) _ps_family_id_for_index(i, fams)];

function _ps_replace_at(list, idx, val) =
    [for (i = [0:1:len(list)-1]) (i == idx) ? val : list[i]];

function _ps_snub_eval_errors_face_df_base(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, d_f_by_family, d_e, angle, handedness=1, edge_reps=undef, eps=1e-12) =
    let(
        d_f_by_face = [
            for (fi = [0:1:len(faces0)-1])
                let(id = face_fid[fi])
                ((id >= 0 && id < len(d_f_by_family)) ? d_f_by_family[id] : d_e)
        ],
        edge_list = is_undef(edge_reps) ? [for (ei = [0:1:len(edges)-1]) ei] : edge_reps,
        req_faces = _ps_snub_required_faces(edge_faces, edge_list),
        face_pts = [
            for (fi = [0:1:len(faces0)-1])
                if (_ps_list_contains(req_faces, fi))
                    let(
                        f = faces0[fi],
                        center = poly_face_center(poly0, fi, 1),
                        ex = poly_face_ex(poly0, fi, 1),
                        ey = poly_face_ey(poly0, fi, 1),
                        n_f = face_n[fi],
                        d_fi = d_f_by_face[fi],
                        p_base = center + d_fi * n_f,
                        pts2d = [
                            for (k = [0:1:len(f)-1])
                                let(p = verts0[f[k]] - center)
                                    [v_dot(p, ex), v_dot(p, ey)]
                        ],
                        inset2d = _ps_face_inset_bisector_2d(f, fi, -d_fi, -d_e, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0),
                        pts2d_rot = [for (p2 = inset2d) _ps_rot2d(p2, handedness * angle)]
                    )
                    [for (p2 = pts2d_rot) p_base + ex * p2[0] + ey * p2[1]]
                else
                    undef
        ],
        per_edge = [
            for (ei = edge_list)
                let(
                    e = edges[ei],
                    fpair = edge_faces[ei],
                    f0 = fpair[0],
                    f1 = fpair[1],
                    v0 = e[0],
                    v1 = e[1],
                    p00 = _ps_snub_face_cached_point(face_pts, faces0, f0, v0),
                    p01 = _ps_snub_face_cached_point(face_pts, faces0, f0, v1),
                    p11 = _ps_snub_face_cached_point(face_pts, faces0, f1, v1),
                    p10 = _ps_snub_face_cached_point(face_pts, faces0, f1, v0),
                    l_face0 = norm(p00 - p01),
                    l_face1 = norm(p10 - p11),
                    l_vert0 = norm(p00 - p10),
                    l_vert1 = norm(p01 - p11),
                    l_diag = min(norm(p00 - p11), norm(p01 - p10)),
                    tris_a = [ [p00, p01, p11], [p00, p11, p10] ],
                    tris_b = [ [p00, p01, p10], [p01, p11, p10] ],
                    err_a = sum([
                        for (tri = tris_a)
                            let(
                                lens = [ norm(tri[0]-tri[1]), norm(tri[1]-tri[2]), norm(tri[2]-tri[0]) ],
                                avg = (lens[0] + lens[1] + lens[2]) / 3,
                                r0 = lens[0] / avg,
                                r1 = lens[1] / avg,
                                r2 = lens[2] / avg
                            )
                                (pow(r0 - 1, 2) + pow(r1 - 1, 2) + pow(r2 - 1, 2))
                    ]) / 2,
                    err_b = sum([
                        for (tri = tris_b)
                            let(
                                lens = [ norm(tri[0]-tri[1]), norm(tri[1]-tri[2]), norm(tri[2]-tri[0]) ],
                                avg = (lens[0] + lens[1] + lens[2]) / 3,
                                r0 = lens[0] / avg,
                                r1 = lens[1] / avg,
                                r2 = lens[2] / avg
                            )
                                (pow(r0 - 1, 2) + pow(r1 - 1, 2) + pow(r2 - 1, 2))
                    ]) / 2
                )
                [ [l_face0, l_face1, l_vert0, l_vert1, l_diag], min(err_a, err_b) ]
        ],
        lens_all = [for (x = per_edge) each x[0]],
        tri_errs = [for (x = per_edge) x[1]],
        good = (len(lens_all) > 0) && (len([for (l = lens_all) if (is_undef(l) || l <= eps) 1]) == 0),
        avg = good ? (sum(lens_all) / len(lens_all)) : undef,
        err_u = (!good || avg <= eps) ? undef : ((max(lens_all) - min(lens_all)) / avg),
        err_t = (len(tri_errs) == 0) ? undef : (sum(tri_errs) / len(tri_errs))
    )
    [err_u, err_t];

function _ps_snub_obj_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_by_family, df_mid, df_max_eff, c, a, handedness=1, edge_reps=undef) =
    let(
        c0 = _ps_clamp(c, 0, 1),
        de = _ps_cantellate_df_from_c_linear(c0, df_mid, df_max_eff),
        d_f_by_family = [for (k = [0:1:len(r_by_family)-1]) r_by_family[k] * de],
        ev = _ps_snub_eval_errors_face_df_base(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, d_f_by_family, de, a, handedness, edge_reps)
    )
    _ps_snub_obj_from_errors(ev);

function _ps_snub_best_ca_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_by_family, df_mid, df_max_eff, c, a, dc, da, c_max, a_max, handedness=1, edge_reps=undef) =
    let(
        cands = [
            for (ic = [-1:1:1])
                for (ia = [-1:1:1])
                    let(
                        c1 = _ps_clamp(c + ic * dc, 0, c_max),
                        a1 = _ps_clamp(a + ia * da, 0, a_max),
                        e1 = _ps_snub_obj_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_by_family, df_mid, df_max_eff, c1, a1, handedness, edge_reps)
                    )
                    [c1, a1, e1]
        ],
        vals = [for (x = cands) x[2]],
        e_min = min(vals),
        idx = [for (i = [0:1:len(vals)-1]) if (abs(vals[i] - e_min) <= 1e-12) i][0]
    )
    cands[idx];

function _ps_snub_local_search_ca_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_by_family, df_mid, df_max_eff, state, steps, c_max, a_max, handedness=1, edge_reps=undef, iter=0, max_iter=12, eps=1e-9) =
    let(
        c = state[0], a = state[1], e = state[2],
        dc = steps[0], da = steps[1],
        done = (iter >= max_iter) || (max([dc, da]) <= 1e-4),
        best = _ps_snub_best_ca_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_by_family, df_mid, df_max_eff, c, a, dc, da, c_max, a_max, handedness, edge_reps),
        improved = best[2] + eps < e,
        next_state = improved ? best : state,
        next_steps = improved ? steps : [dc * 0.6, da * 0.6]
    )
    done ? state : _ps_snub_local_search_ca_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_by_family, df_mid, df_max_eff, next_state, next_steps, c_max, a_max, handedness, edge_reps, iter + 1, max_iter, eps);

function _ps_snub_default_angle_df(poly, df, handedness=1, steps=60, a_max=35, eps=1e-9) =
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        is_reg = _ps_is_regular_base(poly),
        cls = is_reg ? undef : poly_classify(poly, 1),
        edge_reps = is_reg ? [0] : [for (f = cls[1]) f[1][0]],
        angs = [for (i = [0:1:steps]) a_max * i / steps],
        cands = [
            for (a = angs)
                let(err = _ps_snub_uniform_error_base(verts0, faces0, edges, edge_faces, face_n, poly0, df, df, a, handedness, edge_reps))
                if (!is_undef(err)) [a, err]
        ],
        _ = assert(len(cands) > 0, "snub: no valid angle candidates"),
        errs = [for (c = cands) c[1]],
        e_min = min(errs),
        idx = [for (i = [0:1:len(errs)-1]) if (abs(errs[i] - e_min) <= eps) i][0],
        idx20 = [for (i = [0:1:len(cands)-1]) if (abs(cands[i][0] - 20) <= 1e-6) i],
        a0 = cands[0][1],
        a20 = (len(idx20) > 0) ? cands[idx20[0]][1] : undef,
        _0 = echo(str(
            "snub: angle solve (fixed df) min_err=", e_min,
            " at angle=", cands[idx][0],
            " (df=", df, ")",
            " err@angle0=", a0,
            is_undef(a20) ? "" : str(" err@angle20=", a20)
        ))
    )
    cands[idx][0];

// Solve angle for fixed c using representative-edge objective.
function _ps_snub_default_angle_c(poly, c, df=undef, handedness=1, steps=16, a_max=30, eps=1e-9) =
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        is_reg = _ps_is_regular_base(poly),
        cls = is_reg ? undef : poly_classify(poly, 1),
        edge_reps = is_reg ? [0] : [for (f = cls[1]) f[1][0]],
        map = _ps_cantellate_df_map(poly, steps=6),
        df_mid = map[0],
        df_max_eff = map[1],
        de = _ps_cantellate_df_from_c_linear(c, df_mid, df_max_eff),
        d_f = is_undef(df) ? de : df,
        angs = [for (i = [0:1:steps]) a_max * i / steps],
        cands = [
            for (a = angs)
                let(
                    ev = _ps_snub_eval_errors_base(verts0, faces0, edges, edge_faces, face_n, poly0, d_f, de, a, handedness, edge_reps),
                    err = _ps_snub_obj_from_errors(ev)
                )
                [a, err]
        ],
        _ = assert(len(cands) > 0, "snub: no valid angle candidates for fixed c"),
        errs = [for (c0 = cands) c0[1]],
        e_min = min(errs),
        idx = [for (i = [0:1:len(errs)-1]) if (abs(errs[i] - e_min) <= eps) i][0],
        a0 = cands[0][1],
        idx20 = [for (i = [0:1:len(cands)-1]) if (abs(cands[i][0] - 20) <= 1e-6) i],
        a20 = (len(idx20) > 0) ? cands[idx20[0]][1] : undef,
        _0 = echo(str(
            "snub: angle solve (fixed c) min_err=", e_min,
            " at angle=", cands[idx][0],
            " (c=", c, ")",
            " err@angle0=", a0,
            is_undef(a20) ? "" : str(" err@angle20=", a20)
        ))
    )
    cands[idx][0];

// Solve for default snub parameters with a tiered strategy:
// regular -> family representative -> bounded heuristic.
// Returns [df, angle, c, tier].
function _ps_snub_default_params(poly, handedness=1, eps=1e-9) =
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        is_reg = _ps_is_regular_base(poly),
        cls = is_reg ? undef : poly_classify(poly, 1),
        cmap = _ps_cantellate_df_map(poly, steps=6),
        df_mid = cmap[0],
        df_max_eff = cmap[1],
        ff = is_reg ? 1 : len(cls[0]),
        ef = is_reg ? 1 : len(cls[1]),
        vf = is_reg ? 1 : len(cls[2]),
        tier = is_reg ? "regular" : ((ff <= 3 && ef <= 4 && vf <= 4) ? "family" : "heuristic"),
        c_steps = (tier == "family") ? 16 : 10,
        a_steps = (tier == "family") ? 18 : 14,
        c_max = (tier == "family") ? 0.2 : 0.25,
        a_max = 35,
        reg_c_steps = (len(edges) <= 12) ? 2 : 2,
        reg_df_steps = (len(edges) <= 12) ? 2 : 2,
        reg_a_steps = (len(edges) <= 12) ? 3 : 3,
        edge_reps_all = is_reg ? [0] : [for (f = cls[1]) f[1][0]],
        edge_reps = (tier == "heuristic" && len(edge_reps_all) > 12)
            ? [for (i = [0:1:11]) edge_reps_all[i]]
            : edge_reps_all,
        reg_best = is_reg ? _ps_snub_default_params_full(poly, handedness, c_steps=reg_c_steps, df_steps=reg_df_steps, a_steps=reg_a_steps, c_max=0.15, a_max=25, eps=eps, base=base, edge_reps=edge_reps_all) : undef,
        fam_best = (!is_reg && tier == "family") ? _ps_snub_default_params_family(poly, handedness, c_steps=12, a_steps=14, c_max=0.2, a_max=30, eps=eps, base=base, cls=cls) : undef,
        c_vals = (is_reg || tier == "family") ? [] : [for (i = [1:1:c_steps]) c_max * i / (c_steps + 1)],
        angs = (is_reg || tier == "family") ? [] : [for (i = [0:1:a_steps]) a_max * i / a_steps],
        cands = (is_reg || tier == "family") ? [] : [
            for (c = c_vals)
                let(
                    df = _ps_cantellate_df_from_c_linear(c, df_mid, df_max_eff),
                    errs = [
                        for (a = angs)
                            let(err = _ps_snub_uniform_error_base(verts0, faces0, edges, edge_faces, face_n, poly0, df, df, a, handedness, edge_reps))
                            if (!is_undef(err)) [a, err]
                    ],
                    ok = len(errs) > 0,
                    err_vals = ok ? [for (x = errs) x[1]] : [],
                    e_min = ok ? min(err_vals) : undef,
                    idx = ok ? [for (i = [0:1:len(err_vals)-1]) if (abs(err_vals[i] - e_min) <= eps) i][0] : undef,
                    a_best = ok ? errs[idx][0] : undef
                )
                if (ok) [df, a_best, c, e_min]
        ],
        _ = (is_reg || tier == "family") ? 0 : assert(len(cands) > 0, "snub: no valid (c, angle) candidates"),
        all_errs = (is_reg || tier == "family") ? [] : [for (c = cands) c[3]],
        e_min = is_reg ? reg_best[3] : ((tier == "family") ? fam_best[3] : min(all_errs)),
        idx = (is_reg || tier == "family") ? undef : [for (i = [0:1:len(all_errs)-1]) if (abs(all_errs[i] - e_min) <= eps) i][0],
        df_best = is_reg ? reg_best[0] : ((tier == "family") ? fam_best[0] : cands[idx][0]),
        a_best = is_reg ? reg_best[1] : ((tier == "family") ? fam_best[1] : cands[idx][1]),
        c_best = is_reg ? reg_best[2] : ((tier == "family") ? fam_best[2] : cands[idx][2]),
        d_f_by_family = (tier == "family") ? fam_best[4] : undef,
        de_best = _ps_cantellate_df_from_c_linear(c_best, df_mid, df_max_eff),
        face_fid = (tier == "family") ? _ps_family_ids_from_fams(len(faces0), cls[0]) : undef,
        err_a0 = (tier == "family")
            ? let(ev0 = _ps_snub_eval_errors_face_df_base(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, d_f_by_family, de_best, 0, handedness, edge_reps))
              ((is_undef(ev0[0]) || is_undef(ev0[1])) ? undef : max([ev0[0], sqrt(ev0[1])]))
            : _ps_snub_uniform_error_base(verts0, faces0, edges, edge_faces, face_n, poly0, df_best, de_best, 0, handedness, edge_reps),
        err_a20 = (tier == "family")
            ? let(ev20 = _ps_snub_eval_errors_face_df_base(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, d_f_by_family, de_best, 20, handedness, edge_reps))
              ((is_undef(ev20[0]) || is_undef(ev20[1])) ? undef : max([ev20[0], sqrt(ev20[1])]))
            : _ps_snub_uniform_error_base(verts0, faces0, edges, edge_faces, face_n, poly0, df_best, de_best, 20, handedness, edge_reps),
        _0 = echo(str(
            "snub: default auto tier=", tier,
            " families(f/e/v)=", ff, "/", ef, "/", vf,
            " min_err=", e_min,
            " c=", c_best,
            " df=", df_best,
            (tier == "family") ? str(" df_by_family=", d_f_by_family) : "",
            " angle=", a_best,
            " err@angle0=", err_a0,
            is_undef(err_a20) ? "" : str(" err@angle20=", err_a20)
        ))
    )
    is_undef(d_f_by_family) ? [df_best, a_best, c_best, tier] : [df_best, a_best, c_best, tier, d_f_by_family];

// Solve for a c/angle pair by minimizing representative-edge snub uniformity.
// This avoids full-poly rebuilds in the default path.
function _ps_snub_default_params_full(poly, handedness=1, c_steps=10, df_steps=8, a_steps=12, c_max=0.2, a_max=30, eps=1e-9, base=undef, edge_reps=undef) =
    let(
        base0 = is_undef(base) ? _ps_poly_base(poly) : base,
        verts0 = base0[0],
        faces0 = base0[1],
        edges = base0[2],
        edge_faces = base0[3],
        face_n = base0[4],
        poly0 = base0[5],
        cls = is_undef(edge_reps) ? poly_classify(poly, 1) : undef,
        edge_reps_eff = is_undef(edge_reps) ? [for (f = cls[1]) f[1][0]] : edge_reps,
        map = _ps_cantellate_df_map(poly, steps=6),
        df_mid = map[0],
        df_max_eff = map[1],
        // Coarse seed
        cs = [for (i = [1:1:c_steps]) c_max * i / (c_steps + 1)],
        rs = [for (i = [0:1:df_steps]) 0.7 + 0.6 * i / df_steps],
        angs = [for (i = [0:1:a_steps]) a_max * i / a_steps],
        seed_samples = len(cs) * len(rs) * len(angs),
        local_samples_bound = 27 * 24,
        seeds = [
            for (c = cs)
                for (r = rs)
                    for (a = angs)
                        let(e = _ps_snub_obj_regular(verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff, c, r, a, handedness, edge_reps_eff))
                        [c, r, a, e]
        ],
        _ = assert(len(seeds) > 0, "snub: no valid (c,r,angle) seeds"),
        seed_vals = [for (s = seeds) s[3]],
        seed_idx = [for (i = [0:1:len(seed_vals)-1]) if (abs(seed_vals[i] - min(seed_vals)) <= eps) i][0],
        seed = seeds[seed_idx],
        // Adaptive local solve
        dc0 = c_max / (c_steps + 1),
        dr0 = 0.6 / max(1, df_steps),
        da0 = a_max / max(1, a_steps),
        best = _ps_snub_local_search_regular(
            verts0, faces0, edges, edge_faces, face_n, poly0, df_mid, df_max_eff,
            seed, [dc0, dr0, da0], c_max, a_max, handedness, edge_reps_eff, 0, 24, eps
        ),
        c_best = best[0],
        r_best = best[1],
        a_best = best[2],
        e_min2 = best[3],
        de_best = _ps_cantellate_df_from_c_linear(c_best, df_mid, df_max_eff),
        df_best = r_best * de_best,
        _0 = echo(str(
            "snub: default param FULL search min_err=", e_min2,
            " c=", c_best,
            " df=", df_best,
            " angle=", a_best,
            " r=", r_best,
            " reps=", len(edge_reps_eff),
            " samples(seed/local<=)=", seed_samples, "/", local_samples_bound
        ))
    )
    [df_best, a_best, c_best, e_min2];

// Family-tier default solve: global c/angle with per-face-family df scaling.
// Returns [df_avg, angle, c, err, df_by_family].
function _ps_snub_default_params_family(poly, handedness=1, c_steps=12, a_steps=14, c_max=0.2, a_max=30, eps=1e-9, base=undef, cls=undef) =
    let(
        base0 = is_undef(base) ? _ps_poly_base(poly) : base,
        verts0 = base0[0],
        faces0 = base0[1],
        edges = base0[2],
        edge_faces = base0[3],
        face_n = base0[4],
        poly0 = base0[5],
        cls0 = is_undef(cls) ? poly_classify(poly, 1) : cls,
        face_fams = cls0[0],
        edge_reps = [for (f = cls0[1]) f[1][0]],
        ff = len(face_fams),
        face_fid = _ps_family_ids_from_fams(len(faces0), face_fams),
        r0 = [for (k = [0:1:ff-1]) 1],
        map = _ps_cantellate_df_map(poly, steps=6),
        df_mid = map[0],
        df_max_eff = map[1],
        cs = [for (i = [1:1:c_steps]) c_max * i / (c_steps + 1)],
        angs = [for (i = [0:1:a_steps]) a_max * i / a_steps],
        seed_samples = len(cs) * len(angs),
        cands = [
            for (c = cs)
                for (a = angs)
                    let(err = _ps_snub_obj_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r0, df_mid, df_max_eff, c, a, handedness, edge_reps))
                    [c, a, err]
        ],
        _ = assert(len(cands) > 0, "snub: family tier no valid (c, angle) candidates"),
        vals0 = [for (x = cands) x[2]],
        idx0 = [for (i = [0:1:len(vals0)-1]) if (abs(vals0[i] - min(vals0)) <= eps) i][0],
        c0 = cands[idx0][0],
        a0 = cands[idx0][1],
        de0 = _ps_cantellate_df_from_c_linear(c0, df_mid, df_max_eff),
        r_vals = [0.8, 0.9, 1.0, 1.1, 1.2],
        family_r_samples = ff * len(r_vals),
        r_best = [
            for (k = [0:1:ff-1])
                let(
                    candsk = [
                        for (r = r_vals)
                            let(
                                r_try = _ps_replace_at(r0, k, r),
                                d_f_try = [for (j = [0:1:ff-1]) r_try[j] * de0],
                                ev = _ps_snub_eval_errors_face_df_base(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, d_f_try, de0, a0, handedness, edge_reps),
                                e0 = ev[0],
                                e1 = ev[1],
                                err = (is_undef(e0) || is_undef(e1)) ? undef : max([e0, sqrt(e1)])
                            )
                            if (!is_undef(err)) [r, err]
                    ],
                    valsk = [for (x = candsk) x[1]],
                    idxk = [for (i = [0:1:len(valsk)-1]) if (abs(valsk[i] - min(valsk)) <= eps) i][0]
                )
                candsk[idxk][0]
        ],
        e0 = _ps_snub_obj_family(verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_best, df_mid, df_max_eff, c0, a0, handedness, edge_reps),
        dc0 = c_max / (c_steps + 1),
        da0 = a_max / max(1, a_steps),
        local_samples_bound = 9 * 12,
        best_ca = _ps_snub_local_search_ca_family(
            verts0, faces0, edges, edge_faces, face_n, poly0, face_fid, r_best, df_mid, df_max_eff,
            [c0, a0, e0], [dc0, da0], c_max, a_max, handedness, edge_reps, 0, 18, eps
        ),
        c_best = best_ca[0],
        a_best = best_ca[1],
        e_min = best_ca[2],
        de_best = _ps_cantellate_df_from_c_linear(c_best, df_mid, df_max_eff),
        d_f_by_family_best = [for (k = [0:1:ff-1]) r_best[k] * de_best],
        df_avg = sum(d_f_by_family_best) / max(1, len(d_f_by_family_best)),
        _0 = echo(str(
            "snub: family-tier solve min_err=", e_min,
            " c=", c_best,
            " angle=", a_best,
            " df_by_family=", d_f_by_family_best,
            " reps=", len(edge_reps),
            " samples(seed/r/local<=)=", seed_samples, "/", family_r_samples, "/", local_samples_bound
        ))
    )
    [df_avg, a_best, c_best, e_min, d_f_by_family_best];

function _ps_face_max_plane_err(verts, f) =
    let(
        n = ps_face_normal(verts, f),
        v0 = verts[f[0]],
        errs = [for (vid = f) abs(v_dot(n, verts[vid] - v0))]
    )
    (len(errs) == 0) ? 0 : max(errs);

function poly_snub(poly, angle=undef, c=undef, df=undef, de=undef, handedness=1, family_k=undef, eps=1e-8, len_eps=1e-6, df_by_family=undef) =
    let(
        _ = assert(poly_valid(poly, "struct"), "snub: requires structurally valid poly"),
        params = (is_undef(c) && is_undef(df) && is_undef(angle) && is_undef(family_k))
            ? _ps_snub_default_params(poly, handedness)
            : undef,
        _choice = !is_undef(params)
            ? echo(str("snub: using auto defaults tier=", params[3], " c=", params[2], " df=", params[0], " angle=", params[1]))
            : 0,
        fam = is_undef(family_k) ? ps_face_family_mode(poly)[0] : family_k,
        df_by_family_auto = (!is_undef(params) && len(params) > 4) ? params[4] : undef,
        df_by_family_eff = is_undef(df_by_family) ? df_by_family_auto : df_by_family,
        df_base = !is_undef(df)
            ? df
            : is_undef(c)
            ? (is_undef(params)
                ? (is_undef(family_k) ? _ps_cantellate_df_from_c(poly, 0.5) : (0.5 * ps_face_radius_stat(poly, fam)))
                : params[0])
            : _ps_cantellate_df_from_c(poly, c),
        // c controls edge/vertex inset (de). If unspecified, keep de coupled to df.
        de_base = !is_undef(de)
            ? de
            : !is_undef(c)
            ? _ps_cantellate_df_from_c(poly, c)
            : (!is_undef(params) ? _ps_cantellate_df_from_c(poly, params[2]) : df_base),
        // Default to equilateral twist unless an explicit angle is provided.
        angle_eff = is_undef(angle)
            ? let(
                a = !is_undef(params) ? params[1]
                    : (!is_undef(c) ? _ps_snub_default_angle_c(poly, c, df_base, handedness) : _ps_snub_default_angle_df(poly, df_base, handedness)),
                _ = echo(str("snub: angle unspecified, default=", a, " (df=", df_base, is_undef(c) ? ")" : str(", c=", c, ")")))
              ) a
            : angle
    )
    let(
        base = _ps_poly_base(poly),
        verts0 = base[0],
        faces0 = base[1],
        edges = base[2],
        edge_faces = base[3],
        face_n = base[4],
        poly0 = base[5],
        face_fid = is_undef(df_by_family_eff) ? undef : _ps_family_ids_from_fams(len(faces0), poly_classify(poly, 1)[0]),
        face_pts = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    n_f = face_n[fi],
                    center = poly_face_center(poly0, fi, 1),
                    ex = poly_face_ex(poly0, fi, 1),
                    ey = poly_face_ey(poly0, fi, 1),
                    f_id = is_undef(face_fid) ? -1 : face_fid[fi],
                    d_f = (is_undef(df_by_family_eff) || f_id < 0 || f_id >= len(df_by_family_eff)) ? df_base : df_by_family_eff[f_id],
                    d_e = de_base,
                    pts2d = [
                        for (k = [0:1:n-1])
                            let(p = verts0[f[k]] - center)
                                [v_dot(p, ex), v_dot(p, ey)]
                    ],
                    // d_f is outward for callers; _ps_face_inset_bisector_2d uses p0 = center - n_f*d_f.
                    // d_e follows the same outward convention.
                    inset2d = _ps_face_inset_bisector_2d(f, fi, -d_f, -d_e, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0)
                )
                [
                    for (k = [0:1:n-1])
                        let(
                            v = f[k],
                            p_base = center + d_f * n_f,
                            p0_2d = _ps_rot2d(inset2d[k], handedness * angle_eff)
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
        faces_oriented = ps_orient_all_faces_outward(verts, faces_new),
        errs = [for (f = faces_oriented) _ps_face_max_plane_err(verts, f)],
        max_err = (len(errs) == 0) ? 0 : max(errs),
        bad = [for (e = errs) if (e > eps) e],
        _0 = (len(bad) > 0)
            ? echo("snub: non-planar faces", len(bad), "max_plane_err", max_err)
            : 0
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
