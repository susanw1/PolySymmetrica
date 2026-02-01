// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <duals.scad>  // for faces_around_vertex helpers
use <transform.scad>

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

// Index k where edge f[k]->f[k+1] matches (v0->v1), or -1 if not found.
function _ps_face_edge_index(f, v0, v1) =
    let(
        n = len(f),
        hits = [for (k = [0:1:n-1]) if (f[k] == v0 && f[(k+1)%n] == v1) k]
    )
    (len(hits) == 0) ? -1 : hits[0];

function _ps_clamp(x, lo, hi) = (x < lo) ? lo : (x > hi) ? hi : x;

function _ps_face_edge_site(base, k, near_next=false) =
    base + 2*k + (near_next ? 1 : 0);

// Return [s_near_v0, s_near_v1] for edge (v0,v1) on face fidx.
function _ps_face_edge_sites_for_face_edge(faces, fidx, v0, v1, base) =
    let(
        f = faces[fidx],
        k_dir = _ps_face_edge_index(f, v0, v1),
        k = (k_dir >= 0) ? k_dir : _ps_face_edge_index(f, v1, v0),
        flip = (k_dir < 0)
    )
    assert(k >= 0, "cantitruncate: edge not found in face")
    [
        flip ? _ps_face_edge_site(base, k, true) : _ps_face_edge_site(base, k, false),
        flip ? _ps_face_edge_site(base, k, false) : _ps_face_edge_site(base, k, true)
    ];

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

function _ps_map_face_c(face_len, c_by_size, default_c=0) =
    let(
        idxs = [for (i = [0:1:len(c_by_size)-1]) if (c_by_size[i][0] == face_len) i]
    )
    (len(idxs) == 0) ? default_c : c_by_size[idxs[0]][1];

function _ps_map_edge_c(face_len, adj_len, c_by_pair, default_c=0) =
    let(
        a = min(face_len, adj_len),
        b = max(face_len, adj_len),
        idxs = [for (i = [0:1:len(c_by_pair)-1]) if (c_by_pair[i][0] == a && c_by_pair[i][1] == b) i]
    )
    (len(idxs) == 0) ? default_c : c_by_pair[idxs[0]][2];


// Intersect 2D lines n0·x=d0 and n1·x=d1.
function _ps_line2_intersect(n0, d0, n1, d1, eps=1e-12) =
    let(det = n0[0]*n1[1] - n0[1]*n1[0])
    (abs(det) < eps) ? undef
  : [
        (d0*n1[1] - n0[1]*d1) / det,
        (n0[0]*d1 - d0*n1[0]) / det
    ];

// Inset polygon vertices using bisector-plane intersection lines per edge.
// d_f shifts the face plane; d_e shifts the edge-bisector planes along their normals.
function _ps_face_inset_bisector_2d(f, fi, d_f, d_e, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0) =
    let(
        n = len(f),
        p0 = center - n_f * d_f,
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
                    d2 = v_dot(n_edge, verts0[v0]) - v_dot(n_edge, p0) + d_e,
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

// Same as _ps_face_inset_bisector_2d but uses per-edge d_e values.
function _ps_face_inset_bisector_2d_edges(f, fi, d_f, d_e_edges, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0) =
    let(
        n = len(f),
        p0 = center - n_f * d_f,
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
                    d2 = v_dot(n_edge, verts0[v0]) - v_dot(n_edge, p0) + d_e_edges[k],
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

function poly_truncate(poly, t, eps = 1e-8) =
    let(t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t)
    // should this assert? Or should we allow, or silently call rectify?
    assert(t_eff != 0.5, "'t' cannot be 0.5 as this produces degenerate vertices - use poly_rectify() instead")
    (t_eff == 0) 
        ? poly
        :
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        edges = _ps_edges_from_faces(faces),

        // edge points: aligned with edges[ei]=[a,b]
        edge_pts = [
            for (ei = [0:1:len(edges)-1])
                let(
                    a = edges[ei][0],
                    b = edges[ei][1],
                    A = verts[a],
                    B = verts[b],
                    P_a = A + t_eff*(B-A),
                    P_b = B + t_eff*(A-B)
                )
                [P_a, P_b]
        ],

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
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        faces0 = ps_orient_all_faces_outward(verts, faces),
        poly0 = make_poly(verts, faces0, poly_e_over_ir(poly)),
        edges = _ps_edges_from_faces(faces0),

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
        edge_faces = ps_edge_faces_table(faces0, edges),
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

function poly_cantellate(poly, df, eps = 1e-8, len_eps = 1e-6) =
    let(
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        poly0 = make_poly(verts0, faces0, poly_e_over_ir(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),

        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],
        // offset face corners: one point per (face, vertex) incidence
        face_pts = _ps_face_offset_pts(verts0, faces0, face_n, df),
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
                let(
                    fc = faces_around_vertex(poly0, vi, edges, edge_faces)
                )
                [
                    for (fi = fc)
                        let(pos = _ps_index_of(faces0[fi], vi))
                            [1, face_offsets[fi] + pos]
                ]
        ],
        cycles_all = concat(face_cycles, edge_cycles, vert_cycles)
    )
    ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps);

// Chamfer: face faces + edge faces (no vertex faces).
// t is a signed face-plane offset as a fraction of each face's collapse distance.
// Positive t = inward chamfer; negative t = anti-chamfer.
function poly_chamfer(poly, t, eps = 1e-8, len_eps = 1e-6) =
    let(
        t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t,
        _t_ok = assert(abs(t_eff) != 1, "poly_chamfer: |t| must not be 1 (t=±1 collapses faces; use cleanup/hyper mode)"),
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],

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

// Cantitruncation: truncation + cantellation (two parameters).
// t controls face-plane shift (like chamfer), c controls edge/vertex expansion (like cantellate).
function poly_cantitruncate(poly, t, c, eps = 1e-8, len_eps = 1e-6) =
    let(
        t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t,
        c_eff = is_undef(c) ? 0 : c,
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        poly0 = make_poly(verts0, faces0, poly_e_over_ir(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],
        // scale by inter-radius for consistent parameterization
        ir = min([for (e = edges) norm((verts0[e[0]] + verts0[e[1]]) / 2)]),
        d_f = -c_eff * ir,
        d_e = c_eff * ir,
        face_offsets = _ps_face_offsets(faces0),
        // edge points (truncation-style): two per edge
        edge_pts = [
            for (ei = [0:1:len(edges)-1])
                let(
                    a = edges[ei][0],
                    b = edges[ei][1],
                    A = verts0[a],
                    B = verts0[b]
                )
                [ A + t_eff * (B - A), B + t_eff * (A - B) ]
        ],
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
        face_cycles = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    base = edge_site_offset + face_edge_offsets[fi]
                )
                [
                    for (k = [0:1:n-1])
                        let(
                            s0 = _ps_face_edge_site(base, k, false),
                            s1 = _ps_face_edge_site(base, k, true)
                        )
                        each [[1, s0], [1, s1]]
                ]
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
                    b0 = edge_site_offset + face_edge_offsets[f0],
                    b1 = edge_site_offset + face_edge_offsets[f1],
                    s_pair0 = _ps_face_edge_sites_for_face_edge(faces0, f0, v0, v1, b0),
                    s_pair1 = _ps_face_edge_sites_for_face_edge(faces0, f1, v1, v0, b1),
                    s0 = s_pair0[0],
                    s1 = s_pair0[1],
                    s2 = s_pair1[0],
                    s3 = s_pair1[1]
                )
                [
                    [1, s0],
                    [1, s1],
                    [1, s2],
                    [1, s3]
                ]
        ],
        vert_cycles = [
            for (vi = [0:1:len(verts0)-1])
                let(fc = faces_around_vertex(poly0, vi, edges, edge_faces))
                [
                    for (fi = fc)
                        let(
                            f = faces0[fi],
                            n = len(f),
                            pos = _ps_index_of(f, vi),
                            k_prev = (pos - 1 + n) % n,
                            base = edge_site_offset + face_edge_offsets[fi],
                            s_prev = _ps_face_edge_site(base, k_prev, true),
                            s_next = _ps_face_edge_site(base, pos, false)
                        )
                        each [[1, s_prev], [1, s_next]]
                ]
        ],
        // Debug: echo one decagon face cycle points and their angles (first face with n>=5*2)
        cycles_all = concat(face_cycles, edge_cycles, vert_cycles)
    )
    ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps);

// Cantitruncation with per-face-family c values (indexed by face size).
// c_by_size: list of [face_size, c] pairs; default_c used if size not found.
function poly_cantitruncate_families(poly, t, c_by_size, default_c=0, c_edge_by_pair=undef, eps=1e-8, len_eps=1e-6) =
    let(
        t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t,
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        poly0 = make_poly(verts0, faces0, poly_e_over_ir(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],
        ir = min([for (e = edges) norm((verts0[e[0]] + verts0[e[1]]) / 2)]),
        face_offsets = _ps_face_offsets(faces0),
        // edge points (truncation-style): two per edge
        edge_pts = [
            for (ei = [0:1:len(edges)-1])
                let(a = edges[ei][0], b = edges[ei][1], A = verts0[a], B = verts0[b])
                [ A + t_eff * (B - A), B + t_eff * (A - B) ]
        ],
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
                    inset2d = _ps_face_inset_bisector_2d_edges(f, fi, d_f, d_e_edges, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0)
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
        face_cycles = [
            for (fi = [0:1:len(faces0)-1])
                let(n = len(faces0[fi]), base = edge_site_offset + face_edge_offsets[fi])
                [
                    for (k = [0:1:n-1])
                        let(
                            s0 = _ps_face_edge_site(base, k, false),
                            s1 = _ps_face_edge_site(base, k, true)
                        )
                        each [[1, s0], [1, s1]]
                ]
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
                    b0 = edge_site_offset + face_edge_offsets[f0],
                    b1 = edge_site_offset + face_edge_offsets[f1],
                    s_pair0 = _ps_face_edge_sites_for_face_edge(faces0, f0, v0, v1, b0),
                    s_pair1 = _ps_face_edge_sites_for_face_edge(faces0, f1, v1, v0, b1),
                    s0 = s_pair0[0],
                    s1 = s_pair0[1],
                    s2 = s_pair1[0],
                    s3 = s_pair1[1]
                )
                [
                    [1, s0],
                    [1, s1],
                    [1, s2],
                    [1, s3]
                ]
        ],
        vert_cycles = [
            for (vi = [0:1:len(verts0)-1])
                let(fc = faces_around_vertex(poly0, vi, edges, edge_faces))
                [
                    for (fi = fc)
                        let(
                            f = faces0[fi],
                            n = len(f),
                            pos = _ps_index_of(f, vi),
                            k_prev = (pos - 1 + n) % n,
                            base = edge_site_offset + face_edge_offsets[fi],
                            s_prev = _ps_face_edge_site(base, k_prev, true),
                            s_next = _ps_face_edge_site(base, pos, false)
                        )
                        each [[1, s_prev], [1, s_next]]
                ]
        ],
        cycles_all = concat(face_cycles, edge_cycles, vert_cycles)
    )
    ps_poly_transform_from_sites(verts0, sites, site_points, cycles_all, eps, len_eps);

// Dominant-family trig solver (no grid search). Returns [t, c_by_size].
function solve_cantitruncate_dominant(poly, dominant_size, edge_idx=undef) =
    let(
        verts = poly_verts(poly),
        faces = ps_orient_all_faces_outward(verts, poly_faces(poly)),
        sizes = [for (f = faces) len(f)],
        dom_face_idx = [for (i = [0:1:len(faces)-1]) if (sizes[i] == dominant_size) i][0],
        f = faces[dom_face_idx],
        n = len(f),
        v0 = f[0],
        v_prev = f[(n-1)%n],
        v_next = f[1],
        a0 = v_norm(verts[v_prev] - verts[v0]),
        a1 = v_norm(verts[v_next] - verts[v0]),
        phi = acos(_ps_clamp(v_dot(a0, a1), -1, 1)),
        t = 1 / (2 * (1 + sin(phi/2))),
        edges = _ps_edges_from_faces(faces),
        edge_faces = ps_edge_faces_table(faces, edges),
        ir = min([for (e = edges) norm((verts[e[0]] + verts[e[1]]) / 2)]),
        a = norm(verts[v_next] - verts[v0]),
        // compute target sums for each secondary family via edges with dominant
        fam_sizes = _ps_sort([for (s = sizes) s]),
        uniq_sizes = [for (i = [0:1:len(fam_sizes)-1]) if (i==0 || fam_sizes[i] != fam_sizes[i-1]) fam_sizes[i]],
        other_sizes = [for (s = uniq_sizes) if (s != dominant_size) s],
        targets = [
            for (s = other_sizes)
                let(
                    edges_s = [
                        for (ei = [0:1:len(edges)-1])
                            let(
                                fpair = edge_faces[ei],
                                s0 = len(faces[fpair[0]]),
                                s1 = len(faces[fpair[1]])
                            )
                            if ((s0 == dominant_size && s1 == s) || (s1 == dominant_size && s0 == s)) ei
                    ],
                    vals = [
                        for (ei = edges_s)
                            let(
                                fpair = edge_faces[ei],
                                n0 = ps_face_normal(verts, faces[fpair[0]]),
                                n1 = ps_face_normal(verts, faces[fpair[1]]),
                                alpha = acos(_ps_clamp(v_dot(n0, n1), -1, 1))
                            )
                            (1 - 2*t) * a / (2 * sin(alpha/2))
                    ]
                )
                [s, (len(vals) == 0) ? undef : sum(vals)/len(vals)]
        ],
        // choose reference family (first with target) to set d_f_dom
        has_other = (len(targets) > 0),
        ref_idx = has_other ? [for (i = [0:1:len(targets)-1]) if (!is_undef(targets[i][1])) i][0] : undef,
        ref_target = is_undef(ref_idx) ? undef : targets[ref_idx][1],
        d_f_dom = is_undef(ref_target) ? 0 : (ref_target / 2),
        c_by_size = concat(
            [[dominant_size, abs(d_f_dom)/ir]],
            [for (tgt = targets)
                let(
                    s = tgt[0],
                    targ = tgt[1],
                    d_f_s = is_undef(targ) ? d_f_dom : (targ - d_f_dom)
                )
                [s, abs(d_f_s)/ir]
            ]
        )
    )
    [t, c_by_size];

// Build unique [a,b] pairs from a list of pairs (order is preserved).
function _ps_unique_pairs(pairs, i=0, acc=[]) =
    (i >= len(pairs)) ? acc
  : let(
        p = pairs[i],
        has = (len([for (u = acc) if (u[0] == p[0] && u[1] == p[1]) 1]) > 0)
    )
    _ps_unique_pairs(pairs, i+1, has ? acc : concat(acc, [p]));

// Dominant-family trig solver with per-edge-family c values.
// Returns [t, c_by_size, c_edge_by_pair].
function solve_cantitruncate_dominant_edges(poly, dominant_size, edge_idx=undef) =
    let(
        sol = solve_cantitruncate_dominant(poly, dominant_size, edge_idx),
        t = sol[0],
        c_by_size = sol[1],
        verts = poly_verts(poly),
        faces = ps_orient_all_faces_outward(verts, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces),
        edge_faces = ps_edge_faces_table(faces, edges),
        pair_keys = _ps_unique_pairs([
            for (ei = [0:1:len(edges)-1])
                let(
                    fpair = edge_faces[ei],
                    s0 = len(faces[fpair[0]]),
                    s1 = len(faces[fpair[1]])
                )
                [min(s0, s1), max(s0, s1)]
        ]),
        c_edge_by_pair = [
            for (p = pair_keys)
                let(
                    c0 = _ps_map_face_c(p[0], c_by_size, 0),
                    c1 = _ps_map_face_c(p[1], c_by_size, 0),
                    c_edge = (c0 + c1) / 2
                )
                [p[0], p[1], c_edge]
        ]
    )
    [t, c_by_size, c_edge_by_pair];

// Measure how square an edge face is (edge length spread).
function _ps_face_edge_spread(verts, face) =
    let(
        n = len(face),
        ls = (n < 2) ? [] : [for (i = [0:1:n-1]) norm(verts[face[i]] - verts[face[(i+1)%n]])]
    )
    (n == 4) ? (max(ls) - min(ls)) : undef;

function _ps_edge_len_spread(poly) =
    let(
        verts = poly_verts(poly),
        edges = _ps_edges_from_faces(poly_faces(poly)),
        ls = [for (e = edges) norm(verts[e[0]] - verts[e[1]])]
    )
    (len(ls) == 0) ? 1e30 : (max(ls) - min(ls));

function _ps_square_face_spread(poly, face_k=4) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        spreads = [ for (f = faces) if (len(f) == face_k) _ps_face_edge_spread(verts, f) ]
    )
    (len(spreads) == 0) ? 1e30 : max(spreads);

// Fast metrics for cantitruncate: measure cycle edge lengths without building full poly.
function _ps_cantitruncate_edge_metrics(poly, t, c) =
    let(
        t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t,
        c_eff = is_undef(c) ? 0 : c,
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],
        ir = min([for (e = edges) norm((verts0[e[0]] + verts0[e[1]]) / 2)]),
        d_f = -c_eff * ir,
        d_e = c_eff * ir,
        // edge points (truncation-style): two per edge
        edge_pts = [
            for (ei = [0:1:len(edges)-1])
                let(
                    a = edges[ei][0],
                    b = edges[ei][1],
                    A = verts0[a],
                    B = verts0[b]
                )
                [ A + t_eff * (B - A), B + t_eff * (A - B) ]
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
        face_offsets = _ps_face_offsets(faces0),
        face_pts_flat = [ for (fi = [0:1:len(faces0)-1]) for (p = face_pts3d[fi]) p ],
        face_edge_pts_flat = [ for (fi = [0:1:len(faces0)-1]) for (p = face_edge_pts3d[fi]) p ],
        edge_site_offset = 0,
        face_site_offset = len(face_edge_pts_flat),
        // cycles: face cycles (edge points), edge cycles (face points), vert cycles (face points)
        face_cycles = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f),
                    base = edge_site_offset + face_edge_offsets[fi]
                )
                [
                    for (k = [0:1:n-1])
                        let(
                            s0 = _ps_face_edge_site(base, k, false),
                            s1 = _ps_face_edge_site(base, k, true)
                        )
                        each [[1, s0], [1, s1]]
                ]
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
                    [1, face_site_offset + face_offsets[f0] + pos0],
                    [1, face_site_offset + face_offsets[f0] + pos1],
                    [1, face_site_offset + face_offsets[f1] + pos2],
                    [1, face_site_offset + face_offsets[f1] + pos3]
                ]
        ],
        vert_cycles = [
            for (vi = [0:1:len(verts0)-1])
                let(fc = faces_around_vertex(make_poly(verts0, faces0, poly_e_over_ir(poly)), vi, edges, edge_faces))
                [
                    for (fi = fc)
                        let(pos = _ps_index_of(faces0[fi], vi))
                            [1, face_site_offset + face_offsets[fi] + pos]
                ]
        ],
        // point lookup for cycle entries (only site points in this metric)
        _pt = function(c)
            (c[0] == 0) ? verts0[c[1]] :
            (c[1] < face_site_offset) ? face_edge_pts_flat[c[1]] : face_pts_flat[c[1] - face_site_offset],
        cycle_lengths = [
            for (cy = concat(face_cycles, edge_cycles, vert_cycles))
                let(
                    m = len(cy),
                    pts = [ for (c = cy) _pt(c) ]
                )
                [
                    for (i = [0:1:m-1])
                        let(p0 = pts[i], p1 = pts[(i+1)%m])
                        norm(p1 - p0)
                ]
        ],
        all_lengths = [ for (ls = cycle_lengths) for (l = ls) l ],
        min_len = (len(all_lengths) == 0) ? 0 : min(all_lengths),
        deg_penalty = (min_len < 1e-6) ? 1e6 : 0,
        edge_spread = (len(all_lengths) == 0) ? 1e30 : (max(all_lengths) - min(all_lengths)) + deg_penalty,
        square_spread = (len(edge_cycles) == 0) ? 1e30 :
            max([for (ls = [for (cy = edge_cycles) [ for (i = [0:1:len(cy)-1]) let(p0=_pt(cy[i]), p1=_pt(cy[(i+1)%len(cy)])) norm(p1-p0) ] ]) max(ls) - min(ls)]) + deg_penalty
    )
    [edge_spread, square_spread];

// Solve t,c to minimize edge length spread and square-face spread.
function solve_cantitruncate_uniform(poly, t_min=0, t_max=1, c_min=0, c_max=1,
                                     steps_t=12, steps_c=12, square_face_k=4,
                                     w_edge=1, w_square=1, eps=1e-9, len_eps=1e-6) =
    let(
        ts = [ for (i = [0:1:steps_t]) t_min + (t_max - t_min) * i / steps_t ],
        cs = [ for (j = [0:1:steps_c]) c_min + (c_max - c_min) * j / steps_c ],
        cand = [
            for (t = ts)
                for (c = cs)
                    let(
                        q = poly_cantitruncate(poly, t, c, eps, len_eps),
                        e_spread = _ps_edge_len_spread(q),
                        s_spread = _ps_square_face_spread(q, square_face_k),
                        score = w_edge * e_spread + w_square * s_spread
                    )
                    [score, t, c]
        ],
        best = cand[_ps_index_of_min([for (x = cand) x[0]])]
    )
    [best[1], best[2], best[0]];

// Fast solver: uses edge-face quads only (no full poly build).
function solve_cantitruncate_uniform_fast(poly, t_min=0, t_max=1, c_min=0, c_max=1,
                                          steps_t=12, steps_c=12, w_edge=1, w_square=1) =
    let(
        ts = [ for (i = [0:1:steps_t]) t_min + (t_max - t_min) * i / steps_t ],
        cs = [ for (j = [0:1:steps_c]) c_min + (c_max - c_min) * j / steps_c ],
        cand = [
            for (t = ts)
                for (c = cs)
                    let(
                        m = _ps_cantitruncate_edge_metrics(poly, t, c),
                        e_spread = m[0],
                        s_spread = m[1],
                        score = w_edge * e_spread + w_square * s_spread
                    )
                    [score, t, c]
        ],
        best = cand[_ps_index_of_min([for (x = cand) x[0]])]
    )
    [best[1], best[2], best[0]];

// Fast solver with coarse-to-fine refinement.
function solve_cantitruncate_uniform_fast_refine(poly, t_min=0, t_max=1, c_min=0, c_max=1,
                                                 steps=6, rounds=3, w_edge=1, w_square=1) =
    let(
        t0 = t_min, t1 = t_max,
        c0 = c_min, c1 = c_max,
        t_span = t1 - t0,
        c_span = c1 - c0,
        _solve_round = function(t_lo, t_hi, c_lo, c_hi, span_t, span_c, r) 
            let(
                sol = solve_cantitruncate_uniform_fast(poly, t_lo, t_hi, c_lo, c_hi, steps, steps, w_edge, w_square),
                t_best = sol[0],
                c_best = sol[1],
                dt = span_t / steps,
                dc = span_c / steps,
                t_lo2 = max(t_min, t_best - dt),
                t_hi2 = min(t_max, t_best + dt),
                c_lo2 = max(c_min, c_best - dc),
                c_hi2 = min(c_max, c_best + dc),
                span_t2 = t_hi2 - t_lo2,
                span_c2 = c_hi2 - c_lo2
            )
            (r <= 1) ? sol : _solve_round(t_lo2, t_hi2, c_lo2, c_hi2, span_t2, span_c2, r - 1)
    )
    _solve_round(t0, t1, c0, c1, t_span, c_span, rounds);

// Trig-based solver for regular bases (one edge type).
// Uses face interior angle and dihedral to compute t and c directly.
function solve_cantitruncate_trig(poly, face_idx=0, edge_idx=undef) =
    let(
        verts = poly_verts(poly),
        faces = ps_orient_all_faces_outward(verts, poly_faces(poly)),
        f = faces[face_idx],
        n = len(f),
        v0 = f[0],
        v_prev = f[(n-1)%n],
        v_next = f[1],
        a0 = v_norm(verts[v_prev] - verts[v0]),
        a1 = v_norm(verts[v_next] - verts[v0]),
        phi = acos(_ps_clamp(v_dot(a0, a1), -1, 1)),
        t = 1 / (2 * (1 + sin(phi/2))),
        edges = _ps_edges_from_faces(faces),
        ei = is_undef(edge_idx) ? ps_find_edge_index(edges, v0, v_next) : edge_idx,
        fpair = ps_edge_faces_table(faces, edges)[ei],
        f0 = fpair[0],
        f1 = fpair[1],
        n0 = ps_face_normal(verts, faces[f0]),
        n1 = ps_face_normal(verts, faces[f1]),
        // alpha is angle between outward normals
        alpha = acos(_ps_clamp(v_dot(n0, n1), -1, 1)),
        a = norm(verts[v_next] - verts[v0]),
        ir = min([for (e = edges) norm((verts[e[0]] + verts[e[1]]) / 2)]),
        // across-face distance uses sin(alpha/2); for cube alpha=90 so sin=cos
        d_f = (1 - 2*t) * a / (2 * sin(alpha/2)),
        c = abs(d_f) / ir
    )
    [t, c];

function poly_cantitruncate_uniform(poly, t_min=0, t_max=1, c_min=0, c_max=1,
                                   steps_t=12, steps_c=12, square_face_k=4,
                                   w_edge=1, w_square=1, eps=1e-9, len_eps=1e-6) =
    let(sol = solve_cantitruncate_uniform(poly, t_min, t_max, c_min, c_max,
                                          steps_t, steps_c, square_face_k,
                                          w_edge, w_square, eps, len_eps))
    poly_cantitruncate(poly, sol[0], sol[1], eps, len_eps);


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
    let(
        c0 = (c < 0) ? 0 : ((c > 1) ? 1 : c),
        verts = poly_verts(poly),
        edges = _ps_edges_from_faces(poly_faces(poly)),
        e0 = edges[0],
        edge_len = norm(verts[e0[1]] - verts[e0[0]]),
        ir = edge_len / poly_e_over_ir(poly),
        df_max_eff = is_undef(df_max) ? (2 * ir) : df_max,
        df_mid = cantellate_square_df(poly, 0, df_max_eff, steps, family_edge_idx),
        df = (c0 <= 0.5)
            ? (2 * c0 * df_mid)
            : (df_mid + (c0 - 0.5) * 2 * (df_max_eff - df_mid))
    )
    poly_cantellate(poly, df, eps, len_eps);




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
