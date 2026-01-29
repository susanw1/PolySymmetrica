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

// Intersect 2D lines n0·x=d0 and n1·x=d1.
function _ps_line2_intersect(n0, d0, n1, d1, eps=1e-12) =
    let(det = n0[0]*n1[1] - n0[1]*n1[0])
    (abs(det) < eps) ? undef
  : [
        (d0*n1[1] - n0[1]*d1) / det,
        (n0[0]*d1 - d0*n1[0]) / det
    ];

// Inset polygon vertices using bisector-plane intersection lines per edge.
function _ps_face_inset_bisector_2d(f, fi, d_f, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0) =
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
                    d2 = v_dot(n_edge, verts0[v0]) - v_dot(n_edge, p0),
                    d2_use = (n2d_len < 1e-8) ? (v_dot(n2d_use, pts2d[k]) + d_f) : (d2 / n2d_len)
                )
                [n2d_use, d2_use]
        ],
        inset2d = [
            for (k = [0:1:n-1])
                _ps_line2_intersect(
                    lines[(k-1+n)%n][0], lines[(k-1+n)%n][1],
                    lines[k][0], lines[k][1]
                )
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

// Inradius of 2D polygon (CW), measured from origin to nearest edge line.
function _ps_face_inradius_2d(pts2d) =
    let(
        n = len(pts2d),
        ds = [
            for (i = [0:1:n-1])
                let(
                    p0 = pts2d[i],
                    p1 = pts2d[(i+1)%n],
                    e = [p1[0]-p0[0], p1[1]-p0[1]],
                    n_in = v_norm([e[1], -e[0]])
                )
                v_dot(n_in, p0)
        ]
    )
    min(ds);

// Find d_f at which bisector-based inset polygon collapses (min area).
function _ps_face_bisector_collapse_d(f, fi, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0, steps=40) =
    let(
        inr = _ps_face_inradius_2d(pts2d),
        ds = [ for (i = [0:1:steps]) inr * i / steps ],
        areas = [
            for (d = ds)
                let(poly2d = _ps_face_inset_bisector_2d(f, fi, d, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0))
                    _ps_poly_area_abs_2d(poly2d)
        ],
        min_i = _ps_index_of_min(areas)
    )
    ds[min_i];

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
        face_offsets = [
            for (fi = [0:1:len(faces0)-1])
                sum([for (j = [0:1:fi-1]) len(faces0[j])])
        ],
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

        face_offsets = [
            for (fi = [0:1:len(faces0)-1])
                sum([for (j = [0:1:fi-1]) len(faces0[j])])
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
                    inset2d = _ps_face_inset_bisector_2d(f, fi, d_f, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0),
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
