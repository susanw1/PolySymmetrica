// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <duals.scad>  // for faces_around_vertex helpers
use <transform.scad>

CHAMFER_DEBUG = true;

// --- internal helpers ---

// For edge ei=[a,b], edge_pts[ei]=[P_a,P_b] (near a and near b)
function _ps_edge_point_near(edges, edge_pts, a, b, near_v) =
    let(
        ei = ps_find_edge_index(edges, a, b),
        e  = edges[ei]
    )
    (near_v == e[0]) ? edge_pts[ei][0] : edge_pts[ei][1];

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

// Inset a 2D polygon (CW) by distance inset, using edge-offset intersections.
function _ps_face_inset_pts2d(pts2d, inset) =
    let(
        n = len(pts2d),
        edges = [
            for (i = [0:n-1])
                let(
                    p0 = pts2d[i],
                    p1 = pts2d[(i+1)%n],
                    e = [p1[0]-p0[0], p1[1]-p0[1]],
                    n_in = v_norm([e[1], -e[0]])
                )
                [n_in, v_dot(n_in, p0) + inset]
        ]
    )
    [
        for (i = [0:n-1])
            let(
                e_prev = edges[(i-1+n)%n],
                e_cur  = edges[i],
                p = _ps_line2_intersect(e_prev[0], e_prev[1], e_cur[0], e_cur[1])
            )
            p
    ];

// Inset polygon vertices using bisector-plane intersection lines per edge.
function _ps_face_inset_bisector_2d(f, fi, d_f, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0) =
    let(
        n = len(f),
        p0 = center - n_f * d_f,
        lines = [
            for (k = [0:n-1])
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
            for (k = [0:n-1])
                _ps_line2_intersect(
                    lines[(k-1+n)%n][0], lines[(k-1+n)%n][1],
                    lines[k][0], lines[k][1]
                )
        ]
    )
    inset2d;


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
            for (ei = [0:len(edges)-1])
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

        // truncated original faces -> 2n-gons
        trunc_faces_pts = [
            for (f = faces)
                let(n = len(f))
                [
                    for (k = [0:n-1])
                        let(
                            v      = f[k],
                            v_next = f[(k+1)%n],
                            v_prev = f[(k-1+n)%n],
                            P_prev = _ps_edge_point_near(edges, edge_pts, v_prev, v, v),
                            P_next = _ps_edge_point_near(edges, edge_pts, v, v_next, v)
                        )
                        // append the two points at corner v
                        // (ordering yields a consistent 2n winding when f is consistent)
                        each [P_prev, P_next]
                ]
        ],

        // vertex faces: one per original vertex
        // We use faces_around_vertex to get cyclic face order, then derive neighbor vertices
        // by intersecting consecutive faces at v; simplest is to use edges_incident_to_vertex + angular order later.
        // For now: derive from incident edges and the same "walk" logic we already have.
        edge_faces = ps_edge_faces_table(faces, edges),

        vert_faces_pts = [
            for (vi = [0:len(verts)-1])
                let(
                    // get faces around vertex in cyclic order
                    fc = faces_around_vertex(poly, vi, edges, edge_faces),

                    // turn that into an ordered list of neighbor vertices around vi:
                    // take each face in the cycle, find the two neighbors of vi in that face,
                    // and pick the "next" neighbor consistently. Bit fiddly.
                    neigh = [
                        for (idx = [0:len(fc)-1])
                            let(
                                f = faces[ fc[idx] ],
                                m = len(f),
                                pos = [ for (k=[0:m-1]) if (f[k]==vi) k ][0],
                                v_next = f[(pos+1)%m]   // neighbor after vi in that face winding
                            )
                            v_next
                    ],

                    // now map each neighbor to the trunc point on edge (vi, neighbor) near vi
                    pts = [
                        for (vn = neigh)
                            _ps_edge_point_near(edges, edge_pts, vi, vn, vi)
                    ]
                )
                pts
        ],

        // collect all face point lists
        faces_pts_all = concat(trunc_faces_pts, vert_faces_pts)
    )
    _ps_poly_from_face_points(faces_pts_all, eps);

// Rectification: replace each vertex with the midpoint of each incident edge.
function poly_rectify(poly) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        edges = _ps_edges_from_faces(faces),

        edge_mid = [
            for (e = edges)
                (verts[e[0]] + verts[e[1]]) / 2
        ],

        // Faces corresponding to original faces: walk edges in face order.
        face_faces = [
            for (f = faces)
                let(n = len(f))
                [
                    for (k = [0:n-1])
                        let(
                            a = f[k],
                            b = f[(k+1)%n],
                            ei = ps_find_edge_index(edges, a, b)
                        )
                        ei
                ]
        ],

        // Faces corresponding to original vertices: cycle around the vertex.
        edge_faces = ps_edge_faces_table(faces, edges),
        vert_faces = [
            for (vi = [0:len(verts)-1])
                let(
                    fc = faces_around_vertex(poly, vi, edges, edge_faces),
                    neigh = [
                        for (idx = [0:len(fc)-1])
                            let(
                                f = faces[fc[idx]],
                                m = len(f),
                                pos = [ for (k = [0:m-1]) if (f[k]==vi) k ][0],
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
        faces_out = ps_orient_all_faces_outward(edge_mid, faces_idx),

        edges_new = _ps_edges_from_faces(faces_out),
        e0 = edges_new[0],
        vA = edge_mid[e0[0]],
        vB = edge_mid[e0[1]],
        unit_e = norm(vB - vA),
        mid = (vA + vB) / 2,
        ir  = norm(mid),
        e_over_ir = unit_e / ir
    )
    make_poly(edge_mid / unit_e, faces_out, e_over_ir);

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

function _ps_vert_faces_from_offsets(verts0, faces0, face_pts, edges, edge_faces) =
    [
        for (vi = [0:1:len(verts0)-1])
            let(
                fc = faces_around_vertex([verts0, faces0, 1], vi, edges, edge_faces),
                idxs = [for (k = [0:1:len(fc)-1]) k]
            )
            [
                for (k = idxs)
                    let(
                        fi = fc[k],
                        f = faces0[fi],
                        pos = [ for (j=[0:1:len(f)-1]) if (f[j]==vi) j ][0]
                    )
                    face_pts[fi][pos]
            ]
    ];

function _ps_edge_faces_from_offsets(faces0, face_pts, edges, edge_faces) =
    [
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
                face_pts[f0][pos0],
                face_pts[f0][pos1],
                face_pts[f1][pos2],
                face_pts[f1][pos3]
            ]
    ];

function poly_cantellate(poly, df, eps = 1e-8, len_eps = 1e-6) =
    let(
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),

        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],
        // offset face corners: one point per (face, vertex) incidence
        face_pts = _ps_face_offset_pts(verts0, faces0, face_n, df),
        // Faces from original faces (n-gons)
        face_faces_pts = face_pts,
        // Faces from original edges (quads)
        edge_faces_pts = _ps_edge_faces_from_offsets(faces0, face_pts, edges, edge_faces),
        // Faces from original vertices (valence-gons)
        vert_faces_pts = _ps_vert_faces_from_offsets(verts0, faces0, face_pts, edges, edge_faces),

        faces_pts_all = concat(face_faces_pts, edge_faces_pts, vert_faces_pts)
    )
    _ps_poly_from_face_points(faces_pts_all, eps, len_eps);

// Chamfer: face faces + edge faces (no vertex faces).
// t is interpreted as an inset fraction of mean edge length per face.
function poly_chamfer(poly, t, eps = 1e-8, len_eps = 1e-6) =
    let(
        t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t,
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],
        debug = (!is_undef(CHAMFER_DEBUG) && CHAMFER_DEBUG),

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
                        for (k = [0:n-1])
                            let(p = verts0[f[k]] - center)
                                [v_dot(p, ex), v_dot(p, ey)]
                    ],
                    edge_lens = [
                        for (k = [0:n-1])
                            let(p0 = pts2d[k], p1 = pts2d[(k+1)%n])
                                norm([p1[0]-p0[0], p1[1]-p0[1]])
                    ],
                    d_f0 = t_eff * (sum(edge_lens) / n),
                    inset2d_a = _ps_face_inset_bisector_2d(f, fi, d_f0, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0),
                    r0 = sum([for (p = pts2d) norm(p)]) / n,
                    r1 = sum([for (p = inset2d_a) norm(p)]) / n,
                    d_f = (r1 > r0) ? -d_f0 : d_f0,
                    inset2d = (r1 > r0)
                        ? _ps_face_inset_bisector_2d(f, fi, d_f, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0)
                        : inset2d_a,
                    p0 = center - n_f * d_f
                )
                [
                    for (k = [0:n-1])
                        p0 + ex * inset2d[k][0] + ey * inset2d[k][1]
                ]
        ],

        sites = [
            for (fi = [0:1:len(faces0)-1])
                let(f = faces0[fi], n = len(f))
                for (k = [0:n-1])
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
                [ for (k = [0:n-1]) [1, face_offsets[fi] + k] ]
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

        cycles_all = concat(face_cycles, edge_cycles),
        _dbg = debug ? echo("chamfer_debug edge0", edge_cycles[0]) : 0,
        _dbg2 = debug ? echo("chamfer_debug face0", face_cycles[0]) : 0,
        _dbg3 = debug ? echo("chamfer_debug face0_verts", faces0[0]) : 0,
        _dbg4 = debug ? echo("chamfer_debug face0_pts3d", face_pts3d[0]) : 0,
        _dbg5 = debug ? echo("chamfer_debug edge0_e", edges[0], "edge0_faces", edge_faces[0]) : 0
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
