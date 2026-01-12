// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <duals.scad>  // for faces_around_vertex helpers

// --- internal helpers ---

// For edge ei=[a,b], edge_pts[ei]=[P_a,P_b] (near a and near b)
function _ps_edge_point_near(edges, edge_pts, a, b, near_v) =
    let(
        ei = find_edge_index(edges, a, b),
        e  = edges[ei]
    )
    (near_v == e[0]) ? edge_pts[ei][0] : edge_pts[ei][1];


// index of point p in list (or -1)
function _ps_find_point(list, p, eps, i=0) =
    (i >= len(list)) ? -1 :
    (point_eq(list[i], p, eps) ? i : _ps_find_point(list, p, eps, i+1));

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

// --- main truncation ---

function poly_truncate(poly, t, eps = 1e-8) =
    let(t_eff = is_undef(t) ? _ps_truncate_default_t(poly) : t)
    assert(t_eff >= 0 && t_eff < 0.5, "'t' out of range")
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
        edge_faces = edge_faces_table(faces, edges),

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
        faces_pts_all = concat(trunc_faces_pts, vert_faces_pts),

        // flatten all points for uniq build
        all_pts = [ for (fp = faces_pts_all) for (p = fp) p ],
        uniq_verts = _ps_unique_points(all_pts, eps),

        // remap faces
        faces_idx = [ for (fp = faces_pts_all) _ps_face_points_to_indices(uniq_verts, fp, eps) ],

        // orient outward on the constructed geometry
        faces_out = orient_all_faces_outward(uniq_verts, faces_idx),

        // compute unit_edge and e_over_ir from first edge
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
                            ei = find_edge_index(edges, a, b)
                        )
                        ei
                ]
        ],

        // Faces corresponding to original vertices: cycle around the vertex.
        edge_faces = edge_faces_table(faces, edges),
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
                        let(ei = find_edge_index(edges, vi, vn))
                        ei
                ]
        ],

        faces_idx = concat(face_faces, vert_faces),
        faces_out = orient_all_faces_outward(edge_mid, faces_idx),

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
function poly_cantellate(poly, df, eps = 1e-8, len_eps = 1e-6) =
    let(
        verts0 = poly_verts(poly),
        faces0 = orient_all_faces_outward(verts0, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = edge_faces_table(faces0, edges),

        face_n = [ for (f = faces0) face_normal(verts0, f) ],
        // offset face corners: one point per (face, vertex) incidence
        face_pts = [
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
        ],

        // Faces from original faces (n-gons)
        face_faces_pts = face_pts,

        // Faces from original vertices (valence-gons)
        vert_faces_pts = [
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
        ],

        // Faces from original edges (quads): connect corresponding offset edges
        edge_faces_pts = [
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
        ],

        faces_pts_all = concat(face_faces_pts, edge_faces_pts, vert_faces_pts),
        all_pts = [ for (fp = faces_pts_all) for (p = fp) p ],
        uniq_verts = _ps_unique_points(all_pts, len_eps),
        faces_idx = [ for (fp = faces_pts_all) _ps_face_points_to_indices(uniq_verts, fp, len_eps) ],
        faces_out = orient_all_faces_outward(uniq_verts, faces_idx),

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
        faces0 = orient_all_faces_outward(poly_verts(poly), poly_faces(poly)),
        edges0 = _ps_edges_from_faces(faces0),
        edge_faces0 = edge_faces_table(faces0, edges0),
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
