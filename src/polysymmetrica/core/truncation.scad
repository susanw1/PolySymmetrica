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

// Cantellation/expansion: offset faces and vertices, insert edge faces.
// df: face offset (outward along face normals)
// dv: vertex offset (outward along vertex directions)
// de: edge offset along edge-bisector plane (default 0)
function poly_cantellate(poly, df, dv, de = 0, eps = 1e-8, len_eps = 1e-6) =
    let(
        verts0 = poly_verts(poly),
        faces0 = orient_all_faces_outward(verts0, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = edge_faces_table(faces0, edges),

        // face planes (n·x = d)
        face_n = [ for (f = faces0) face_normal(verts0, f) ],
        face_d = [
            for (fi = [0:len(faces0)-1])
                v_dot(face_n[fi], verts0[faces0[fi][0]]) + df
        ],

        // vertex planes (n·x = d)
        vert_n = [ for (v = verts0) v_norm(v) ],
        vert_d = [ for (v = verts0) norm(v) + dv ],

        // edge planes (bisector between adjacent faces)
        edge_n = [
            for (ei = [0:len(edges)-1])
                let(fpair = edge_faces[ei])
                    v_norm(face_n[fpair[0]] + face_n[fpair[1]])
        ],
        edge_d = [
            for (ei = [0:len(edges)-1])
                let(
                    e = edges[ei],
                    mid = (verts0[e[0]] + verts0[e[1]]) / 2
                )
                v_dot(edge_n[ei], mid) + de
        ],

        // intersection point for (edge, face, vertex) triple
        edge_pts = [
            for (ei = [0:len(edges)-1])
                let(
                    e = edges[ei],
                    fpair = edge_faces[ei],
                    f0 = fpair[0],
                    f1 = fpair[1],
                    v0 = e[0],
                    v1 = e[1],

                    n_edge = edge_n[ei],
                    d_edge = edge_d[ei],

                    n_f0 = face_n[f0],
                    n_f1 = face_n[f1],
                    d_f0 = face_d[f0],
                    d_f1 = face_d[f1],

                    n_v0 = vert_n[v0],
                    n_v1 = vert_n[v1],
                    d_v0 = vert_d[v0],
                    d_v1 = vert_d[v1],

                    // solve for each combination (face, vertex)
                    p_f0_v0 = _ps_solve3([n_f0, n_v0, n_edge], [d_f0, d_v0, d_edge], eps),
                    p_f0_v1 = _ps_solve3([n_f0, n_v1, n_edge], [d_f0, d_v1, d_edge], eps),
                    p_f1_v0 = _ps_solve3([n_f1, n_v0, n_edge], [d_f1, d_v0, d_edge], eps),
                    p_f1_v1 = _ps_solve3([n_f1, n_v1, n_edge], [d_f1, d_v1, d_edge], eps)
                )
                [[p_f0_v0, p_f0_v1], [p_f1_v0, p_f1_v1]]
        ],

        _ = assert(min([
            for (ei = [0:len(edge_pts)-1])
                min([
                    for (fi = [0:1])
                        for (vi = [0:1])
                            is_undef(edge_pts[ei][fi][vi]) ? 0 : 1
                ])
        ]) == 1, "poly_cantellate: plane intersection failed"),

        // Faces from original faces
        face_faces_pts = [
            for (fi = [0:1:len(faces0)-1])
                let(
                    f = faces0[fi],
                    n = len(f)
                )
                [
                    for (k = [0:1:n-1])
                        let(
                            v_prev = f[(k-1+n)%n],
                            v = f[k],
                            v_next = f[(k+1)%n],
                            ei_prev = find_edge_index(edges, v_prev, v),
                            ei_next = find_edge_index(edges, v, v_next),
                            fpos_prev = (edge_faces[ei_prev][0] == fi) ? 0 : 1,
                            fpos_next = (edge_faces[ei_next][0] == fi) ? 0 : 1,
                            vpos_prev = (edges[ei_prev][0] == v) ? 0 : 1,
                            vpos_next = (edges[ei_next][0] == v) ? 0 : 1
                        )
                        each [
                            edge_pts[ei_prev][fpos_prev][vpos_prev],
                            edge_pts[ei_next][fpos_next][vpos_next]
                        ]
                ]
        ],

        // Faces from original vertices:
        // intersect vertex plane with the two incident edge planes per face
        vert_faces_pts = [
            for (vi = [0:1:len(verts0)-1])
                let(
                    fc = faces_around_vertex([verts0, faces0, 1], vi, edges, edge_faces)
                )
                [
                    for (fidx = fc)
                        let(
                            f = faces0[fidx],
                            m = len(f),
                            pos = [ for (k = [0:1:m-1]) if (f[k]==vi) k ][0],
                            v_prev = f[(pos-1+m)%m],
                            v_next = f[(pos+1)%m],
                            ei_prev = find_edge_index(edges, v_prev, vi),
                            ei_next = find_edge_index(edges, vi, v_next),
                            n_v = vert_n[vi],
                            d_v = vert_d[vi],
                            n_e0 = edge_n[ei_prev],
                            d_e0 = edge_d[ei_prev],
                            n_e1 = edge_n[ei_next],
                            d_e1 = edge_d[ei_next]
                        )
                        _ps_solve3([n_v, n_e0, n_e1], [d_v, d_e0, d_e1], eps)
                ]
        ],

        // Faces from original edges
        edge_faces_pts = [
            for (ei = [0:len(edges)-1])
                let(
                    e = edges[ei],
                    fpair = edge_faces[ei],
                    f0 = fpair[0],
                    f1 = fpair[1],
                    v0 = e[0],
                    v1 = e[1],
                    f0pos = 0,
                    f1pos = 1,
                    v0pos = 0,
                    v1pos = 1
                )
                [
                    edge_pts[ei][f0pos][v0pos],
                    edge_pts[ei][f0pos][v1pos],
                    edge_pts[ei][f1pos][v1pos],
                    edge_pts[ei][f1pos][v0pos]
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
