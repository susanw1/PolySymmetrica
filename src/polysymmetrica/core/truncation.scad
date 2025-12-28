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
