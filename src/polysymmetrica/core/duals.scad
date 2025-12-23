// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <placement.scad>


//function poly_face_centers_unscaled(poly) =
//    let(
//        faces = poly_faces(poly),
//        verts = poly_verts(poly)
//    )
//    [ for (fi = [0 : len(faces)-1])
//        face_centroid(verts, faces[fi])
//    ];


function vertex_incident_faces(poly, vi) =
    let(faces = poly_faces(poly))
    [
        for (fi = [0 : len(faces)-1])
            let(f = faces[fi])
            if (sum([
                for (k = [0 : len(f)-1])
                    f[k] == vi ? 1 : 0
            ]) > 0)
                fi
    ];


//function edges_incident_to_vertex(edges, v) =
//    [
//        for (ei = [0 : 1 : len(edges)-1])
//            let(e = edges[ei])
//            if (e[0] == v || e[1] == v) ei
//    ];


function find_edge_index(edges, a, b) =
    let(
        e = (a < b) ? [a,b] : [b,a],
        idxs = [for (i = [0 : len(edges)-1]) if (edge_equal(edges[i], e)) i]
    )
    idxs[0];   // assume the edge exists


function next_face_around_vertex(v, f_cur, f_prev, faces, edges, edge_faces) =
    let(
        f = faces[f_cur],
        n = len(f),

        // position of v in this face
        pos = [for (k = [0 : n-1]) if (f[k] == v) k],
        k0  = pos[0],                 // there should be exactly one

        // the two neighbours of v in this face
        k_prev = (k0 - 1 + n) % n,
        k_next = (k0 + 1) % n,
        v_prev = f[k_prev],
        v_next = f[k_next],

        // edges (v -> v_next) and (v_prev -> v)
        ei1 = find_edge_index(edges, v,      v_next),
        ei2 = find_edge_index(edges, v_prev, v     ),

        ef1 = edge_faces[ei1],
        ef2 = edge_faces[ei2],

        cand1 = (ef1[0] == f_cur ? ef1[1] : ef1[0]),
        cand2 = (ef2[0] == f_cur ? ef2[1] : ef2[0]),

        candidates = [cand1, cand2],
        filtered   = [for (cf = candidates) if (cf != f_prev) cf]
    )
    filtered[0];   // in a convex poly this is unique


function faces_around_vertex_rec(v, f_cur, f_prev, f_start,
                                 faces, edges, edge_faces, acc = []) =
    let(next = next_face_around_vertex(v, f_cur, f_prev, faces, edges, edge_faces))
        (next == f_start)
            ? concat(acc, [f_cur])
            : faces_around_vertex_rec(
                  v,
                  next,
                  f_cur,
                  f_start,
                  faces, edges, edge_faces,
                  concat(acc, [f_cur])
          );


function faces_around_vertex(poly, v, edges, edge_faces) =
    let(
        faces = poly_faces(poly),
        inc   = vertex_incident_faces(poly, v),
        start = inc[0]
    )
    faces_around_vertex_rec(v, start, -1, start, faces, edges, edge_faces);


function dual_faces(poly, centers) =
    let(
        faces      = poly_faces(poly),
        edges      = _ps_edges_from_faces(faces),
        edge_faces = edge_faces_table(faces, edges),
        verts      = poly_verts(poly)
    )
    [
        for (vi = [0 : len(verts)-1])
            faces_around_vertex(poly, vi, edges, edge_faces)
    ];


function dual_unit_edge_and_e_over_ir(verts, faces) =
    let(
        edges    = _ps_edges_from_faces(faces),
        e0       = edges[0],
        vA       = verts[e0[0]],
        vB       = verts[e0[1]],
        unit_e   = norm(vB - vA),
        mid      = (vA + vB) / 2,
        ir       = norm(mid),
        e_over_ir = unit_e / ir
    )
    [unit_e, e_over_ir];


//function poly_verts_world(poly, edge_len) =
//    let(scale = edge_len / poly_unit_edge(poly))
//    poly_verts(poly) * scale;


// Polar dual vertices from face planes.
// verts: vertex positions (unit-edge OR world-space; doesn't matter as long as origin is inside)
// faces: outward oriented faces
function ps_face_polar_verts(verts, faces) =
    [
        for (fi = [0 : len(faces)-1])
            let(
                f = faces[fi],
                n = face_normal(verts, f),        // unit outward normal (given your face_normal)
                d = v_dot(n, verts[f[0]])         // plane offset along n
            )
            assert(d > 0, str("ps_face_polar_verts: d<=0 at face ", fi))
            (n / d)
    ];


//function poly_face_polar_verts(poly) =
//    ps_face_polar_verts(poly_verts(poly), poly_faces(poly));
//


// Build polar dual in the *same world units* as the given verts/faces.
// verts_world: already scaled positions (not unit-edge)
// faces: oriented outward
function poly_dual_polar_vf(verts, faces) =
    let(
        dual_verts  = ps_face_polar_verts(verts, faces),
        faces_raw   = dual_faces([verts, faces, 1, 1], dual_verts),
        faces_orient = orient_all_faces_outward(dual_verts, faces_raw)
    )
    [dual_verts, faces_orient];


//function poly_dual_polar_world(verts_world, faces) =
//    poly_dual_polar_vf(verts_world, faces);

//// Scale a poly descriptor's verts by k (faces unchanged; unit_edge/e_over_ir recomputed later)
//function _ps_poly_scale_verts(poly, k) =
//    [ poly_verts(poly) * k, poly_faces(poly), poly_unit_edge(poly), poly_e_over_ir(poly) ];



function ps_mean(a) =
    assert(len(a) > 0, "ps_mean: empty")
    sum(a) / len(a);

//function ps_vertex_radii(poly) =
//    [ for (v = poly_verts(poly)) norm(v) ];
//
//function ps_vertex_radius_stat(poly, mode) =
//    let(r = ps_vertex_radii(poly))
//    assert(len(r) > 0, "ps_vertex_radius_stat: no verts")
//    (mode == "vertex_max") ? max(r) : ps_mean(r);

//// “Apparent vertex radius per unit IR” under place_on_*_ir scaling
//function ps_apparent_radius_per_ir(poly, mode) =
//    let(
//        r_stat = ps_vertex_radius_stat(poly, mode),
//        s_ir   = poly_e_over_ir(poly) / poly_unit_edge(poly)
//    )
//    r_stat * s_ir;

//// Scale factor to apply to IR for 'dual' so it overlays nicely with 'poly'
//function ps_scale_dual_pair(poly, dual, mode="vertex_mean") =
//    let(
//        rp = ps_apparent_radius_per_ir(poly, mode),
//        rd = ps_apparent_radius_per_ir(dual, mode)
//    )
//    assert(rp > 0 && rd > 0, "ps_scale_dual_pair: nonpositive radii")
//    (rp / rd);

//// Return k such that: (your dual verts) q = k * (raw polar verts) p
//function ps_dual_polar_k(poly, dual) =
//    let(
//        verts = poly_verts(poly),
//        faces = orient_all_faces_outward(verts, poly_faces(poly)),
//        dv    = poly_verts(dual),
//
//        ks = [
//            for (fi = [0 : len(faces)-1])
//                let(
//                    f = faces[fi],
//                    n = face_normal(verts, f),          // unit outward normal
//                    d = v_dot(n, verts[f[0]]),          // plane offset
//                    q = dv[fi],                         // dual vertex corresponding to face fi
//                    _ = assert(d > 0, str("ps_dual_polar_k: d<=0 at face ", fi)),
//                    k_i = d * v_dot(q, n)
//                )
//                k_i
//        ]
//    )
//    ps_mean(ks);

//// Factor to multiply IR when placing 'dual' to get polar/reciprocal overlay with 'poly'
//function ps_scale_dual_polar(poly, dual) =
//    let(k = ps_dual_polar_k(poly, dual))
//    assert(k != 0, "ps_scale_dual_polar: k==0")
//    (1 / k);
//
//
//// argmax index for a numeric list (no sort required)
//function _ps_argmax(vals, i=0, best_i=0, best_v=-1e100) =
//    (i >= len(vals)) ? best_i :
//    (vals[i] > best_v)
//        ? _ps_argmax(vals, i+1, i,      vals[i])
//        : _ps_argmax(vals, i+1, best_i, best_v);

//// Find dual vertex that best matches face normal direction
//function _ps_match_dual_vert_idx(dverts, n_unit) =
//    let(
//        scores = [
//            for (j = [0 : len(dverts)-1])
//                let(q = dverts[j])
//                    v_dot(n_unit, v_norm(q))   // cosine similarity
//        ]
//    )
//    _ps_argmax(scores);

//// IR multiplier so that when overlaid, each face plane of 'poly' passes through
//// its corresponding dual vertex of 'dual' (direction-matched, order-independent).
//function ps_scale_dual_overlay_on_faces(poly, dual) =
//    let(
//        pv = poly_verts(poly),
//        pf = orient_all_faces_outward(pv, poly_faces(poly)),
//        dv = poly_verts(dual),
//
//        sp = poly_e_over_ir(poly) / poly_unit_edge(poly),
//        sd = poly_e_over_ir(dual) / poly_unit_edge(dual),
//
//        _0 = assert(len(dv) > 0, "ps_scale_dual_overlay_on_faces: dual has no verts"),
//        _1 = assert(len(pf) > 0, "ps_scale_dual_overlay_on_faces: poly has no faces"),
//
//        ms = [
//            for (fi = [0 : len(pf)-1])
//                let(
//                    f  = pf[fi],
//                    n  = face_normal(pv, f),          // unit outward normal
//                    d  = v_dot(n, pv[f[0]]),          // plane offset
//                    j  = _ps_match_dual_vert_idx(dv, n),
//                    q  = dv[j],
//                    nq = v_dot(n, q),
//                    m_i = (sp * d) / (sd * nq),
//                    _ = echo("face", fi, " -> dual vert", j, " d=", d, " n·q=", nq, " m_i=", m_i)
//                )
//                assert(d > 0,  str("ps_scale_dual_overlay_on_faces: d<=0 at face ", fi))
//                assert(abs(nq) > 1e-12, str("ps_scale_dual_overlay_on_faces: n·q ~ 0 at face ", fi))
//                m_i
//        ],
//        _3 = echo("overlay multiplier (mean)", sum(ms)/len(ms), "min", min(ms), "max", max(ms))
//    )
//        ps_mean(ms);

// ---- Edge midradius helpers ----

// Return list of edge midpoint radii for a poly in its own coordinate system (unit-edge verts)
function ps_edge_midradius_list(poly) =
    let(
        verts = poly_verts(poly),
        edges = _ps_edges_from_faces(poly_faces(poly)),
        rs = [
            for (e = edges)
                let(
                    m = (verts[e[0]] + verts[e[1]]) / 2
                )
                norm(m)
        ]
    )
    assert(len(rs) > 0, "ps_edge_midradius_list: poly has no edges")
    rs;

function ps_edge_midradius_stat(poly, mode="min") =
    let(rs = ps_edge_midradius_list(poly))
    (mode == "min") ? min(rs)
  : (mode == "max") ? max(rs)
  : ps_mean(rs);



// ---- IR overlay multiplier so dual's edges line up with poly's edges (mid-sphere style) ----
// Returns multiplier 'm' such that using IR*m for the dual tends to align edge crossings.
function scale_dual(poly, dual, mode="min") =
    let(
        // scaling from unit-edge coords to "per-IR world coords":
        // world = IR * (e_over_ir/unit_edge) * verts_unit
        sp = poly_e_over_ir(poly)  / poly_unit_edge(poly),
        sd = poly_e_over_ir(dual)  / poly_unit_edge(dual),

        rp = ps_edge_midradius_stat(poly,  mode),
        rd = ps_edge_midradius_stat(dual,  mode),

        _0 = assert(abs(rd) > 1e-12, "ps_scale_dual_overlay_on_edges: dual edge midradius ~ 0")
    )
    (sp * rp) / (sd * rd);


//function ps_scale_dual_polar_correct(poly, dual) =
//    let(
//        ir_poly = poly_unit_edge(poly) / poly_e_over_ir(poly),
//        ir_dual = poly_unit_edge(dual) / poly_e_over_ir(dual),
//        // For polar reciprocity: ir_poly * (ir_dual * scale) = 1
//        // So: scale = 1 / (ir_poly * ir_dual)
//        scale = 1.0 / (ir_poly * ir_dual)
//    )
//    scale;




// Public: polar dual, returned as a normalised poly descriptor.
// Default behaviour: returns a descriptor with unit_edge = 1 (library convention).
function poly_dual(poly) =
    let(
        // Ensure input faces are outward for correct polar normals
        verts0 = poly_verts(poly),
        faces0 = orient_all_faces_outward(verts0, poly_faces(poly)),

        // Build raw polar dual in same coordinate system as verts0
        dual_vf_raw = poly_dual_polar_vf(verts0, faces0),
        dv_raw = dual_vf_raw[0],
        df_raw = dual_vf_raw[1],

        // Compute metrics for raw dual
        ue_eir_raw = dual_unit_edge_and_e_over_ir(dv_raw, df_raw),
        unit_e_raw = ue_eir_raw[0],
        e_over_ir_raw = ue_eir_raw[1],

        // Renormalise so returned descriptor follows convention unit_edge = 1
        k = 1 / unit_e_raw,
        dv = dv_raw * k,

        // Recompute metrics after renormalisation
        ue_eir = dual_unit_edge_and_e_over_ir(dv, df_raw),
        unit_e = ue_eir[0],          // should be ~1
        e_over_ir = ue_eir[1]
    )
    make_poly(dv, df_raw, unit_e, e_over_ir);
