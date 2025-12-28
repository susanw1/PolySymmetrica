// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <placement.scad>


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


// Polar dual vertices from face planes.
// verts: vertex positions (unit-edge OR world-space; doesn't matter as long as origin is inside)
// faces: outward oriented faces
function ps_face_polar_verts(verts, faces) =
    [
        for (fi = [0 : len(faces)-1])
            let(
                f = faces[fi],
                n = face_normal(verts, f),        // unit outward normal
                d = v_dot(n, verts[f[0]])         // plane offset along n
            )
            assert(d > 0, str("ps_face_polar_verts: d<=0 at face ", fi))
            (n / d)
    ];


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
  : sum(rs) / len(rs);



// ---- IR overlay multiplier so dual's edges line up with poly's edges (mid-sphere style) ----
// Returns multiplier 'm' such that using IR*m for the dual tends to align edge crossings.
function scale_dual(poly, dual, mode="min") =
    let(
        // scaling from unit-edge coords to "per-IR world coords":
        // world = IR * (e_over_ir/unit_edge) * verts_unit
        sp = poly_e_over_ir(poly),
        sd = poly_e_over_ir(dual),

        rp = ps_edge_midradius_stat(poly,  mode),
        rd = ps_edge_midradius_stat(dual,  mode),

        _0 = assert(abs(rd) > 1e-12, "scale_dual: dual edge midradius ~ 0")
    )
    (sp * rp) / (sd * rd);



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
        e_over_ir = ue_eir[1]
    )
    make_poly(dv, df_raw, e_over_ir);
