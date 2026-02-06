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
        ei1 = ps_find_edge_index(edges, v,      v_next),
        ei2 = ps_find_edge_index(edges, v_prev, v     ),

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
        edge_faces = ps_edge_faces_table(faces, edges),
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
                n = ps_face_normal(verts, f),        // unit outward normal
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
        faces_orient = ps_orient_all_faces_outward(dual_verts, faces_raw)
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

function ps_edge_midradius_stat(poly) =
    let(rs = ps_edge_midradius_list(poly))
        min(rs);

// ---- Face radius helpers ----

// Mean vertex distance from face centroid for each face (unit-edge coords).
function ps_face_radius_list(poly) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly)
    )
    [
        for (f = faces)
            let(
                c = ps_face_centroid(verts, f),
                rs = [ for (vid = f) norm(verts[vid] - c) ]
            )
            sum(rs) / len(rs)
    ];

function ps_face_radius_stat(poly, face_k=undef) =
    let(
        faces = poly_faces(poly),
        rs_all = ps_face_radius_list(poly),
        rs = is_undef(face_k)
            ? rs_all
            : [ for (i = [0:len(faces)-1]) if (len(faces[i]) == face_k) rs_all[i] ],
        _0 = assert(len(rs) > 0, "ps_face_radius_stat: no faces of that size")
    )
    min(rs);

// ---- Face-family helpers ----

// Return [[k, count], ...] for each face size k, sorted by k ascending.
function ps_face_family_list(poly) =
    let(
        faces = poly_faces(poly),
        sizes = [ for (f = faces) len(f) ],
        uniq = [
            for (i = [0:len(sizes)-1])
                let(k = sizes[i])
                    if (sum([ for (j = [0:1:i-1]) sizes[j] == k ? 1 : 0 ]) == 0) k
        ],
        ks = sort(uniq)
    )
    [ for (k = ks) [k, sum([ for (s = sizes) s == k ? 1 : 0 ])] ];

// Return [k, count] for the face size that appears most frequently.
// Ties are resolved by choosing the smallest k.
function ps_face_family_mode(poly) =
    let(
        faces = poly_faces(poly),
        sizes = [ for (f = faces) len(f) ],
        uniq = [
            for (i = [0:len(sizes)-1])
                let(k = sizes[i])
                    if (sum([ for (j = [0:1:i-1]) sizes[j] == k ? 1 : 0 ]) == 0) k
        ],
        counts = [ for (k = uniq) sum([ for (s = sizes) s == k ? 1 : 0 ]) ],
        max_count = max(counts),
        best = [
            for (i = [0:len(uniq)-1])
                if (counts[i] == max_count) uniq[i]
        ],
        k = min(best)
    )
    [k, max_count];

// Return [k, count] for the largest face size.
function ps_face_family_max(poly) =
    let(
        faces = poly_faces(poly),
        sizes = [ for (f = faces) len(f) ],
        k = max(sizes),
        count = sum([ for (s = sizes) s == k ? 1 : 0 ])
    )
    [k, count];

// ---- Edge-crossing scale helpers ----

function _ps_vertex_valence_list(verts, edges) =
    [
        for (vi = [0:len(verts)-1])
            sum([
                for (e = edges)
                    (e[0] == vi || e[1] == vi) ? 1 : 0
            ])
    ];

function _ps_edge_signature_full(edges, faces, edge_faces, valences, ei) =
    let(
        fpair = edge_faces[ei],
        k0 = len(faces[fpair[0]]),
        k1 = len(faces[fpair[1]]),
        ks = (k0 < k1) ? [k0, k1] : [k1, k0],
        e = edges[ei],
        v0 = valences[e[0]],
        v1 = valences[e[1]],
        vs = (v0 < v1) ? [v0, v1] : [v1, v0]
    )
    [ks[0], ks[1], vs[0], vs[1]];

// Return edge index for an edge on face face_idx at edge_pos.
function ps_edge_from_face(poly, face_idx, edge_pos) =
    let(
        faces = poly_faces(poly),
        f = faces[face_idx],
        n = len(f),
        a = f[edge_pos % n],
        b = f[(edge_pos + 1) % n],
        edges = _ps_edges_from_faces(faces)
    )
    ps_find_edge_index(edges, a, b);

// Compute a scale that makes dual edges intersect the selected edge family.
// face_idx selects a face; edge_pos selects the reference edge within that face.
// The edge family is defined by matching adjacent face sizes and endpoint valences.
function scale_dual_edge_cross(poly, dual, face_idx, edge_pos=0, eps=1e-12, len_eps=1e-6) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        edges = _ps_edges_from_faces(faces),
        edge_faces = ps_edge_faces_table(faces, edges),
        valences = _ps_vertex_valence_list(verts, edges),

        ref_ei = ps_edge_from_face(poly, face_idx, edge_pos),
        ref_sig = _ps_edge_signature_full(edges, faces, edge_faces, valences, ref_ei),
        ref_len = norm(verts[edges[ref_ei][1]] - verts[edges[ref_ei][0]]),

        // gather edges matching the signature
        edge_ids = [
            for (ei = [0:len(edges)-1])
                let(
                    sig = _ps_edge_signature_full(edges, faces, edge_faces, valences, ei),
                    len_e = norm(verts[edges[ei][1]] - verts[edges[ei][0]])
                )
                if (sig == ref_sig && abs(len_e - ref_len) <= len_eps) ei
        ],

        dverts = poly_verts(dual),
        sp = poly_e_over_ir(poly),
        sd = poly_e_over_ir(dual),

        scales = [
            for (ei = edge_ids)
                let(
                    e = edges[ei],
                    A = verts[e[0]],
                    B = verts[e[1]],
                    d_f = edge_faces[ei],
                    D0 = dverts[d_f[0]],
                    D1 = dverts[d_f[1]],
                    E = D1 - D0,
                    M = [
                        [ (B-A)[0], -D0[0], -E[0] ],
                        [ (B-A)[1], -D0[1], -E[1] ],
                        [ (B-A)[2], -D0[2], -E[2] ]
                    ],
                    sol = _ps_solve3(M, -A, eps)
                )
                is_undef(sol) ? undef
              : let(
                    v = sol[0],
                    s = sol[1],
                    y = sol[2],
                    u = (abs(s) < eps) ? undef : y / s,
                    ok = (s > 0) && (v >= 0) && (v <= 1) && (!is_undef(u)) && (u >= 0) && (u <= 1)
                )
                ok ? (s * sp / sd) : undef
        ],

        s_vals = [ for (s = scales) if (!is_undef(s)) s ],
        _0 = assert(len(s_vals) > 0, "scale_dual_edge_cross: no valid edges found"),
        sorted = _ps_sort(s_vals),
        n = len(sorted)
    )
    (n % 2 == 1)
        ? sorted[(n-1)/2]
        : (sorted[n/2 - 1] + sorted[n/2]) / 2;


// ---- IR overlay multiplier so dual's edges line up with poly's edges (mid-sphere style) ----
// Returns multiplier 'm' such that using IR*m for the dual tends to align edge crossings.
function scale_dual(poly, dual) =
    let(
        // scaling from unit-edge coords to "per-IR world coords":
        // world = IR * (e_over_ir/unit_edge) * verts_unit
        sp = poly_e_over_ir(poly),
        sd = poly_e_over_ir(dual),

        rp = ps_edge_midradius_stat(poly),
        rd = ps_edge_midradius_stat(dual),

        _0 = assert(abs(rd) > 1e-12, "scale_dual: dual edge midradius ~ 0")
    )
    (sp * rp) / (sd * rd);

// ---- Face-radius scaling to align face overlays ----
// Returns multiplier 'm' so that the unit-edge face radius of poly matches the dual's.
// Use face_k (vertex count) to select a face family on the poly (e.g., 4 for squares).
// For world-space overlays, multiply this by scale_dual(poly, dual).
function scale_dual_face_radius(poly, dual, face_k=undef, dual_face_k=undef) =
    let(
        rp = ps_face_radius_stat(poly, face_k),
        rd = ps_face_radius_stat(dual, dual_face_k),
        _0 = assert(abs(rd) > 1e-12, "scale_dual_face_radius: dual face radius ~ 0")
    )
    rp / rd * scale_dual(poly, dual);



// Public: polar dual, returned as a normalised poly descriptor.
// Default behaviour: returns a descriptor with unit_edge = 1 (library convention).
function poly_dual(poly) =
    let(
        // Ensure input faces are outward for correct polar normals
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),

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
