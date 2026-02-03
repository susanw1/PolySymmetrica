// ---------------------------------------------------------------------------
// PolySymmetrica - Cantellation / Expansion
// Shared cantellation helpers and operators.

use <funcs.scad>
use <duals.scad>  // faces_around_vertex
use <transform.scad>
use <transform_util.scad>

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

function _ps_cantellate_df_from_c(poly, c, df_max=undef, steps=16, family_edge_idx=0) =
    let(
        verts = poly_verts(poly),
        edges = _ps_edges_from_faces(poly_faces(poly)),
        e0 = edges[0],
        edge_len = norm(verts[e0[1]] - verts[e0[0]]),
        ir = edge_len / poly_e_over_ir(poly),
        df_max_eff = is_undef(df_max) ? (2 * ir) : df_max,
        df_mid = cantellate_square_df(poly, 0, df_max_eff, steps, family_edge_idx),
        df = (c <= 0.5)
            ? (2 * c * df_mid)
            : (df_mid + (c - 0.5) * 2 * (df_max_eff - df_mid))
    )
    df;

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
