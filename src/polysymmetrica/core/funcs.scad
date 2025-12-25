// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT


///////////////////////////////////////
// ---- Poly descriptor accessors ----
function poly_verts(poly)      = poly[0];
function poly_faces(poly)      = poly[1];
function poly_e_over_ir(poly)  = poly[2];

// Constructor with full validation and optional auto-computation
function make_poly(verts, faces, e_over_ir=undef) =
    let(
        // Validation
        _0 = assert(len(verts) >= 4, "Polyhedron must have at least 4 vertices"),
        _1 = assert(len(faces) >= 4, "Polyhedron must have at least 4 faces"),
        _2 = assert(all_faces_valid(verts, faces), "Invalid face indices"),

        // Auto-compute if not provided
        edges = _ps_edges_from_faces(faces),
        _3 = assert(len(edges) >= 6, "Polyhedron must have at least 6 edges"),

        // compute ir from min edge-midradius, not just the first edge
        mids = [
            for (e = edges)
                norm((verts[e[0]] + verts[e[1]]) / 2)
        ],
        ir = min(mids),
        _ir_ok = assert(ir > 0, "make_poly: inter-radius (min edge-midradius) must be positive"),

        // choose an edge achieving that min (first one that matches)
        ei_ir = [ for (i = [0:len(edges)-1]) if (abs(mids[i] - ir) < 1e-12) i ][0],
        e_ir  = edges[ei_ir],

        computed_e_over_ir = is_undef(e_over_ir)
            ? norm(verts[e_ir[1]] - verts[e_ir[0]]) / ir
            : e_over_ir,
        
        _5 = assert(computed_e_over_ir > 0, "e_over_ir must be positive")
    )
    [verts, faces, computed_e_over_ir];

// Helper validation
function all_faces_valid(verts, faces) =
    len([
        for (f = faces)
            if (len(f) >= 3 && all_indices_in_range(f, len(verts)))
                1
    ]) == len(faces);

function all_indices_in_range(face, max_idx) =
    len([for (vi = face) if (vi >= 0 && vi < max_idx) 1]) == len(face);



///////////////////////////////////////
// Handy vector functions (and aliases)

function v_add(a, b)   = a + b;
function v_sub(a, b)   = a - b;
function v_scale(a, k) = a * k;             // scalar multiplication
function v_dot(a, b)   = a * b;             // dot product
function v_cross(a, b) = cross(a, b);       // OpenSCAD built-in
function v_len(a)      = norm(a);           // built-in length
function v_norm(a)     = let(L = norm(a)) (L == 0 ? [0,0,0] : a / L);

// Edge equality
function edge_equal(e1, e2) = (e1[0] == e2[0] && e1[1] == e2[1]);

/** Sum of vector */
function sum(a, i = 0) =
    i >= len(a) ? 0 : a[i] + sum(a, i + 1);

/** Sum of vector list (3D) */
function v_sum(list) = [
    sum([for (v = list) v[0]]),
    sum([for (v = list) v[1]]),
    sum([for (v = list) v[2]])
];

///////////////////////////////////////
// Polygon helpers

/** Calculate polygon edge, given N and radius */
function calc_edge(n_vertex, rad) = 2 * rad * sin(180 / n_vertex);

/** Calculate polygon radius, given N and edge length */
function calc_radius(n_vertex, edge_len) = edge_len / (2 * sin(180 / n_vertex));


///////////////////////////////////////
// Geometry helpers

/** Face centroid from verts + index list */
function face_centroid(verts, f) =
    len(f) == 0
        ? [0,0,0]
        : v_scale(v_sum([for (vid = f) verts[vid]]), 1 / len(f));

// Face normal (not scaled, just direction)
function face_normal(verts, f) =
    v_norm(v_cross(
        verts[f[1]] - verts[f[0]],
        verts[f[2]] - verts[f[0]]
    ));


// Return index of some neighbour vertex of vi (from the face list)
function poly_vertex_neighbor(poly, vi) =
    let(
        faces = poly_faces(poly),
        // collect "next" vertex after vi in any face that contains it
        candidates = [
            for (f = faces)
                for (k = [0 : len(f)-1])
                    if (f[k] == vi) f[(k+1) % len(f)]
        ]
    ) candidates[0];  // first one is enough


// Internal: build unique undirected edges from a face list
function _ps_edges_from_faces(faces) =
    let(
        raw_edges = [
            for (fi = [0 : len(faces)-1])
                let(f = faces[fi])
                    for (k = [0 : len(f)-1])
                        let(
                            a = f[k],
                            b = f[(k+1) % len(f)],
                            e = (a < b) ? [a,b] : [b,a]
                        ) e
        ],
        uniq_edges = [
            for (i = [0 : len(raw_edges)-1])
                let(ei = raw_edges[i])
                    if (sum([
                            for (j = [0 : 1 : i-1])
                                edge_equal(raw_edges[j], ei) ? 1 : 0
                        ]) == 0) ei
        ]
    )
    uniq_edges;


// For each edge, list the indices of the faces incident to it (should be 2)
function edge_faces_table(faces, edges) =
    [
        for (ei = [0 : len(edges)-1])
            let(e = edges[ei])
            [
                for (fi = [0 : len(faces)-1])
                    if (face_has_edge(faces[fi], e[0], e[1])) fi
            ]
    ];



// Does face f contain undirected edge {a,b}?
function face_has_edge(f, a, b) =
    sum([
        for (k = [0 : len(f)-1])
            let(
                x = f[k],
                y = f[(k+1) % len(f)]
            )
            ((x==a && y==b) || (x==b && y==a)) ? 1 : 0
    ]) > 0;


// ---- Generic face-frame helpers ----
function poly_face_center(poly, fi, scale) =
    let(
        f   = poly_faces(poly)[fi],
        vs  = poly_verts(poly),
        xs  = [ for (vid = f) vs[vid][0] * scale ],
        ys  = [ for (vid = f) vs[vid][1] * scale ],
        zs  = [ for (vid = f) vs[vid][2] * scale ]
    )
    [
        sum(xs) / len(f),
        sum(ys) / len(f),
        sum(zs) / len(f)
    ];


function poly_face_ex(poly, fi, scale) =
    let(f      = poly_faces(poly)[fi],
        vs     = poly_verts(poly),
        center = poly_face_center(poly, fi, scale),
        v0     = vs[f[0]] * scale)
    v_norm(v0 - center);   // local +X points to vertex 0


function poly_face_ey(poly, fi, scale) =
    v_cross(
        poly_face_ez(poly, fi, scale),
        poly_face_ex(poly, fi, scale)
    );


function poly_face_ez(poly, fi, scale) =
    let(f  = poly_faces(poly)[fi],
        vs = poly_verts(poly),
        v0 = vs[f[0]] * scale,
        v1 = vs[f[1]] * scale,
        v2 = vs[f[2]] * scale)
    v_norm(v_cross(v1 - v0, v2 - v0));   // outward normal


function frame_matrix(center, ex, ey, ez) = [
    [ex[0], ey[0], ez[0], center[0]],
    [ex[1], ey[1], ez[1], center[1]],
    [ex[2], ey[2], ez[2], center[2]],
    [0,      0,     0,     1]
];


// Ensure face orientation so normal points outward (centroid·normal > 0)
function orient_face_outward(verts, f) =
    let(
        c = face_centroid(verts, f),
        n = face_normal(verts, f)
    )
    (v_dot(c, n) >= 0)
        ? f
        : [ for (i = [len(f)-1 : -1 : 0]) f[i] ];  // reversed

function orient_all_faces_outward(verts, faces) =
    [ for (f = faces) orient_face_outward(verts, f) ];

