// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

// Shared-helper quick index. Prefer reusing these before adding new local copies:
// - Poly descriptors: poly_verts, poly_faces, poly_e_over_ir, poly_edges, make_poly
// - Basic validation: ps_faces_valid, ps_indices_in_range
// - Vector math: v_*, ps_rot_axis
// - Face/edge topology: ps_find_edge_index, ps_edge_faces_table, ps_face_has_edge, faces_around_vertex
// - Face geometry/frames: ps_face_centroid, ps_face_normal, ps_face_frame_normal, poly_face_center/ex/ey/ez
// - Orientation: poly_fix_winding, ps_orient_face_outward, ps_orient_all_faces_outward

///////////////////////////////////////
// ---- Poly descriptor API ----
/**
 * Function: Return the vertex list from a poly descriptor.
 * Params: poly ([verts, faces, e_over_ir] descriptor).
 * Returns: `poly[0]`.
 */
function poly_verts(poly)      = poly[0];
/**
 * Function: Return the face list from a poly descriptor.
 * Params: poly ([verts, faces, e_over_ir] descriptor).
 * Returns: `poly[1]`.
 */
function poly_faces(poly)      = poly[1];
/**
 * Function: Return the stored edge-to-inter-radius scale from a poly descriptor.
 * Params: poly ([verts, faces, e_over_ir] descriptor).
 * Returns: `poly[2]`.
 */
function poly_e_over_ir(poly)  = poly[2];
/**
 * Function: Derive the unique undirected edge list for a poly descriptor.
 * Params: poly ([verts, faces, e_over_ir] descriptor).
 * Returns: Canonical undirected edges `[[a, b] ...]`.
 * Limitations: Derived from faces each time; prefer `_ps_edges_from_faces(faces)` when faces are already available.
 */
function poly_edges(poly)      = _ps_edges_from_faces(poly_faces(poly));

/**
 * Function: Clamp a scalar to an inclusive range.
 * Params: x (value), lo (lower bound), hi (upper bound).
 * Returns: `x` clamped to `[lo, hi]`.
 */
function ps_clamp(x, lo, hi) = min(max(x, lo), hi);

/**
 * Function: Construct a poly descriptor and auto-compute `e_over_ir` when omitted.
 * Params: verts (`[[x,y,z] ...]`), faces (`[[i0,i1,...] ...]`), e_over_ir (optional override).
 * Returns: Normalized poly descriptor `[verts_centered, faces, e_over_ir]`.
 * Limitations: Structural validation only; does not itself guarantee manifoldness or closure.
 */
function make_poly(verts, faces, e_over_ir=undef) =
    let(
        // Validation
        _0 = assert(len(verts) >= 3, "Polyhedron must have at least 3 vertices"),
        _1 = assert(len(faces) >= 1, "Polyhedron must have at least 1 face"),
        _2 = assert(ps_faces_valid(verts, faces), "Invalid face indices"),

        // Auto-compute if not provided
        edges = _ps_edges_from_faces(faces),
        _3 = assert(len(edges) >= 1, "Polyhedron must have at least 1 edge"),
        center = _ps_poly_mid_center(verts, faces),
        verts_centered = [for (v = verts) v - center],

        // compute ir from min edge-midradius, not just the first edge
        mids = [
            for (e = edges)
                norm((verts_centered[e[0]] + verts_centered[e[1]]) / 2)
        ],
        ir = min(mids),
        _ir_ok = assert(ir > 0, "make_poly: inter-radius (min edge-midradius) must be positive"),

        // choose an edge achieving that min (first one that matches)
        ei_ir = [ for (i = [0:len(edges)-1]) if (abs(mids[i] - ir) < 1e-12) i ][0],
        e_ir  = edges[ei_ir],

        computed_e_over_ir = is_undef(e_over_ir)
            ? norm(verts_centered[e_ir[1]] - verts_centered[e_ir[0]]) / ir
            : e_over_ir,
        
        _5 = assert(computed_e_over_ir > 0, "e_over_ir must be positive")
    )
    [verts_centered, faces, computed_e_over_ir];

/**
 * Function: Reorient faces so shared undirected edges appear with opposite directions.
 * Params: poly (poly descriptor to repair).
 * Returns: Poly descriptor with corrected face winding.
 * Limitations: Topological winding repair only; does not decide outwardness from geometry.
 */
function poly_fix_winding(poly) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        fixed = _ps_fix_winding_all(faces)
    )
    [verts, fixed, poly_e_over_ir(poly)];

///////////////////////////////////////
// ---- Basic validation helpers ----
/**
 * Function: Check that every face has arity >= 3 and all indices are in range.
 * Params: verts (vertex list), faces (face index loops).
 * Returns: `true` when every face satisfies the basic structural checks.
 */
function ps_faces_valid(verts, faces) =
    len([
        for (f = faces)
            if (len(f) >= 3 && ps_indices_in_range(f, len(verts)))
                1
    ]) == len(faces);

/**
 * Function: Check that all indices in one face lie within `[0, max_idx)`.
 * Params: face (face index loop), max_idx (exclusive upper bound).
 * Returns: `true` when all indices are valid.
 */
function ps_indices_in_range(face, max_idx) =
    len([for (vi = face) if (vi >= 0 && vi < max_idx) 1]) == len(face);

///////////////////////////////////////
// ---- Winding/orientation helpers (private) ----
function _ps_face_edges_dir(f) =
    let(n = len(f))
    [ for (i = [0:1:n-1]) [f[i], f[(i+1)%n]] ];

function _ps_face_edge_dir(f, a, b) =
    let(
        n = len(f),
        vals = [
            for (i = [0:1:n-1])
                let(u = f[i], v = f[(i+1)%n])
                (u == a && v == b) ? 1 : (u == b && v == a) ? -1 : 0
        ]
    )
    (max(vals) == 1) ? 1 : (min(vals) == -1) ? -1 : 0;

function _ps_adjacent_faces_for_edge(faces, a, b, fi) =
    [
        for (fj = [0:1:len(faces)-1])
            if (fj != fi && _ps_face_edge_dir(faces[fj], a, b) != 0) fj
    ];

function _ps_list_set(list, idx, val) =
    [ for (i = [0:1:len(list)-1]) (i == idx) ? val : list[i] ];

function _ps_index_of_undef(list) =
    let(idx = [for (i = [0:1:len(list)-1]) if (is_undef(list[i])) i])
    (len(idx) == 0) ? -1 : idx[0];

function _ps_reverse(list) =
    [ for (i = [len(list)-1 : -1 : 0]) list[i] ];

function _ps_identity_map(n) =
    [for (i = [0:1:n-1]) i];

function _ps_distinct_count(list) =
    len([for (i = [0:1:len(list)-1]) if (_ps_index_of(list, list[i]) == i) 1]);

function _ps_face_strip_adjacent_dups(f) =
    let(n = len(f))
    (n == 0) ? [] :
    [for (i = [0:1:n-1]) if (f[i] != f[(i-1+n)%n]) f[i]];

function _ps_face_trim_closing_dup(f) =
    (len(f) >= 2 && f[0] == f[len(f)-1]) ? [for (i = [0:1:len(f)-2]) f[i]] : f;

function _ps_face_clean_cycle(f) =
    _ps_face_trim_closing_dup(_ps_face_strip_adjacent_dups(f));

function _ps_faces_clean_cycles(faces) =
    [for (f = faces) _ps_face_clean_cycle(f)];

function _ps_faces_remap(faces, old_to_new) =
    [for (f = faces) [for (vi = f) old_to_new[vi]]];

function _ps_neighbor_face(neighbors, idx) =
    let(hit = [for (p = neighbors) if (p[0] == idx) p[1]])
    (len(hit) == 0) ? undef : hit[0];

function _ps_apply_neighbors(fixed, neighbors) =
    [
        for (i = [0:1:len(fixed)-1])
            is_undef(fixed[i]) ? _ps_neighbor_face(neighbors, i) : fixed[i]
    ];

function _ps_new_neighbor_indices(fixed, neighbors) =
    [ for (p = neighbors) if (is_undef(fixed[p[0]])) p[0] ];

function _ps_fix_winding_queue(faces, fixed, queue) =
    (len(queue) == 0) ? fixed :
    let(
        fi = queue[0],
        fcur = fixed[fi],
        edges = _ps_face_edges_dir(fcur),
        neighbors = [
            for (e = edges)
                let(
                    a = e[0],
                    b = e[1],
                    nbrs = _ps_adjacent_faces_for_edge(faces, a, b, fi)
                )
                for (nj = nbrs)
                    let(
                        dir_cur = _ps_face_edge_dir(fcur, a, b),
                        dir_n = _ps_face_edge_dir(faces[nj], a, b),
                        desired = (dir_n == dir_cur) ? _ps_reverse(faces[nj]) : faces[nj]
                    )
                    [nj, desired]
        ],
        updated = _ps_apply_neighbors(fixed, neighbors),
        new_queue = concat([for (i = [1:1:len(queue)-1]) queue[i]], _ps_new_neighbor_indices(fixed, neighbors))
    )
    _ps_fix_winding_queue(faces, updated, new_queue);

function _ps_fix_winding_all(faces, fixed=undef) =
    let(
        init = is_undef(fixed) ? [for (i = [0:1:len(faces)-1]) undef] : fixed,
        seed = _ps_index_of_undef(init)
    )
    (seed < 0) ? init :
    _ps_fix_winding_all(
        faces,
        _ps_fix_winding_queue(faces, _ps_list_set(init, seed, faces[seed]), [seed])
    );

///////////////////////////////////////
// ---- List helpers (private) ----
function _ps_list_contains(list, v) =
    len(search(v, list)) > 0;

function _ps_index_of(list, v) =
    let(idx = search(v, list))
    (len(idx) == 0) ? -1 : idx[0];

///////////////////////////////////////
// ---- Vector math ----
/**
 * Function: Add two vectors component-wise.
 * Params: a, b (same-dimension vectors).
 * Returns: `a + b`.
 */
function v_add(a, b)   = a + b;
/**
 * Function: Subtract one vector from another component-wise.
 * Params: a, b (same-dimension vectors).
 * Returns: `a - b`.
 */
function v_sub(a, b)   = a - b;
/**
 * Function: Scale a vector by a scalar.
 * Params: a (vector), k (scalar multiplier).
 * Returns: `a * k`.
 */
function v_scale(a, k) = a * k;             // scalar multiplication
/**
 * Function: Compute the dot product of two vectors.
 * Params: a, b (same-dimension vectors).
 * Returns: Scalar dot product.
 */
function v_dot(a, b)   = a * b;             // dot product
/**
 * Function: Compute the 3D cross product of two vectors.
 * Params: a, b (3D vectors).
 * Returns: `cross(a, b)`.
 */
function v_cross(a, b) = cross(a, b);       // OpenSCAD built-in
/**
 * Function: Compute the Euclidean norm of a vector.
 * Params: a (vector).
 * Returns: `norm(a)`.
 */
function v_len(a)      = norm(a);           // built-in length
/**
 * Function: Normalize a vector safely.
 * Params: a (vector).
 * Returns: Unit vector, or `[0,0,0]` when `a` is zero.
 * Limitations: Zero-vector fallback is silent; assert first if that is not acceptable.
 */
function v_norm(a)     = let(L = norm(a)) (L == 0 ? [0,0,0] : a / L);

function _ps_ordered_pair(a, b) = (a < b) ? [a,b] : [b,a];


///////////////////////////////////////
// ---- Edge/list primitives ----
/**
 * Function: Compare two canonical undirected edge pairs for exact equality.
 * Params: e1, e2 (edge pairs `[a, b]`).
 * Returns: `true` when the pairs match exactly.
 * Limitations: Does not normalize orientation; callers should canonicalize first when needed.
 */
function edge_equal(e1, e2) = (e1[0] == e2[0] && e1[1] == e2[1]);

/**
 * Function: Sum a list of scalars recursively.
 * Params: a (scalar list), i (internal recursion index).
 * Returns: Sum of all entries in `a`.
 */
function sum(a, i = 0) =
    i >= len(a) ? 0 : a[i] + sum(a, i + 1);

/**
 * Function: Sum a list of vectors component-wise.
 * Params: list (vector list).
 * Returns: One vector with the component-wise sum.
 * Limitations: Assumes all vectors in the list have the same dimension.
 */
function v_sum(list) =
    (len(list) == 0) ? [] :
    let(n = len(list[0]))
    [ for (i = [0:1:n-1]) sum([for (v = list) v[i]]) ];

/**
 * Function: Rotate a 3D vector about an axis by an angle in degrees.
 * Params: v (vector), axis (rotation axis), ang (degrees).
 * Returns: Rotated vector.
 * Limitations: Zero axis collapses to the `v_norm(...)` zero-vector behavior.
 */
function ps_rot_axis(v, axis, ang) =
    let(
        a = v_norm(axis),
        c = cos(ang),
        s = sin(ang),
        term1 = v_scale(v, c),
        term2 = v_scale(v_cross(a, v), s),
        term3 = v_scale(a, v_dot(a, v) * (1 - c))
    )
    v_add(v_add(term1, term2), term3);

/**
 * Function: Find the index of an undirected edge in a canonical edge list.
 * Params: edges (canonical undirected edge list), a,b (edge endpoints).
 * Returns: Index of `{a, b}` within `edges`.
 * Limitations: Assumes the edge exists; intentionally does not return `undef` on miss.
 */
function ps_find_edge_index(edges, a, b) =
    let(
        e = _ps_ordered_pair(a, b),
        idxs = [for (i = [0 : len(edges)-1]) if (edge_equal(edges[i], e)) i]
    )
    idxs[0];   // assume the edge exists

/**
 * Function: Compare two points within an absolute epsilon.
 * Params: p,q (points/vectors), eps (equality tolerance).
 * Returns: `true` when `norm(p - q) <= eps`.
 */
function ps_point_eq(p,q,eps) = norm(p-q) <= eps;

function _ps_list_min(list, i=0, cur=undef) =
    (i >= len(list)) ? cur :
    let(v = list[i])
    _ps_list_min(list, i+1, is_undef(cur) ? v : (v < cur ? v : cur));

function _ps_remove_first(list, v, i=0) =
    (i >= len(list)) ? [] :
    (list[i] == v) ? [for (j = [i+1:1:len(list)-1]) list[j]]
                  : concat([list[i]], _ps_remove_first(list, v, i+1));

function _ps_sort(list, acc=[]) =
    (len(list) == 0) ? acc :
    let(mn = _ps_list_min(list))
    _ps_sort(_ps_remove_first(list, mn), concat(acc, [mn]));

///////////////////////////////////////
// ---- Polygon/polygram helpers ----
// Euclidean gcd for integer validation.
function _ps_gcd(a, b) =
    let(ai = abs(round(a)), bi = abs(round(b)))
    (bi == 0) ? ai : _ps_gcd(bi, ai % bi);

// Validate Schläfli-like polygon params {n,p} for single-cycle polygrams.
function _ps_validate_np(n, p, who) =
    let(
        n_i = round(n),
        p_i = round(p),
        _n_int = assert(abs(n - n_i) < 1e-9, str(who, ": n must be an integer")),
        _p_int = assert(abs(p - p_i) < 1e-9, str(who, ": p must be an integer")),
        _n_ok = assert(n_i >= 3, str(who, ": n must be >= 3")),
        _p_ok = assert(p_i >= 1 && (2 * p_i) < n_i, str(who, ": p must satisfy 1 <= p < n/2")),
        _cop = assert(_ps_gcd(n_i, p_i) == 1, str(who, ": n and p must be coprime"))
    )
    [n_i, p_i];

// Circumradius for regular/star polygon {n,p} with edge/chord length `edge`.
function _ps_polygram_radius(n, p, edge) =
    edge / (2 * sin(180 * p / n));

// Backward-compatible regular n-gon helper ({n,1}).
function _ps_ngon_radius(n, edge) =
    _ps_polygram_radius(n, 1, edge);

// Single-cycle vertex order for {n,p}.
function _ps_polygram_cycle(n, p) =
    [for (k = [0:1:n-1]) (k * p) % n];

// n-gon/polygram support ring at fixed z.
function _ps_ngon_ring(n, radius, z, phase=0) =
    [for (k = [0:1:n-1]) [radius * cos(360 * k / n + phase), radius * sin(360 * k / n + phase), z]];

// Inter-radius from the minimum edge-midpoint radius of an explicit verts/faces mesh.
function _ps_poly_mid_center(verts, faces) =
    let(
        edges = _ps_edges_from_faces(faces),
        mids = [for (e = edges) (verts[e[0]] + verts[e[1]]) / 2],
        _ok = assert(len(mids) > 0, "poly: edge-midpoint center requires at least one edge")
    )
    v_scale(v_sum(mids), 1 / len(mids));

function _ps_poly_ir(verts, faces) =
    let(
        edges = _ps_edges_from_faces(faces),
        center = _ps_poly_mid_center(verts, faces),
        verts_centered = [for (v = verts) v - center],
        mids = [for (e = edges) norm((verts_centered[e[0]] + verts_centered[e[1]]) / 2)],
        ir = min(mids),
        _ok = assert(ir > 0, "poly: inter-radius must be > 0")
    ) ir;

///////////////////////////////////////
// ---- Linear algebra helpers ----
function _ps_det3(m) =
    m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
    m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
    m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

function _ps_replace_col(m, col, b) =
    [
        [ col == 0 ? b[0] : m[0][0], col == 1 ? b[0] : m[0][1], col == 2 ? b[0] : m[0][2] ],
        [ col == 0 ? b[1] : m[1][0], col == 1 ? b[1] : m[1][1], col == 2 ? b[1] : m[1][2] ],
        [ col == 0 ? b[2] : m[2][0], col == 1 ? b[2] : m[2][1], col == 2 ? b[2] : m[2][2] ]
    ];

// Solve 3x3 linear system M*x = b via Cramer's rule; returns [x0,x1,x2] or undef.
function _ps_solve3(m, b, eps=1e-12) =
    let(det = _ps_det3(m))
    (abs(det) < eps) ? undef
  : [
        _ps_det3(_ps_replace_col(m, 0, b)) / det,
        _ps_det3(_ps_replace_col(m, 1, b)) / det,
        _ps_det3(_ps_replace_col(m, 2, b)) / det
    ];

///////////////////////////////////////
// ---- Prefix offsets / face offsets ----
// Prefix offsets for a list of counts: [0, c0, c0+c1, ...]
function _ps_prefix_offsets(counts, acc=[]) =
    (len(counts) == 0) ? acc :
    let(last = (len(acc) == 0) ? 0 : acc[len(acc)-1])
    _ps_prefix_offsets(
        [for (i = [1:1:len(counts)-1]) counts[i]],
        concat(acc, [last + counts[0]])
    );

function _ps_face_offsets(faces) =
    let(counts = [for (f = faces) len(f)])
    _ps_prefix_offsets(counts, [0]);

function _ps_face_edge_offsets(faces) =
    let(counts = [for (f = faces) 2 * len(f)])
    _ps_prefix_offsets(counts, [0]);


///////////////////////////////////////
// ---- Polygon helpers ----
/**
 * Function: Compute the edge length of a regular `n`-gon from its circumradius.
 * Params: n_vertex (vertex count), rad (circumradius).
 * Returns: Edge length.
 */
function ps_calc_edge(n_vertex, rad) = 2 * rad * sin(180 / n_vertex);

/**
 * Function: Compute the circumradius of a regular `n`-gon from its edge length.
 * Params: n_vertex (vertex count), edge_len (edge length).
 * Returns: Circumradius.
 */
function ps_calc_radius(n_vertex, edge_len) = edge_len / (2 * sin(180 / n_vertex));

function _ps_orient2(a, b, c) =
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);

// Intersection of two 2D lines given as dot(n, p) = d.
function _ps_line_intersect_2d(n1, d1, n2, d2, eps=1e-9) =
    let(det = n1[0] * n2[1] - n1[1] * n2[0])
    (abs(det) < eps)
        ? [0, 0]
        : [
            (d1 * n2[1] - n1[1] * d2) / det,
            (n1[0] * d2 - d1 * n2[0]) / det
        ];

function _ps_poly_turns2(poly) =
    let(n = len(poly))
    [
        for (i = [0:1:n-1])
            _ps_orient2(poly[(i - 1 + n) % n], poly[i], poly[(i + 1) % n])
    ];

function _ps_poly_is_convex2(poly, eps=1e-9) =
    let(
        turns = [for (t = _ps_poly_turns2(poly)) if (abs(t) > eps) t]
    )
    (len(turns) <= 1) ||
    (min(turns) >= -eps) ||
    (max(turns) <= eps);


///////////////////////////////////////
// ---- Geometry helpers ----
/**
 * Function: Compute the simple centroid of one face from its vertices.
 * Params: verts (vertex list), f (face index loop).
 * Returns: Arithmetic mean of the face vertices.
 * Limitations: Vertex-average centroid, not an area-weighted polygon centroid.
 */
function ps_face_centroid(verts, f) =
    len(f) == 0
        ? [0,0,0]
        : v_scale(v_sum([for (vid = f) verts[vid]]), 1 / len(f));

/**
 * Function: Compute the topological outward normal implied by one face's winding.
 * Params: verts (vertex list), f (face index loop).
 * Returns: Unit normal direction.
 * Limitations: Uses OpenSCAD's left-hand-rule winding convention (clockwise from outside).
 */
function ps_face_normal(verts, f) =
    // OpenSCAD expects LHR (clockwise from outside), so flip cross product.
    v_norm(v_cross(
        verts[f[2]] - verts[f[0]],
        verts[f[1]] - verts[f[0]]
    ));

/**
 * Function: Compute the placement normal used for face-local frames.
 * Params: verts (vertex list), f (face index loop), eps (degeneracy tolerance).
 * Returns: Unit face-frame normal.
 * Limitations: Non-planar faces use a Newell-style best-fit normal aligned to `ps_face_normal(...)`.
 */
function ps_face_frame_normal(verts, f, eps=1e-12) =
    let(
        n = len(f),
        nx = (n < 3) ? 0 : sum([
            for (i = [0:1:n-1])
                let(
                    j = (i + 1) % n,
                    pi = verts[f[i]],
                    pj = verts[f[j]]
                )
                (pi[1] - pj[1]) * (pi[2] + pj[2])
        ]),
        ny = (n < 3) ? 0 : sum([
            for (i = [0:1:n-1])
                let(
                    j = (i + 1) % n,
                    pi = verts[f[i]],
                    pj = verts[f[j]]
                )
                (pi[2] - pj[2]) * (pi[0] + pj[0])
        ]),
        nz = (n < 3) ? 0 : sum([
            for (i = [0:1:n-1])
                let(
                    j = (i + 1) % n,
                    pi = verts[f[i]],
                    pj = verts[f[j]]
                )
                (pi[0] - pj[0]) * (pi[1] + pj[1])
        ]),
        n_newell = [nx, ny, nz],
        n_topo = ps_face_normal(verts, f),
        n_raw = (norm(n_newell) > eps) ? n_newell : n_topo,
        n_aligned = (v_dot(n_raw, n_topo) < 0) ? [-n_raw[0], -n_raw[1], -n_raw[2]] : n_raw
    )
    v_norm(n_aligned);

// Magnitude of polygon area via triangle fan (works for planar faces).
function _ps_face_area_mag(verts, f) =
    (len(f) < 3) ? 0 :
    sum([
        for (i = [1:1:len(f)-2])
            norm(v_cross(verts[f[i]] - verts[f[0]], verts[f[i+1]] - verts[f[0]])) / 2
    ]);

// Max signed-distance deviation of face vertices from face plane.
function _ps_face_planarity_err(verts, f, eps=1e-12) =
    (len(f) < 3) ? 0 :
    let(
        n_raw = ps_face_normal(verts, f),
        n_len = norm(n_raw),
        n = (n_len <= eps) ? [0,0,1] : (n_raw / n_len),
        d = v_dot(n, verts[f[0]]),
        errs = [for (vi = f) abs(v_dot(n, verts[vi]) - d)]
    )
    (len(errs) == 0) ? 0 : max(errs);

function _ps_faces_max_planarity_err(verts, faces, eps=1e-12) =
    (len(faces) == 0) ? 0 : max([for (f = faces) _ps_face_planarity_err(verts, f, eps)]);
/**
 * Function: Return one adjacent vertex of the requested vertex.
 * Params: poly (poly descriptor), vi (vertex index).
 * Returns: Index of one neighboring vertex.
 * Limitations: Returns the first available neighbor only; not a full adjacency query.
 */
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
///////////////////////////////////////
// ---- Topology helpers ----
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


/**
 * Function: Build the incident-face table for a known edge list.
 * Params: faces (face index loops), edges (canonical undirected edge list).
 * Returns: One face-index list per edge.
 * Limitations: Closed manifolds usually yield two faces per edge; open shells may yield fewer.
 */
function ps_edge_faces_table(faces, edges) =
    [
        for (ei = [0 : len(edges)-1])
            let(e = edges[ei])
            [
                for (fi = [0 : len(faces)-1])
                    if (ps_face_has_edge(faces[fi], e[0], e[1])) fi
            ]
    ];

/**
 * Function: Test whether a face contains an undirected edge.
 * Params: f (face index loop), a,b (edge endpoints).
 * Returns: `true` when `{a, b}` appears along the face boundary.
 */
function ps_face_has_edge(f, a, b) =
    sum([
        for (k = [0 : len(f)-1])
            let(
                x = f[k],
                y = f[(k+1) % len(f)]
            )
            ((x==a && y==b) || (x==b && y==a)) ? 1 : 0
    ]) > 0;

/**
 * Function: Return the set of faces incident to one vertex.
 * Params: poly (poly descriptor), vi (vertex index).
 * Returns: Unordered incident face indices.
 */
function vertex_incident_faces(poly, vi) =
    let(faces = poly_faces(poly))
    [
        for (fi = [0 : len(faces)-1])
            let(f = faces[fi])
            if (sum([for (k = [0 : len(f)-1]) f[k] == vi ? 1 : 0]) > 0) fi
    ];

// Given current face around vertex v, pick the next incident face in the cycle.
function next_face_around_vertex(v, f_cur, f_prev, faces, edges, edge_faces) =
    let(
        f = faces[f_cur],
        n = len(f),
        pos = [for (k = [0 : n-1]) if (f[k] == v) k],
        k0 = pos[0],
        k_prev = (k0 - 1 + n) % n,
        k_next = (k0 + 1) % n,
        v_prev = f[k_prev],
        v_next = f[k_next],
        ei1 = ps_find_edge_index(edges, v, v_next),
        ei2 = ps_find_edge_index(edges, v_prev, v),
        ef1 = edge_faces[ei1],
        ef2 = edge_faces[ei2],
        cand1 = (ef1[0] == f_cur ? ef1[1] : ef1[0]),
        cand2 = (ef2[0] == f_cur ? ef2[1] : ef2[0]),
        candidates = [cand1, cand2],
        filtered = [for (cf = candidates) if (cf != f_prev) cf]
    )
    filtered[0];

function faces_around_vertex_rec(v, f_cur, f_prev, f_start, faces, edges, edge_faces, acc = []) =
    let(next = next_face_around_vertex(v, f_cur, f_prev, faces, edges, edge_faces))
    (next == f_start)
        ? concat(acc, [f_cur])
        : faces_around_vertex_rec(v, next, f_cur, f_start, faces, edges, edge_faces, concat(acc, [f_cur]));

/**
 * Function: Walk the incident faces around one vertex in cyclic order.
 * Params: poly (poly descriptor), v (vertex index), edges (edge list), edge_faces (incident-face table).
 * Returns: Incident face indices ordered around the vertex.
 * Limitations: Assumes a manifold neighborhood around the vertex.
 */
function faces_around_vertex(poly, v, edges, edge_faces) =
    let(
        faces = poly_faces(poly),
        inc = vertex_incident_faces(poly, v),
        start = inc[0]
    )
    faces_around_vertex_rec(v, start, -1, start, faces, edges, edge_faces);


///////////////////////////////////////
// ---- Frame/placement helpers ----
/**
 * Function: Compute the world-space center used for one face placement frame.
 * Params: poly (poly descriptor), fi (face index), scale (global vertex scale).
 * Returns: Face center in world coordinates.
 * Limitations: Vertex-average face center, matching the placement path.
 */
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

/**
 * Function: Compute the face-local +X axis used by `place_on_faces(...)`.
 * Params: poly (poly descriptor), fi (face index), scale (global vertex scale).
 * Returns: Unit +X direction in world coordinates.
 * Limitations: Projects into the face plane and falls back to the next face edge direction if needed.
 */
function poly_face_ex(poly, fi, scale) =
    let(f      = poly_faces(poly)[fi],
        vs     = poly_verts(poly),
        center = poly_face_center(poly, fi, scale),
        v0     = vs[f[0]] * scale,
        ez     = poly_face_ez(poly, fi, scale),
        ex_raw = v0 - center,
        // Keep face frame orthonormal even for non-planar faces:
        // project candidate x-axis into the local face plane.
        ex_proj = ex_raw - ez * v_dot(ex_raw, ez),
        ex_fallback_raw = (vs[f[1]] * scale) - (vs[f[0]] * scale),
        ex_fallback = ex_fallback_raw - ez * v_dot(ex_fallback_raw, ez))
    (norm(ex_proj) > 1e-12)
        ? v_norm(ex_proj)
        : v_norm(ex_fallback);   // local +X points towards face vertex order

/**
 * Function: Compute the face-local +Y axis used by `place_on_faces(...)`.
 * Params: poly (poly descriptor), fi (face index), scale (global vertex scale).
 * Returns: Unit +Y direction in world coordinates.
 */
function poly_face_ey(poly, fi, scale) =
    v_cross(
        poly_face_ez(poly, fi, scale),
        poly_face_ex(poly, fi, scale)
    );

/**
 * Function: Compute the face-local +Z axis used by `place_on_faces(...)`.
 * Params: poly (poly descriptor), fi (face index), scale (global vertex scale).
 * Returns: Unit +Z direction in world coordinates.
 * Limitations: Placement/frame normal, not necessarily the first-triangle normal for non-planar faces.
 */
function poly_face_ez(poly, fi, scale) =
    let(f  = poly_faces(poly)[fi],
        vs = poly_verts(poly),
        vs_scaled = [for (v = vs) v * scale])
    // Frame +Z for placement; best-fit for non-planar faces, aligned to LHR.
    ps_face_frame_normal(vs_scaled, f);

/**
 * Function: Build a homogeneous transform matrix from an orthonormal frame.
 * Params: center (frame origin), ex/ey/ez (frame basis vectors).
 * Returns: A 4x4 matrix suitable for `multmatrix(...)`.
 */
function ps_frame_matrix(center, ex, ey, ez) = [
    [ex[0], ey[0], ez[0], center[0]],
    [ex[1], ey[1], ez[1], center[1]],
    [ex[2], ey[2], ez[2], center[2]],
    [0,      0,     0,     1]
];
/**
 * Function: Orient one face so its winding points outward relative to its centroid.
 * Params: verts (vertex list), f (face index loop).
 * Returns: Either `f` or its reversal, whichever gives outward orientation.
 * Limitations: Uses the centroid-dot-normal heuristic; not a substitute for full validation.
 */
function ps_orient_face_outward(verts, f) =
    let(
        c = ps_face_centroid(verts, f),
        n = ps_face_normal(verts, f)
    )
    (v_dot(c, n) >= 0)
        ? f
        : _ps_reverse(f);  // reversed

/**
 * Function: Orient every face in a mesh outward using `ps_orient_face_outward(...)`.
 * Params: verts (vertex list), faces (face index loops).
 * Returns: Face list with outward-oriented winding.
 */
function ps_orient_all_faces_outward(verts, faces) =
    [ for (f = faces) ps_orient_face_outward(verts, f) ];
