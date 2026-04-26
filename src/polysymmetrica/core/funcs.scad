// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT


///////////////////////////////////////
// ---- Poly descriptor API ----
/**
 * Function: Return the vertex list from a poly descriptor.
 * Params: poly (`[verts, faces, e_over_ir]`)
 * Returns: vertex list
 */
function poly_verts(poly)      = poly[0];

/**
 * Function: Return the face index loops from a poly descriptor.
 * Params: poly (`[verts, faces, e_over_ir]`)
 * Returns: face list
 */
function poly_faces(poly)      = poly[1];

/**
 * Function: Return descriptor scale ratio `edge / inter-radius`.
 * Params: poly (`[verts, faces, e_over_ir]`)
 * Returns: positive scalar ratio
 */
function poly_e_over_ir(poly)  = poly[2];

/**
 * Function: Derive unique undirected edges from a poly descriptor.
 * Params: poly (`[verts, faces, e_over_ir]`)
 * Returns: edge list `[[a,b], ...]`; O(face-edge-count)
 */
function poly_edges(poly)      = _ps_edges_from_faces(poly_faces(poly));

/**
 * Function: Clamp a scalar to a closed interval.
 * Params: x (value), lo (lower bound), hi (upper bound)
 * Returns: `min(max(x, lo), hi)`
 */
function ps_clamp(x, lo, hi) = min(max(x, lo), hi);

/**
 * Function: Build a normalized poly descriptor from vertices and faces.
 * Params: verts (3D vertex list), faces (face index loops), e_over_ir (optional scale ratio)
 * Returns: `[centered_verts, faces, e_over_ir]`
 * Limitations/Gotchas: recenters by mean edge-midpoint center and computes default scale from minimum edge-midradius
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
 * Function: Make adjacent face windings consistent across shared edges.
 * Params: poly (poly descriptor)
 * Returns: poly descriptor with original vertices/scale and fixed face order
 * Limitations/Gotchas: fixes topological consistency only; it does not prove outward orientation
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
 * Function: Validate all face loops against a vertex list.
 * Params: verts (vertex list), faces (face index loops)
 * Returns: `true` when every face has at least 3 valid vertex indices
 */
function ps_faces_valid(verts, faces) =
    len([
        for (f = faces)
            if (len(f) >= 3 && ps_indices_in_range(f, len(verts)))
                1
    ]) == len(faces);

/**
 * Function: Check whether every index in a face is within `[0, max_idx)`.
 * Params: face (index list), max_idx (exclusive upper bound)
 * Returns: boolean
 */
function ps_indices_in_range(face, max_idx) =
    len([for (vi = face) if (vi >= 0 && vi < max_idx) 1]) == len(face);

/**
 * Function: Build successive pairs from a circular list.
 * Params: list (item list)
 * Returns: `[[list[i], list[i+1]], ... , [last, first]]`, or `[]` for lists shorter than 2
 */
function ps_cyclic_pairs(list) =
    let(n = len(list))
    (n < 2) ? [] :
    [ for (i = [0:1:n-1]) [list[i], list[(i+1)%n]] ];

///////////////////////////////////////
// ---- Winding/orientation helpers (private) ----
/**
 * Function: Return directed cyclic edges from a face loop.
 * Params: f (face index loop)
 * Returns: `[[a,b], ...]`
 */
function _ps_face_edges_dir(f) =
    ps_cyclic_pairs(f);

/**
 * Function: Determine whether a directed edge appears in a face.
 * Params: f (face index loop), a,b (edge endpoints)
 * Returns: `1` for `a->b`, `-1` for `b->a`, `0` when absent
 */
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

/**
 * Function: Find faces adjacent to an edge, excluding one face.
 * Params: faces (face list), a,b (edge endpoints), fi (face to exclude)
 * Returns: face indices containing undirected edge `{a,b}`
 */
function _ps_adjacent_faces_for_edge(faces, a, b, fi) =
    [
        for (fj = [0:1:len(faces)-1])
            if (fj != fi && _ps_face_edge_dir(faces[fj], a, b) != 0) fj
    ];

/**
 * Function: Replace one list element.
 * Params: list (input list), idx (target index), val (replacement value)
 * Returns: copy of `list` with `list[idx] = val`
 */
function _ps_list_set(list, idx, val) =
    [ for (i = [0:1:len(list)-1]) (i == idx) ? val : list[i] ];

/**
 * Function: Find the first undefined element in a list.
 * Params: list (input list)
 * Returns: first index with `is_undef(...)`, or `-1`
 */
function _ps_index_of_undef(list) =
    let(idx = [for (i = [0:1:len(list)-1]) if (is_undef(list[i])) i])
    (len(idx) == 0) ? -1 : idx[0];

/**
 * Function: Reverse a list.
 * Params: list (input list)
 * Returns: reversed list
 */
function _ps_reverse(list) =
    [ for (i = [len(list)-1 : -1 : 0]) list[i] ];

/**
 * Function: Build identity index map.
 * Params: n (size)
 * Returns: `[0, 1, ..., n-1]`
 */
function _ps_identity_map(n) =
    [for (i = [0:1:n-1]) i];

/**
 * Function: Count distinct values in a list.
 * Params: list (input list)
 * Returns: number of first occurrences
 */
function _ps_distinct_count(list) =
    len([for (i = [0:1:len(list)-1]) if (_ps_index_of(list, list[i]) == i) 1]);

/**
 * Function: Remove adjacent duplicate vertex ids from a cyclic face.
 * Params: f (face index loop)
 * Returns: face loop without adjacent repeats
 */
function _ps_face_strip_adjacent_dups(f) =
    let(n = len(f))
    (n == 0) ? [] :
    [for (i = [0:1:n-1]) if (f[i] != f[(i-1+n)%n]) f[i]];

/**
 * Function: Remove duplicated closing vertex from a face loop.
 * Params: f (face index loop)
 * Returns: `f` without final element when `last == first`
 */
function _ps_face_trim_closing_dup(f) =
    (len(f) >= 2 && f[0] == f[len(f)-1]) ? [for (i = [0:1:len(f)-2]) f[i]] : f;

/**
 * Function: Normalize one face loop by removing trivial duplicate vertices.
 * Params: f (face index loop)
 * Returns: cleaned face loop
 */
function _ps_face_clean_cycle(f) =
    _ps_face_trim_closing_dup(_ps_face_strip_adjacent_dups(f));

/**
 * Function: Clean every face loop in a face list.
 * Params: faces (face list)
 * Returns: cleaned face list
 */
function _ps_faces_clean_cycles(faces) =
    [for (f = faces) _ps_face_clean_cycle(f)];

/**
 * Function: Remap all vertex ids in face loops.
 * Params: faces (face list), old_to_new (index map)
 * Returns: remapped face list
 */
function _ps_faces_remap(faces, old_to_new) =
    [for (f = faces) [for (vi = f) old_to_new[vi]]];

/**
 * Function: Look up a pending neighbor face assignment.
 * Params: neighbors (`[[idx, value], ...]`), idx (face index)
 * Returns: assigned value or `undef`
 */
function _ps_neighbor_face(neighbors, idx) =
    let(hit = [for (p = neighbors) if (p[0] == idx) p[1]])
    (len(hit) == 0) ? undef : hit[0];

/**
 * Function: Apply pending neighbor face assignments into a fixed-face list.
 * Params: fixed (possibly-undef face list), neighbors (`[[idx, face], ...]`)
 * Returns: updated fixed list
 */
function _ps_apply_neighbors(fixed, neighbors) =
    [
        for (i = [0:1:len(fixed)-1])
            is_undef(fixed[i]) ? _ps_neighbor_face(neighbors, i) : fixed[i]
    ];

/**
 * Function: Extract newly assigned face indices.
 * Params: fixed (previous fixed list), neighbors (`[[idx, face], ...]`)
 * Returns: indices whose previous value was `undef`
 */
function _ps_new_neighbor_indices(fixed, neighbors) =
    [ for (p = neighbors) if (is_undef(fixed[p[0]])) p[0] ];

/**
 * Function: Breadth-first face-winding propagation over connected faces.
 * Params: faces (source face list), fixed (assigned oriented faces), queue (face indices)
 * Returns: fixed face list for the connected component reachable from queue
 */
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

/**
 * Function: Fix winding across all connected face components.
 * Params: faces (face list), fixed (optional in-progress oriented faces)
 * Returns: oriented face list with each shared edge opposite-directed
 */
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
/**
 * Function: Test whether a list contains a value.
 * Params: list (input list), v (value)
 * Returns: boolean
 */
function _ps_list_contains(list, v) =
    len(search(v, list)) > 0;

/**
 * Function: Find first index of a value.
 * Params: list (input list), v (value)
 * Returns: first index, or `-1`
 */
function _ps_index_of(list, v) =
    let(idx = search(v, list))
    (len(idx) == 0) ? -1 : idx[0];

///////////////////////////////////////
// ---- Vector math ----
/**
 * Function: Add two vectors.
 * Params: a,b (vectors)
 * Returns: `a + b`
 */
function v_add(a, b)   = a + b;

/**
 * Function: Subtract two vectors.
 * Params: a,b (vectors)
 * Returns: `a - b`
 */
function v_sub(a, b)   = a - b;

/**
 * Function: Scale a vector.
 * Params: a (vector), k (scalar)
 * Returns: `a * k`
 */
function v_scale(a, k) = a * k;             // scalar multiplication

/**
 * Function: Dot product.
 * Params: a,b (vectors)
 * Returns: scalar dot product
 */
function v_dot(a, b)   = a * b;             // dot product

/**
 * Function: Cross product.
 * Params: a,b (3D vectors)
 * Returns: `cross(a,b)`
 */
function v_cross(a, b) = cross(a, b);       // OpenSCAD built-in

/**
 * Function: Vector length.
 * Params: a (vector)
 * Returns: Euclidean norm
 */
function v_len(a)      = norm(a);           // built-in length

/**
 * Function: Normalize a vector.
 * Params: a (vector)
 * Returns: unit vector, or zero-like vector when input length is zero
 */
function v_norm(a)     = let(L = norm(a)) (L == 0 ? [0,0,0] : a / L);

/**
 * Function: Sort two endpoint indices into canonical undirected-edge order.
 * Params: a,b (indices)
 * Returns: `[min,max]`
 */
function _ps_ordered_pair(a, b) = (a < b) ? [a,b] : [b,a];


///////////////////////////////////////
// ---- Edge/list primitives ----
/**
 * Function: Test equality of ordered edge records.
 * Params: e1,e2 (`[a,b]`)
 * Returns: boolean
 */
function edge_equal(e1, e2) = (e1[0] == e2[0] && e1[1] == e2[1]);

/**
 * Function: Sum scalar list entries recursively.
 * Params: a (number list), i (start index)
 * Returns: scalar sum from `i` to end
 */
function sum(a, i = 0) =
    i >= len(a) ? 0 : a[i] + sum(a, i + 1);

/**
 * Function: Sum a list of equal-dimension vectors.
 * Params: list (vector list)
 * Returns: component-wise vector sum; `[]` for empty input
 */
function v_sum(list) =
    (len(list) == 0) ? [] :
    let(n = len(list[0]))
    [ for (i = [0:1:n-1]) sum([for (v = list) v[i]]) ];

/**
 * Function: Compute the simple centroid of a 2D point list.
 * Params: points (2D point list)
 * Returns: centroid `[x, y]`, or `[0, 0]` for an empty list
 */
function ps_centroid2d(points) =
    (len(points) == 0) ? [0, 0] :
    v_scale(v_sum(points), 1 / len(points));

/**
 * Function: Compute the midpoint of a 2D segment.
 * Params: seg2d (`[[x0,y0],[x1,y1]]`)
 * Returns: midpoint `[x, y]`
 */
function ps_segment_midpoint2d(seg2d) =
    [(seg2d[0][0] + seg2d[1][0]) / 2, (seg2d[0][1] + seg2d[1][1]) / 2];

/**
 * Function: Project points to the XY plane.
 * Params: points (2D/3D/ND point list)
 * Returns: `[[x, y], ...]`
 */
function ps_xy(points) =
    [for (p = points) [p[0], p[1]]];

// Rotate vector v around axis by ang (degrees).
/**
 * Function: Rotate a vector around an axis.
 * Params: v (3D vector), axis (3D axis vector), ang (degrees)
 * Returns: rotated vector
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
 * Function: Find an undirected edge in an edge list.
 * Params: edges (canonical edge list), a,b (edge endpoints)
 * Returns: edge index
 * Limitations/Gotchas: assumes the edge exists
 */
function ps_find_edge_index(edges, a, b) =
    let(
        e = _ps_ordered_pair(a, b),
        idxs = [for (i = [0 : len(edges)-1]) if (edge_equal(edges[i], e)) i]
    )
    idxs[0];   // assume the edge exists

/**
 * Function: Compare points with tolerance.
 * Params: p,q (points), eps (distance tolerance)
 * Returns: boolean
 */
function ps_point_eq(p,q,eps) = norm(p-q) <= eps;

/**
 * Function: Minimum scalar in a list.
 * Params: list (number list), i (scan index), cur (current minimum)
 * Returns: minimum value, or `undef` for empty input
 */
function _ps_list_min(list, i=0, cur=undef) =
    (i >= len(list)) ? cur :
    let(v = list[i])
    _ps_list_min(list, i+1, is_undef(cur) ? v : (v < cur ? v : cur));

/**
 * Function: Remove the first matching value from a list.
 * Params: list (input list), v (value), i (scan index)
 * Returns: list with first occurrence removed
 */
function _ps_remove_first(list, v, i=0) =
    (i >= len(list)) ? [] :
    (list[i] == v) ? [for (j = [i+1:1:len(list)-1]) list[j]]
                  : concat([list[i]], _ps_remove_first(list, v, i+1));

/**
 * Function: Sort a scalar list by selection recursion.
 * Params: list (input list), acc (accumulator)
 * Returns: ascending sorted list
 */
function _ps_sort(list, acc=[]) =
    (len(list) == 0) ? acc :
    let(mn = _ps_list_min(list))
    _ps_sort(_ps_remove_first(list, mn), concat(acc, [mn]));

///////////////////////////////////////
// ---- Polygon/polygram helpers ----
/**
 * Function: Compute Euclidean gcd for integer-like values.
 * Params: a,b (numbers)
 * Returns: non-negative integer gcd after rounding
 */
function _ps_gcd(a, b) =
    let(ai = abs(round(a)), bi = abs(round(b)))
    (bi == 0) ? ai : _ps_gcd(bi, ai % bi);

/**
 * Function: Validate Schläfli-like polygon/polygram parameters.
 * Params: n (vertex count), p (step), who (caller label)
 * Returns: rounded `[n, p]`
 * Limitations/Gotchas: requires `n,p` coprime so the polygram is a single cycle
 */
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

/**
 * Function: Circumradius for regular/star polygon `{n,p}`.
 * Params: n (vertex count), p (step), edge (chord length)
 * Returns: radius scalar
 */
function _ps_polygram_radius(n, p, edge) =
    edge / (2 * sin(180 * p / n));

/**
 * Function: Circumradius for regular polygon `{n,1}`.
 * Params: n (vertex count), edge (edge length)
 * Returns: radius scalar
 */
function _ps_ngon_radius(n, edge) =
    _ps_polygram_radius(n, 1, edge);

/**
 * Function: Vertex order for a single-cycle polygram `{n,p}`.
 * Params: n (vertex count), p (step)
 * Returns: index cycle `[(k*p)%n, ...]`
 */
function _ps_polygram_cycle(n, p) =
    [for (k = [0:1:n-1]) (k * p) % n];

/**
 * Function: Generate a regular support ring at fixed Z.
 * Params: n (vertex count), radius (ring radius), z (Z coordinate), phase (degrees)
 * Returns: 3D point ring
 */
function _ps_ngon_ring(n, radius, z, phase=0) =
    [for (k = [0:1:n-1]) [radius * cos(360 * k / n + phase), radius * sin(360 * k / n + phase), z]];

/**
 * Function: Compute mean edge-midpoint center for a mesh.
 * Params: verts (3D vertices), faces (face loops)
 * Returns: center point
 */
function _ps_poly_mid_center(verts, faces) =
    let(
        edges = _ps_edges_from_faces(faces),
        mids = [for (e = edges) (verts[e[0]] + verts[e[1]]) / 2],
        _ok = assert(len(mids) > 0, "poly: edge-midpoint center requires at least one edge")
    )
    v_scale(v_sum(mids), 1 / len(mids));

/**
 * Function: Compute inter-radius from minimum centered edge-midpoint radius.
 * Params: verts (3D vertices), faces (face loops)
 * Returns: positive inter-radius
 */
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
/**
 * Function: Determinant of a 3x3 matrix.
 * Params: m (3x3 matrix)
 * Returns: determinant scalar
 */
function _ps_det3(m) =
    m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
    m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
    m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

/**
 * Function: Replace one column in a 3x3 matrix.
 * Params: m (3x3 matrix), col (column index 0..2), b (replacement column)
 * Returns: updated 3x3 matrix
 */
function _ps_replace_col(m, col, b) =
    [
        [ col == 0 ? b[0] : m[0][0], col == 1 ? b[0] : m[0][1], col == 2 ? b[0] : m[0][2] ],
        [ col == 0 ? b[1] : m[1][0], col == 1 ? b[1] : m[1][1], col == 2 ? b[1] : m[1][2] ],
        [ col == 0 ? b[2] : m[2][0], col == 1 ? b[2] : m[2][1], col == 2 ? b[2] : m[2][2] ]
    ];

/**
 * Function: Solve a 3x3 linear system with Cramer's rule.
 * Params: m (3x3 matrix), b (right-hand vector), eps (singularity tolerance)
 * Returns: solution `[x0,x1,x2]`, or `undef` for near-singular systems
 */
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
/**
 * Function: Build prefix offsets from a count list.
 * Params: counts (non-negative counts), acc (offset accumulator)
 * Returns: offsets such as `[0, c0, c0+c1, ...]` when seeded with `[0]`
 */
function _ps_prefix_offsets(counts, acc=[]) =
    (len(counts) == 0) ? acc :
    let(last = (len(acc) == 0) ? 0 : acc[len(acc)-1])
    _ps_prefix_offsets(
        [for (i = [1:1:len(counts)-1]) counts[i]],
        concat(acc, [last + counts[0]])
    );

/**
 * Function: Prefix offsets for concatenated face vertices.
 * Params: faces (face list)
 * Returns: offsets by face
 */
function _ps_face_offsets(faces) =
    let(counts = [for (f = faces) len(f)])
    _ps_prefix_offsets(counts, [0]);

/**
 * Function: Prefix offsets for concatenated directed face edges.
 * Params: faces (face list)
 * Returns: offsets by face, counting two directed half-edges per edge
 */
function _ps_face_edge_offsets(faces) =
    let(counts = [for (f = faces) 2 * len(f)])
    _ps_prefix_offsets(counts, [0]);


///////////////////////////////////////
// ---- Polygon helpers ----
/**
 * Function: Compute regular polygon edge length from radius.
 * Params: n_vertex (vertex count), rad (circumradius)
 * Returns: edge length
 */
function ps_calc_edge(n_vertex, rad) = 2 * rad * sin(180 / n_vertex);

/**
 * Function: Compute regular polygon circumradius from edge length.
 * Params: n_vertex (vertex count), edge_len (edge length)
 * Returns: radius scalar
 */
function ps_calc_radius(n_vertex, edge_len) = edge_len / (2 * sin(180 / n_vertex));


///////////////////////////////////////
// ---- Geometry helpers ----
/**
 * Function: Compute mean vertex centroid for one face.
 * Params: verts (3D vertex list), f (face index loop)
 * Returns: 3D centroid, or `[0,0,0]` for empty face
 */
function ps_face_centroid(verts, f) =
    len(f) == 0
        ? [0,0,0]
        : v_scale(v_sum([for (vid = f) verts[vid]]), 1 / len(f));

/**
 * Function: Compute topological face normal using OpenSCAD LHR winding.
 * Params: verts (3D vertex list), f (face index loop)
 * Returns: unit normal direction
 * Limitations/Gotchas: uses first three vertices; use `ps_face_frame_normal(...)` for non-planar placement frames
 */
function ps_face_normal(verts, f) =
    // OpenSCAD expects LHR (clockwise from outside), so flip cross product.
    v_norm(v_cross(
        verts[f[2]] - verts[f[0]],
        verts[f[1]] - verts[f[0]]
    ));

/**
 * Function: Compute placement frame normal for a face.
 * Params: verts (3D vertex list), f (face index loop), eps (degeneracy tolerance)
 * Returns: unit normal direction aligned with `ps_face_normal(...)`
 * Limitations/Gotchas: uses Newell-style best-fit normal for non-planar faces
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

/**
 * Function: Compute face area magnitude by triangle fan.
 * Params: verts (3D vertex list), f (face index loop)
 * Returns: non-negative area
 * Limitations/Gotchas: intended for planar faces
 */
function _ps_face_area_mag(verts, f) =
    (len(f) < 3) ? 0 :
    sum([
        for (i = [1:1:len(f)-2])
            norm(v_cross(verts[f[i]] - verts[f[0]], verts[f[i+1]] - verts[f[0]])) / 2
    ]);

/**
 * Function: Compute maximum vertex deviation from a face plane.
 * Params: verts (3D vertex list), f (face index loop), eps (normal tolerance)
 * Returns: maximum absolute signed-distance error
 */
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

/**
 * Function: Compute maximum planarity error over a face list.
 * Params: verts (3D vertex list), faces (face list), eps (normal tolerance)
 * Returns: maximum face planarity error
 */
function _ps_faces_max_planarity_err(verts, faces, eps=1e-12) =
    (len(faces) == 0) ? 0 : max([for (f = faces) _ps_face_planarity_err(verts, f, eps)]);


/**
 * Function: Return one neighbor vertex of a vertex.
 * Params: poly (poly descriptor), vi (vertex index)
 * Returns: first adjacent vertex found in face traversal
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


///////////////////////////////////////
// ---- Topology helpers ----
/**
 * Function: Build unique undirected edges from a face list.
 * Params: faces (face index loops)
 * Returns: canonical edge list `[[min,max], ...]`
 */
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
 * Function: Map each edge to incident faces.
 * Params: faces (face list), edges (canonical edge list)
 * Returns: `[[face_idx, ...], ...]` matching `edges`
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
 * Params: f (face index loop), a,b (edge endpoints)
 * Returns: boolean
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
 * Function: Find faces incident to a vertex.
 * Params: poly (poly descriptor), vi (vertex index)
 * Returns: unordered face-index list
 */
function vertex_incident_faces(poly, vi) =
    let(faces = poly_faces(poly))
    [
        for (fi = [0 : len(faces)-1])
            let(f = faces[fi])
            if (sum([for (k = [0 : len(f)-1]) f[k] == vi ? 1 : 0]) > 0) fi
    ];

/**
 * Function: Pick next incident face around a vertex.
 * Params: v (vertex index), f_cur (current face), f_prev (previous face), faces/edges/edge_faces (topology tables)
 * Returns: next face index in the local fan
 * Limitations/Gotchas: assumes manifold local topology
 */
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

/**
 * Function: Recursively traverse incident faces around a vertex.
 * Params: v (vertex index), f_cur/f_prev/f_start (face indices), faces/edges/edge_faces (topology tables), acc (accumulator)
 * Returns: cyclically ordered face-index list
 */
function faces_around_vertex_rec(v, f_cur, f_prev, f_start, faces, edges, edge_faces, acc = []) =
    let(next = next_face_around_vertex(v, f_cur, f_prev, faces, edges, edge_faces))
    (next == f_start)
        ? concat(acc, [f_cur])
        : faces_around_vertex_rec(v, next, f_cur, f_start, faces, edges, edge_faces, concat(acc, [f_cur]));

/**
 * Function: Return incident faces around a vertex in cyclic order.
 * Params: poly (poly descriptor), v (vertex index), edges (edge list), edge_faces (edge-face table)
 * Returns: ordered face-index list
 * Limitations/Gotchas: assumes manifold vertex neighborhood
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
 * Function: Compute scaled face center for placement.
 * Params: poly (poly descriptor), fi (face index), scale (scale factor)
 * Returns: 3D center point
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
 * Function: Compute local face-frame X axis for placement.
 * Params: poly (poly descriptor), fi (face index), scale (scale factor)
 * Returns: unit 3D X-axis vector in world/poly coordinates
 * Limitations/Gotchas: candidate axis is projected into the face plane to remain orthonormal for non-planar faces
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
 * Function: Compute local face-frame Y axis for placement.
 * Params: poly (poly descriptor), fi (face index), scale (scale factor)
 * Returns: unit 3D Y-axis vector
 */
function poly_face_ey(poly, fi, scale) =
    v_cross(
        poly_face_ez(poly, fi, scale),
        poly_face_ex(poly, fi, scale)
    );


/**
 * Function: Compute local face-frame Z axis for placement.
 * Params: poly (poly descriptor), fi (face index), scale (scale factor)
 * Returns: unit 3D Z-axis vector
 */
function poly_face_ez(poly, fi, scale) =
    let(f  = poly_faces(poly)[fi],
        vs = poly_verts(poly),
        vs_scaled = [for (v = vs) v * scale])
    // Frame +Z for placement; best-fit for non-planar faces, aligned to LHR.
    ps_face_frame_normal(vs_scaled, f);


/**
 * Function: Build a 4x4 transform matrix from frame axes and center.
 * Params: center (translation), ex/ey/ez (basis vectors)
 * Returns: OpenSCAD `multmatrix`-compatible matrix
 */
function ps_frame_matrix(center, ex, ey, ez) = [
    [ex[0], ey[0], ez[0], center[0]],
    [ex[1], ey[1], ez[1], center[1]],
    [ex[2], ey[2], ez[2], center[2]],
    [0,      0,     0,     1]
];


/**
 * Function: Orient one face so its normal points away from origin.
 * Params: verts (3D vertex list), f (face index loop)
 * Returns: original or reversed face loop
 * Limitations/Gotchas: origin-relative heuristic; appropriate for centered radial polyhedra
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
 * Function: Orient all faces so normals point away from origin.
 * Params: verts (3D vertex list), faces (face list)
 * Returns: oriented face list
 */
function ps_orient_all_faces_outward(verts, faces) =
    [ for (f = faces) ps_orient_face_outward(verts, f) ];
