use <funcs.scad>

// Polyhedral element classification by families (faces / edges / vertices).
// detail=0: topology only
// detail>=1: adds neighbour refinement; geometry is controlled by include_geom.

/**
 * Quantize numeric value to epsilon grid.
 *
 * Purpose:
 * - stabilize geometric keys against floating point noise.
 *
 * Args:
 * - `v`: scalar value.
 * - `eps`: quantization step (<=0 disables quantization).
 *
 * Returns:
 * - rounded scalar.
 */
function _ps_round(v, eps) =
    (eps <= 0) ? v : (round(v / eps) * eps);

/**
 * Compute Euclidean length of an edge.
 *
 * Args:
 * - `verts`: vertex array.
 * - `e`: edge pair `[a,b]`.
 *
 * Returns:
 * - edge length.
 */
function _ps_edge_len(verts, e) =
    norm(verts[e[0]] - verts[e[1]]);

/**
 * Keep first occurrence of each key (stable order).
 *
 * Approach:
 * - linear first-seen filter preserving input order.
 *
 * Args:
 * - `keys`: list of comparable values.
 *
 * Returns:
 * - list of unique keys.
 */
function _ps_unique_keys(keys) =
    [
        for (i = [0:1:len(keys)-1])
            if (sum([for (j = [0:1:i-1]) keys[j] == keys[i] ? 1 : 0]) == 0)
                keys[i]
    ];

/**
 * Convert a key to canonical comparable string.
 *
 * Args:
 * - `k`: key value (possibly nested list).
 *
 * Returns:
 * - string representation.
 */
function _ps_key_str(k) = str(k);

/**
 * Group element indices by identical key.
 *
 * Approach:
 * - stringify keys, unique string set, then collect member indices.
 *
 * Args:
 * - `keys`: per-element key list.
 *
 * Returns:
 * - family list `[[key, idxs], ...]`.
 */
function _ps_group_by_key(keys) =
    let(
        key_strs = [for (k = keys) _ps_key_str(k)],
        uniq_strs = _ps_unique_keys(key_strs)
    )
    [
        for (ks = uniq_strs)
            let(idxs = [for (i = [0:1:len(keys)-1]) if (key_strs[i] == ks) i])
                [keys[idxs[0]], idxs]
    ];

/**
 * Build direct element->family-id map.
 *
 * Args:
 * - `n`: element count.
 * - `fams`: family list `[[key, idxs], ...]`.
 *
 * Returns:
 * - ids array of length `n`; `-1` if unmatched.
 */
// Map each element index -> family id.
function _ps_family_ids(n, fams) =
    [
        for (i = [0:1:n-1])
            let(hit = [for (fi = [0:1:len(fams)-1]) if (search(i, fams[fi][1]) != []) fi])
            (len(hit) == 0) ? -1 : hit[0]
    ];

/**
 * Cyclic left rotation helper.
 *
 * Args:
 * - `list`: sequence.
 * - `k`: rotation offset.
 *
 * Returns:
 * - rotated sequence.
 */
function _ps_rotate(list, k) =
    let(n = len(list))
    [for (i = [0:1:n-1]) list[(i + k) % n]];

/**
 * Lexicographic strict-greater comparison.
 *
 * Args:
 * - `a`, `b`: equal-length sequences.
 *
 * Returns:
 * - true iff `a > b` lexicographically.
 */
function _ps_lex_gt(a, b) =
    (len(a) == 0) ? false :
    (a[0] > b[0]) ? true :
    (a[0] < b[0]) ? false :
    _ps_lex_gt([for (i = [1:1:len(a)-1]) a[i]], [for (i = [1:1:len(b)-1]) b[i]]);

/**
 * Lexicographic strict-less comparison.
 *
 * Args:
 * - `a`, `b`: equal-length sequences.
 *
 * Returns:
 * - true iff `a < b` lexicographically.
 */
function _ps_lex_lt(a, b) =
    (len(a) == 0) ? false :
    (a[0] < b[0]) ? true :
    (a[0] > b[0]) ? false :
    _ps_lex_lt([for (i = [1:1:len(a)-1]) a[i]], [for (i = [1:1:len(b)-1]) b[i]]);

/**
 * Pick lexicographic extremum from a sequence set.
 *
 * Args:
 * - `list`: candidate sequences.
 * - `mode`: `"max"` or `"min"`.
 * - `acc`, `idx`: recursion internals.
 *
 * Returns:
 * - best sequence by chosen order.
 */
function _ps_lex_best(list, mode="max", acc=undef, idx=0) =
    (len(list) == 0) ? [] :
    (idx >= len(list)) ? acc :
    let(
        cur = list[idx],
        best = is_undef(acc)
            ? cur
            : ((mode == "min") ? (_ps_lex_lt(cur, acc) ? cur : acc)
                               : (_ps_lex_gt(cur, acc) ? cur : acc))
    )
    _ps_lex_best(list, mode, best, idx + 1);

/**
 * Canonical cyclic representation (necklace/bracelet form).
 *
 * Approach:
 * - enumerate rotations (and optional reflected rotations),
 *   then pick lexicographic extremum.
 *
 * Args:
 * - `seq`: cyclic sequence.
 * - `mode`: `"max"` or `"min"` for canonical representative.
 * - `allow_reflect`: include reversed sequence rotations.
 *
 * Returns:
 * - canonicalized sequence.
 */
// Canonical cyclic rotation (necklace). If allow_reflect, treats reversed
// sequences as equivalent (bracelet).
function _ps_cyclic_canonical(seq, mode="max", allow_reflect=false) =
    let(
        n = len(seq),
        rots = [for (i = [0:1:n-1]) _ps_rotate(seq, i)],
        rseq = allow_reflect ? [for (i = [0:1:n-1]) _ps_rotate(_ps_reverse(seq), i)] : [],
        all = concat(rots, rseq),
        best = _ps_lex_best(all, mode)
    )
    best;

/**
 * Build base face keys.
 *
 * Args:
 * - `verts`, `faces`: mesh.
 * - `detail`: 0 => `[n]`; >0 => `[n, avg_edge_len_rounded]`.
 * - `eps`: geometric rounding epsilon.
 *
 * Returns:
 * - per-face key list.
 */
function _ps_face_keys_from(verts, faces, detail, eps) =
    [
        for (f = faces)
            let(
                n = len(f),
                lens = [
                    for (k = [0:1:n-1])
                        _ps_edge_len(verts, [f[k], f[(k+1)%n]])
                ],
                avg_len = (n == 0) ? 0 : (sum(lens) / n)
            )
            (detail <= 0)
                ? [n]
                : [n, _ps_round(avg_len, eps)]
    ];

/**
 * Build base edge keys.
 *
 * Args:
 * - `verts`, `faces`, `edges`, `edge_faces`: mesh + adjacency.
 * - `detail`: 0 => `[adj_face_n0, adj_face_n1]`;
 *             >0 => append rounded edge length.
 * - `eps`: geometric rounding epsilon.
 *
 * Returns:
 * - per-edge key list.
 */
function _ps_edge_keys_from(verts, faces, edges, edge_faces, detail, eps) =
    [
        for (ei = [0:1:len(edges)-1])
            let(
                e = edges[ei],
                fpair = edge_faces[ei],
                k0 = len(faces[fpair[0]]),
                k1 = len(faces[fpair[1]]),
                ks = (k0 < k1) ? [k0, k1] : [k1, k0],
                el = _ps_edge_len(verts, e)
            )
            (detail <= 0)
                ? [ks[0], ks[1]]
                : [ks[0], ks[1], _ps_round(el, eps)]
    ];

/**
 * Encode undirected edges as sortable scalar keys.
 *
 * Args:
 * - `edges`: edge pairs `[a,b]`.
 * - `nv`: vertex count (base for pair packing).
 *
 * Returns:
 * - numeric key list.
 */
function _ps_edge_keys_list(edges, nv) =
    [for (e = edges) min(e[0], e[1]) * nv + max(e[0], e[1])];

/**
 * Lookup edge index by packed numeric key.
 *
 * Args:
 * - `edge_keys`: output of `_ps_edge_keys_list`.
 * - `a`, `b`: edge endpoints.
 * - `nv`: vertex count used for packing.
 *
 * Returns:
 * - edge index.
 */
function _ps_edge_index(edge_keys, a, b, nv) =
    let(e = min(a, b) * nv + max(a, b), idxs = search(e, edge_keys))
    idxs[0];

/**
 * Build base vertex keys from valence and cyclic incident-face sizes.
 *
 * Args:
 * - `verts`, `faces`, `edges`, `edge_faces`: mesh + adjacency.
 * - `detail`: 0 => `[valence, cyclic_face_sizes...]`;
 *             >0 => append rounded mean incident edge length.
 * - `eps`: geometric rounding epsilon.
 *
 * Returns:
 * - per-vertex key list.
 */
function _ps_vert_keys_from(verts, faces, edges, edge_faces, detail, eps) =
    let(
        valences = [
            for (vi = [0:1:len(verts)-1])
                sum([for (e = edges) (e[0] == vi || e[1] == vi) ? 1 : 0])
        ]
    )
    [
        for (vi = [0:1:len(verts)-1])
            let(
                fc = faces_around_vertex([verts, faces, 0], vi, edges, edge_faces),
                ks = _ps_cyclic_canonical([for (fi = fc) len(faces[fi])], "max", false),
                elens = [
                    for (e = edges)
                        if (e[0] == vi || e[1] == vi)
                            _ps_edge_len(verts, e)
                ],
                avg_len = (len(elens) == 0) ? 0 : (sum(elens) / len(elens))
            )
            (detail <= 0)
                ? concat([valences[vi]], ks)
                : concat([valences[vi]], ks, [_ps_round(avg_len, eps)])
    ];

/**
 * Refine face keys by neighboring face-family and vertex-family context.
 *
 * Approach:
 * - convert neighbor keys into temporary family ids,
 * - collect cyclic neighbor face ids around each edge of face,
 * - collect cyclic vertex family ids around each face,
 * - append canonicalized signatures to existing key.
 *
 * Args:
 * - `poly`: descriptor.
 * - `face_keys`: current face keys being refined.
 * - `face_nbr_keys`: keys used to map neighboring faces to temporary ids.
 * - `vert_ids`: vertex family id map.
 * - `edges`, `edge_faces`, `edge_keys`, `nv`: adjacency/lookup support.
 *
 * Returns:
 * - refined face key list.
 */
// Refine face keys using neighboring face topology (by face_nbr_keys).
function _ps_refine_face_keys(poly, face_keys, face_nbr_keys, vert_ids, edges, edge_faces, edge_keys, nv) =
    let(
        faces = poly_faces(poly),
        face_nbr_fams = _ps_group_by_key(face_nbr_keys),
        face_nbr_ids = _ps_family_ids(len(faces), face_nbr_fams)
    )
    [
        for (fi = [0:1:len(faces)-1])
            let(
                f = faces[fi],
                nbr_ids = [
                    for (k = [0:1:len(f)-1])
                        let(
                            a = f[k],
                            b = f[(k+1)%len(f)],
                            ei = _ps_edge_index(edge_keys, a, b, nv),
                            fpair = edge_faces[ei],
                            f_other = (fpair[0] == fi) ? fpair[1] : fpair[0]
                        )
                        face_nbr_ids[f_other]
                ],
                v_ids = [for (vi = f) vert_ids[vi]],
                // Face neighborhoods are cyclic up to rotation and reflection.
                nbr_key = _ps_cyclic_canonical(nbr_ids, "max", true),
                v_key = _ps_cyclic_canonical(v_ids, "max", true)
            )
            concat(face_keys[fi], [nbr_key, v_key])
    ];

/**
 * Refine edge keys using current face/vertex family ids.
 *
 * Args:
 * - `poly`: descriptor.
 * - `edge_keys`: current edge keys.
 * - `face_keys`, `vert_keys`: current face/vertex keys.
 * - `edges`, `edge_faces`: adjacency.
 *
 * Returns:
 * - refined edge key list.
 */
function _ps_refine_edge_keys(poly, edge_keys, face_keys, vert_keys, edges, edge_faces) =
    let(
        faces = poly_faces(poly),
        face_fams = _ps_group_by_key(face_keys),
        vert_fams = _ps_group_by_key(vert_keys),
        face_ids = _ps_family_ids(len(faces), face_fams),
        vert_ids = _ps_family_ids(len(poly_verts(poly)), vert_fams)
    )
    [
        for (ei = [0:1:len(edges)-1])
            let(
                e = edges[ei],
                fpair = edge_faces[ei],
                f0 = face_ids[fpair[0]],
                f1 = face_ids[fpair[1]],
                v0 = vert_ids[e[0]],
                v1 = vert_ids[e[1]],
                vk = (v0 < v1) ? [v0, v1] : [v1, v0],
                fk = (f0 < f1) ? [f0, f1] : [f1, f0]
            )
            concat(edge_keys[ei], [fk, vk])
    ];

/**
 * Refine vertex keys using incident face-family ids.
 *
 * Args:
 * - `poly`: descriptor.
 * - `vert_keys`: current vertex keys.
 * - `face_keys`: current face keys.
 * - `edges`, `edge_faces`: adjacency.
 *
 * Returns:
 * - refined vertex key list.
 */
function _ps_refine_vert_keys(poly, vert_keys, face_keys, edges, edge_faces) =
    let(
        faces = poly_faces(poly),
        face_fams = _ps_group_by_key(face_keys),
        face_ids = _ps_family_ids(len(faces), face_fams)
    )
    [
        for (vi = [0:1:len(poly_verts(poly))-1])
            let(
                fc = faces_around_vertex(poly, vi, edges, edge_faces),
                nbr_ids = [for (fi = fc) face_ids[fi]],
                nbr_key = _ps_cyclic_canonical(nbr_ids, "max", false)
            )
            concat(vert_keys[vi], [nbr_key])
    ];

/**
 * Exact elementwise key list equality.
 *
 * Args:
 * - `a`, `b`: key lists.
 *
 * Returns:
 * - true when lengths and all elements are equal.
 */
function _ps_keys_equal(a, b) =
    (len(a) == len(b)) ? (min([for (i = [0:1:len(a)-1]) a[i] == b[i] ? 1 : 0]) == 1) : false;

/**
 * Compare partition equality (family membership), ignoring key payload growth.
 *
 * Args:
 * - `keys_a`, `keys_b`: key lists.
 *
 * Returns:
 * - true when induced family-id maps are equal.
 */
function _ps_families_equal(keys_a, keys_b) =
    let(
        fam_a = _ps_group_by_key(keys_a),
        fam_b = _ps_group_by_key(keys_b),
        ids_a = _ps_family_ids(len(keys_a), fam_a),
        ids_b = _ps_family_ids(len(keys_b), fam_b)
    )
    _ps_keys_equal(ids_a, ids_b);

/**
 * Iterative face-key refinement until stable families or budget exhausted.
 *
 * Args:
 * - `poly`: descriptor.
 * - `keys`: initial face keys.
 * - `vert_ids`: current vertex family ids.
 * - `edges`, `edge_faces`, `edge_keys`, `nv`: adjacency/lookup support.
 * - `max_iter`: iteration budget.
 *
 * Returns:
 * - refined face keys.
 */
function _ps_refine_face_keys_iter(poly, keys, vert_ids, edges, edge_faces, edge_keys, nv, max_iter=6) =
    (max_iter <= 0) ? keys :
    let(next = _ps_refine_face_keys(poly, keys, keys, vert_ids, edges, edge_faces, edge_keys, nv))
    _ps_families_equal(next, keys) ? keys : _ps_refine_face_keys_iter(poly, next, vert_ids, edges, edge_faces, edge_keys, nv, max_iter - 1);

/**
 * Iterative edge-key refinement until stable keys or budget exhausted.
 *
 * Args:
 * - `poly`: descriptor.
 * - `keys`: initial edge keys.
 * - `face_keys`, `vert_keys`: current neighbor keys.
 * - `edges`, `edge_faces`: adjacency.
 * - `max_iter`: iteration budget.
 *
 * Returns:
 * - refined edge keys.
 */
function _ps_refine_edge_keys_iter(poly, keys, face_keys, vert_keys, edges, edge_faces, max_iter=6) =
    (max_iter <= 0) ? keys :
    let(next = _ps_refine_edge_keys(poly, keys, face_keys, vert_keys, edges, edge_faces))
    _ps_keys_equal(next, keys) ? keys : _ps_refine_edge_keys_iter(poly, next, face_keys, vert_keys, edges, edge_faces, max_iter - 1);

/**
 * Iterative vertex-key refinement until stable families or budget exhausted.
 *
 * Args:
 * - `poly`: descriptor.
 * - `keys`: initial vertex keys.
 * - `face_keys`: current face keys.
 * - `edges`, `edge_faces`: adjacency.
 * - `max_iter`: iteration budget.
 *
 * Returns:
 * - refined vertex keys.
 */
function _ps_refine_vert_keys_iter(poly, keys, face_keys, edges, edge_faces, max_iter=6) =
    (max_iter <= 0) ? keys :
    let(next = _ps_refine_vert_keys(poly, keys, face_keys, edges, edge_faces))
    _ps_families_equal(next, keys) ? keys : _ps_refine_vert_keys_iter(poly, next, face_keys, edges, edge_faces, max_iter - 1);

/**
 * Classify polyhedral faces/edges/vertices into structural families.
 *
 * Approach:
 * - build topology base keys for each element type,
 * - optionally refine keys using neighboring-family context,
 * - optionally append geometric signatures (rounded by eps),
 * - group equal keys into families.
 *
 * Args:
 * - `poly`: poly descriptor.
 * - `detail`:
 *    - `0`: topology only.
 *    - `1`: one-pass neighbor refinement.
 *    - `2+`: iterative neighbor refinement (budget controlled by `radius`).
 * - `eps`: geometric rounding epsilon for include_geom signatures.
 * - `radius`: iterative refinement depth budget for `detail>=2`.
 * - `include_geom`: append geometric signatures to keys.
 *
 * Returns:
 * - classification tuple `[face_families, edge_families, vert_families]`
 *   where each families list is `[[key, idxs], ...]`.
 */
function poly_classify(poly, detail=1, eps=1e-6, radius=1, include_geom=false) =
    let(
        radius_eff = is_undef(radius) ? 1 : radius,
        faces = poly_faces(poly),
        edges = _ps_edges_from_faces(faces),
        edge_faces = ps_edge_faces_table(faces, edges),
        verts = poly_verts(poly),
        edge_keys = _ps_edge_keys_list(edges, len(verts)),

        // Topology-only keys
        face_topo = _ps_face_keys_from(verts, faces, 0, eps),
        edge_topo = _ps_edge_keys_from(verts, faces, edges, edge_faces, 0, eps),
        vert_topo = _ps_vert_keys_from(verts, faces, edges, edge_faces, 0, eps),

        // Optional geometry keys (avg edge lengths)
        face_geom = include_geom ? _ps_face_keys_from(verts, faces, 1, eps) : [],
        edge_geom = include_geom ? _ps_edge_keys_from(verts, faces, edges, edge_faces, 1, eps) : [],
        vert_geom = include_geom ? _ps_vert_keys_from(verts, faces, edges, edge_faces, 1, eps) : [],

        // Neighbor refinement (topology only)
        face_ids_topo = _ps_family_ids(len(faces), _ps_group_by_key(face_topo)),
        vert_ids_topo = _ps_family_ids(len(verts), _ps_group_by_key(vert_topo)),

        face_ref = (detail >= 1)
            ? ((detail >= 2)
                ? _ps_refine_face_keys_iter(poly, face_topo, vert_ids_topo, edges, edge_faces, edge_keys, len(verts), radius_eff)
                : _ps_refine_face_keys(poly, face_topo, face_topo, vert_ids_topo, edges, edge_faces, edge_keys, len(verts)))
            : face_topo,
        vert_ref = (detail >= 1)
            ? ((detail >= 2)
                ? _ps_refine_vert_keys_iter(poly, vert_topo, face_ref, edges, edge_faces, radius_eff)
                : _ps_refine_vert_keys(poly, vert_topo, face_topo, edges, edge_faces))
            : vert_topo,
        edge_ref = (detail >= 1)
            ? _ps_refine_edge_keys(poly, edge_topo, face_ref, vert_ref, edges, edge_faces)
            : edge_topo,

        // Add geometry last (optional)
        face_keys = include_geom
            ? [for (i = [0:1:len(face_ref)-1]) concat(face_ref[i], [face_geom[i][1]])]
            : face_ref,
        edge_keys2 = include_geom
            ? [for (i = [0:1:len(edge_ref)-1]) concat(edge_ref[i], [edge_geom[i][2]])]
            : edge_ref,
        vert_keys = include_geom
            ? [for (i = [0:1:len(vert_ref)-1]) concat(vert_ref[i], [vert_geom[i][len(vert_geom[i])-1]])]
            : vert_ref,
        face_fams = _ps_group_by_key(face_keys),
        edge_fams = _ps_group_by_key(edge_keys2),
        vert_fams = _ps_group_by_key(vert_keys)
    )
    [face_fams, edge_fams, vert_fams];

/**
 * Safe accessor for face families from classification tuple.
 *
 * Args:
 * - `cls`: classification tuple.
 *
 * Returns:
 * - face families list or `[]`.
 */
// Access helpers for classification tuples [face_fams, edge_fams, vert_fams].
function ps_classify_face_families(cls) =
    (is_undef(cls) || len(cls) < 1 || is_undef(cls[0])) ? [] : cls[0];

/**
 * Safe accessor for edge families from classification tuple.
 *
 * Args:
 * - `cls`: classification tuple.
 *
 * Returns:
 * - edge families list or `[]`.
 */
function ps_classify_edge_families(cls) =
    (is_undef(cls) || len(cls) < 2 || is_undef(cls[1])) ? [] : cls[1];

/**
 * Safe accessor for vertex families from classification tuple.
 *
 * Args:
 * - `cls`: classification tuple.
 *
 * Returns:
 * - vertex families list or `[]`.
 */
function ps_classify_vert_families(cls) =
    (is_undef(cls) || len(cls) < 3 || is_undef(cls[2])) ? [] : cls[2];

/**
 * Return family counts for all three element types.
 *
 * Args:
 * - `cls`: classification tuple.
 *
 * Returns:
 * - `[face_family_count, edge_family_count, vert_family_count]`.
 */
function ps_classify_counts(cls) =
    let(
        ff = ps_classify_face_families(cls),
        ef = ps_classify_edge_families(cls),
        vf = ps_classify_vert_families(cls)
    )
    [len(ff), len(ef), len(vf)];

/**
 * Build per-face family-id map from classification.
 *
 * Args:
 * - `cls`: classification tuple.
 * - `n_faces`: expected face count.
 *
 * Returns:
 * - id map length `n_faces`.
 */
function ps_classify_face_ids(cls, n_faces) =
    _ps_family_ids(n_faces, ps_classify_face_families(cls));

/**
 * Build per-edge family-id map from classification.
 *
 * Args:
 * - `cls`: classification tuple.
 * - `n_edges`: expected edge count.
 *
 * Returns:
 * - id map length `n_edges`.
 */
function ps_classify_edge_ids(cls, n_edges) =
    _ps_family_ids(n_edges, ps_classify_edge_families(cls));

/**
 * Build per-vertex family-id map from classification.
 *
 * Args:
 * - `cls`: classification tuple.
 * - `n_verts`: expected vertex count.
 *
 * Returns:
 * - id map length `n_verts`.
 */
function ps_classify_vert_ids(cls, n_verts) =
    _ps_family_ids(n_verts, ps_classify_vert_families(cls));

/**
 * Convenience selector: collect face indices by face arity.
 *
 * Args:
 * - `cls`: classification tuple.
 * - `n`: desired face vertex count (e.g. 3 for triangles).
 *
 * Returns:
 * - concatenated face index list across all matching face families.
 */
function ps_classify_face_idxs_by_n(cls, n) =
    let(ff = ps_classify_face_families(cls))
    [for (fam = ff) if (len(fam[0]) > 0 && fam[0][0] == n) each fam[1]];

/**
 * Debug/inspection printer for classification output.
 *
 * Args:
 * - same as `poly_classify(...)`.
 *
 * Returns:
 * - module side effects only (echo output).
 */
// Pretty-print classification info for a poly.
module show_poly(poly, detail=1, eps=1e-6, radius=1, include_geom=false) {
    cls = poly_classify(poly, detail, eps, radius, include_geom);
    face_fams = cls[0];
    edge_fams = cls[1];
    vert_fams = cls[2];

    echo("=== poly_classify ===");
    echo("face_families:", len(face_fams));
    for (i = [0:1:len(face_fams)-1])
        let(k = face_fams[i][0], idxs = face_fams[i][1])
            echo("  face_family#", i, "key=", k, include_geom ? "(n, avg_edge_len)" : "(n)", "count=", len(idxs), "idxs=", idxs);

    echo("edge_families:", len(edge_fams));
    for (i = [0:1:len(edge_fams)-1])
        let(k = edge_fams[i][0], idxs = edge_fams[i][1])
            echo("  edge_family#", i, "key=", k, include_geom ? "(adj face sizes, edge_len)" : "(adj face sizes)", "count=", len(idxs), "idxs=", idxs);

    echo("vert_families:", len(vert_fams));
    for (i = [0:1:len(vert_fams)-1])
        let(k = vert_fams[i][0], idxs = vert_fams[i][1])
            echo("  vert_family#", i, "key=", k, include_geom ? "(valence, face sizes..., avg_edge_len)" : "(valence, face sizes...)", "count=", len(idxs), "idxs=", idxs);
    echo("====================");
}
