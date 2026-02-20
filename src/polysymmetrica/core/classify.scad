use <funcs.scad>
use <duals.scad>

// Polyhedral element classification by families (faces / edges / vertices).
// detail=0: topology only
// detail=1: adds length-based info (rounded by eps)

function _ps_round(v, eps) =
    (eps <= 0) ? v : (round(v / eps) * eps);

function _ps_edge_len(verts, e) =
    norm(verts[e[0]] - verts[e[1]]);

function _ps_unique_keys(keys) =
    [
        for (i = [0:1:len(keys)-1])
            if (sum([for (j = [0:1:i-1]) keys[j] == keys[i] ? 1 : 0]) == 0)
                keys[i]
    ];

function _ps_key_str(k) = str(k);

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

// Map each element index -> family id.
function _ps_family_ids(n, fams) =
    [
        for (i = [0:1:n-1])
            let(hit = [for (fi = [0:1:len(fams)-1]) if (search(i, fams[fi][1]) != []) fi])
            (len(hit) == 0) ? -1 : hit[0]
    ];

function _ps_rotate(list, k) =
    let(n = len(list))
    [for (i = [0:1:n-1]) list[(i + k) % n]];

function _ps_lex_gt(a, b) =
    (len(a) == 0) ? false :
    (a[0] > b[0]) ? true :
    (a[0] < b[0]) ? false :
    _ps_lex_gt([for (i = [1:1:len(a)-1]) a[i]], [for (i = [1:1:len(b)-1]) b[i]]);

function _ps_lex_lt(a, b) =
    (len(a) == 0) ? false :
    (a[0] < b[0]) ? true :
    (a[0] > b[0]) ? false :
    _ps_lex_lt([for (i = [1:1:len(a)-1]) a[i]], [for (i = [1:1:len(b)-1]) b[i]]);

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

function _ps_edge_keys_list(edges, nv) =
    [for (e = edges) min(e[0], e[1]) * nv + max(e[0], e[1])];

function _ps_edge_index(edge_keys, a, b, nv) =
    let(e = min(a, b) * nv + max(a, b), idxs = search(e, edge_keys))
    idxs[0];

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

function _ps_keys_equal(a, b) =
    (len(a) == len(b)) ? (min([for (i = [0:1:len(a)-1]) a[i] == b[i] ? 1 : 0]) == 1) : false;

function _ps_families_equal(keys_a, keys_b) =
    let(
        fam_a = _ps_group_by_key(keys_a),
        fam_b = _ps_group_by_key(keys_b),
        ids_a = _ps_family_ids(len(keys_a), fam_a),
        ids_b = _ps_family_ids(len(keys_b), fam_b)
    )
    _ps_keys_equal(ids_a, ids_b);

function _ps_refine_face_keys_iter(poly, keys, vert_ids, edges, edge_faces, edge_keys, nv, max_iter=6) =
    (max_iter <= 0) ? keys :
    let(next = _ps_refine_face_keys(poly, keys, keys, vert_ids, edges, edge_faces, edge_keys, nv))
    _ps_families_equal(next, keys) ? keys : _ps_refine_face_keys_iter(poly, next, vert_ids, edges, edge_faces, edge_keys, nv, max_iter - 1);

function _ps_refine_edge_keys_iter(poly, keys, face_keys, vert_keys, edges, edge_faces, max_iter=6) =
    (max_iter <= 0) ? keys :
    let(next = _ps_refine_edge_keys(poly, keys, face_keys, vert_keys, edges, edge_faces))
    _ps_keys_equal(next, keys) ? keys : _ps_refine_edge_keys_iter(poly, next, face_keys, vert_keys, edges, edge_faces, max_iter - 1);

function _ps_refine_vert_keys_iter(poly, keys, face_keys, edges, edge_faces, max_iter=6) =
    (max_iter <= 0) ? keys :
    let(next = _ps_refine_vert_keys(poly, keys, face_keys, edges, edge_faces))
    _ps_families_equal(next, keys) ? keys : _ps_refine_vert_keys_iter(poly, next, face_keys, edges, edge_faces, max_iter - 1);

// Return [face_families, edge_families, vert_families].
// Each family is [key, idxs].
// detail:
//   0 = topology only
//   1 = topology + neighbour refinement
//   2 = iterated neighbour refinement (radius)
// include_geom: when true, append geometry info (avg edge lengths).
// radius controls how far neighbour refinement propagates (default 1).
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
