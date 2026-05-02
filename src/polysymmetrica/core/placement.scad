// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
include <segments.scad>
use <classify.scad>

function _ps_cls_opt(classify_opts, i, def) =
    (is_undef(classify_opts) || !is_list(classify_opts) || i >= len(classify_opts) || is_undef(classify_opts[i]))
        ? def
        : classify_opts[i];

/**
 * Function: Test whether a placement element index is selected.
 * Params: idx (element index), indices (`undef`, scalar index, or list of indices)
 * Returns: true when the element should be visited
 */
function _ps_place_idx_selected(idx, indices) =
    is_undef(indices)
        ? true
        : is_list(indices)
            ? _ps_list_contains(indices, idx)
            : idx == indices;

function _ps_resolve_classify(poly, classify=undef, classify_opts=undef) =
    is_undef(classify) && is_undef(classify_opts)
        ? undef
        : !is_undef(classify)
        ? classify
        : poly_classify(
            poly,
            _ps_cls_opt(classify_opts, 0, 1),
            _ps_cls_opt(classify_opts, 1, 1e-6),
            _ps_cls_opt(classify_opts, 2, 1),
            _ps_cls_opt(classify_opts, 3, false)
        );

/**
 * Function: Build face-neighbor indices in face-edge order for one face site.
 * Params: face (face index cycle), fi (face index), faces0/edges/edge_faces (oriented topology tables)
 * Returns: list of adjacent face indices, or undef on boundary edges
 */
function _ps_face_site_neighbors_idx(face, fi, faces0, edges, edge_faces) =
    [
        for (k = [0:1:len(face)-1])
            let(
                v0 = face[k],
                v1 = face[(k+1)%len(face)],
                ei = ps_find_edge_index(edges, v0, v1),
                adj = edge_faces[ei]
            )
            (len(adj) < 2) ? undef : ((adj[0] == fi) ? adj[1] : adj[0])
    ];

/**
 * Function: Build per-edge face dihedrals in face-edge order for one face site.
 * Params: face (face index cycle), fi (face index), faces0/edges/edge_faces (oriented topology tables), face_n (per-face normals)
 * Returns: list of dihedral angles in degrees, aligned with the face edge order
 */
function _ps_face_site_dihedrals(face, fi, faces0, edges, edge_faces, face_n) =
    [
        for (k = [0:1:len(face)-1])
            let(
                v0 = face[k],
                v1 = face[(k+1)%len(face)],
                ei = ps_find_edge_index(edges, v0, v1),
                adj = edge_faces[ei],
                n0 = face_n[fi],
                n1 = (len(adj) < 2) ? n0 : face_n[(adj[0] == fi) ? adj[1] : adj[0]],
                dotn = v_dot(n0, n1),
                c = (dotn > 1) ? 1 : ((dotn < -1) ? -1 : dotn)
            )
            180 - acos(c)
    ];

/**
 * Function: Build full neighbor indices for one vertex site from the edge list.
 * Params: edges (undirected edge list), vi (vertex index)
 * Returns: list of neighboring vertex indices in edge scan order
 */
function _ps_vertex_site_neighbors_idx(edges, vi) =
    [
        for (e = edges)
            if (e[0] == vi) e[1]
            else if (e[1] == vi) e[0]
    ];

/**
 * Function: Pick a stable axis perpendicular to a normal for degenerate local face frames.
 * Params: n (unit normal)
 * Returns: unit vector perpendicular to `n`
 */
function _ps_any_perp(n) =
    abs(n[0]) < 0.9
        ? v_norm(v_cross(n, [1, 0, 0]))
        : v_norm(v_cross(n, [0, 1, 0]));

/**
 * Function: Build a projected local face-frame X axis from target-local poly vertices.
 * Params: verts_local (poly vertices in target-local coordinates), f (face index loop), center (face center), ez (face normal), eps (degeneracy tolerance)
 * Returns: unit X axis in target-local coordinates
 */
function _ps_local_face_ex(verts_local, f, center, ez, eps=1e-12) =
    let(
        ex_raw = verts_local[f[0]] - center,
        ex_proj = ex_raw - ez * v_dot(ex_raw, ez),
        ex_fallback_raw = verts_local[f[1]] - verts_local[f[0]],
        ex_fallback = ex_fallback_raw - ez * v_dot(ex_fallback_raw, ez)
    )
    norm(ex_proj) > eps
        ? v_norm(ex_proj)
        : norm(ex_fallback) > eps
        ? v_norm(ex_fallback)
        : _ps_any_perp(ez);

/**
 * Function: Build a canonical face placement site from vertices already in a local coordinate system.
 * Params: face_idx (face index), faces (face list), verts_local (poly vertices in parent-local coordinates), poly_center_local_parent (optional poly center in parent-local coords), eps (degeneracy tolerance)
 * Returns: face site record matching `ps_face_sites(...)`
 */
function _ps_face_site_from_local_poly(face_idx, faces, verts_local, poly_center_local_parent=undef, eps=1e-12) =
    let(
        f = faces[face_idx],
        center = ps_face_centroid(verts_local, f),
        ez = ps_face_frame_normal(verts_local, f, eps),
        ex = _ps_local_face_ex(verts_local, f, center, ez, eps),
        ey = v_cross(ez, ex),
        edge_lens = [
            for (k = [0:1:len(f)-1])
                norm(verts_local[f[(k + 1) % len(f)]] - verts_local[f[k]])
        ],
        edge_len = len(edge_lens) == 0 ? 0 : sum(edge_lens) / len(edge_lens),
        face_midradius = norm(center),
        rad_vec = [for (vid = f) norm(verts_local[vid] - center)],
        face_radius = len(rad_vec) == 0 ? 0 : sum(rad_vec) / len(rad_vec),
        poly_center_parent = is_undef(poly_center_local_parent) ? [0, 0, 0] : poly_center_local_parent,
        poly_center_delta = poly_center_parent - center,
        poly_center_local_raw = [
            v_dot(poly_center_delta, ex),
            v_dot(poly_center_delta, ey),
            v_dot(poly_center_delta, ez)
        ],
        face_verts_local_raw = [
            for (vid = f)
                let(p = verts_local[vid] - center)
                    [v_dot(p, ex), v_dot(p, ey), v_dot(p, ez)]
        ],
        poly_verts_local_raw = [
            for (p0 = verts_local)
                let(p = p0 - center)
                    [v_dot(p, ex), v_dot(p, ey), v_dot(p, ez)]
        ],
        zvals = [for (p = face_verts_local_raw) p[2]],
        zmean = len(zvals) == 0 ? 0 : sum(zvals) / len(zvals),
        face_planarity_err = len(zvals) == 0 ? 0 : max([for (z = zvals) abs(z - zmean)]),
        face_pts3d_local = [for (p = face_verts_local_raw) [p[0], p[1], p[2] - zmean]],
        poly_center_local = [poly_center_local_raw[0], poly_center_local_raw[1], poly_center_local_raw[2] - zmean],
        poly_verts_local = [for (p = poly_verts_local_raw) [p[0], p[1], p[2] - zmean]],
        face_pts2d = ps_xy(face_pts3d_local),
        edges = _ps_edges_from_faces(faces),
        edge_faces = ps_edge_faces_table(faces, edges),
        face_n = [for (face = faces) ps_face_normal(verts_local, face)],
        face_neighbors_idx = _ps_face_site_neighbors_idx(f, face_idx, faces, edges, edge_faces),
        face_dihedrals = _ps_face_site_dihedrals(f, face_idx, faces, edges, edge_faces, face_n)
    )
    [
        face_idx,
        center,
        ex,
        ey,
        ez,
        edge_len,
        len(face_pts2d),
        face_midradius,
        face_radius,
        poly_center_local,
        face_pts2d,
        face_pts3d_local,
        poly_verts_local,
        faces,
        face_planarity_err,
        face_planarity_err <= eps,
        undef,
        undef,
        undef,
        undef,
        face_neighbors_idx,
        face_dihedrals
    ];

/**
 * Function: Find the canonical global edge index for two vertex ids.
 * Params: edges (canonical edge list), a/b (edge endpoint vertex ids)
 * Returns: edge index, or `-1` when absent
 */
function _ps_edge_idx_from_verts(edges, a, b) =
    _ps_index_of(edges, (a < b) ? [a, b] : [b, a]);

/**
 * Function: Test whether a 3D segment intersects the local Z=0 plane.
 * Params: a/b (3D endpoints), eps (tolerance)
 * Returns: boolean
 */
function _ps_seg3_hits_z0(a, b, eps=1e-8) =
    min(a[2], b[2]) <= eps && max(a[2], b[2]) >= -eps;

/**
 * Function: Deduplicate scalar values while preserving first-seen order.
 * Params: vals (input values), i/acc (recursion state)
 * Returns: unique values
 */
function _ps_unique_values(vals, i=0, acc=[]) =
    (i >= len(vals)) ? acc :
    _ps_unique_values(
        vals,
        i + 1,
        _ps_list_contains(acc, vals[i]) ? acc : concat(acc, [vals[i]])
    );

/**
 * Function: Build a canonical edge placement site from vertices already in a local coordinate system.
 * Params: edge_idx (global edge index), faces (face list), verts_local (poly vertices in parent-local coordinates), poly_center_local_parent (optional poly center in parent-local coords), eps (degeneracy tolerance)
 * Returns: edge site record matching `ps_edge_sites(...)`
 */
function _ps_edge_site_from_local_poly(edge_idx, faces, verts_local, poly_center_local_parent=undef, eps=1e-12) =
    let(
        edges = _ps_edges_from_faces(faces),
        e = edges[edge_idx],
        v0 = verts_local[e[0]],
        v1 = verts_local[e[1]],
        center = (v0 + v1) / 2,
        edge_vec = v1 - v0,
        ex = (norm(edge_vec) <= eps) ? [1, 0, 0] : v_norm(edge_vec),
        edge_faces = ps_edge_faces_table(faces, edges),
        adj_faces_idx = edge_faces[edge_idx],
        face_n = [for (f = faces) ps_face_normal(verts_local, f)],
        poly_center_parent = is_undef(poly_center_local_parent) ? [0, 0, 0] : poly_center_local_parent,
        radial_raw = center - poly_center_parent,
        radial_ref = (norm(radial_raw) <= eps) ? _ps_any_perp(ex) : v_norm(radial_raw),
        bisector_raw =
            (len(adj_faces_idx) < 2)
                ? radial_ref
                : face_n[adj_faces_idx[0]] + face_n[adj_faces_idx[1]],
        bisector_signed =
            (norm(bisector_raw) <= eps)
                ? radial_ref
                : ((v_dot(bisector_raw, radial_ref) < 0) ? -bisector_raw : bisector_raw),
        ez_proj = bisector_signed - ex * v_dot(bisector_signed, ex),
        radial_proj = radial_ref - ex * v_dot(radial_ref, ex),
        ez_dir =
            (norm(ez_proj) <= eps)
                ? radial_proj
                : ez_proj,
        ez = (norm(ez_dir) <= eps) ? _ps_any_perp(ex) : v_norm(ez_dir),
        ey = v_norm(v_cross(ez, ex)),
        edge_midradius = norm(center - poly_center_parent),
        edge_len_actual = norm(edge_vec),
        poly_center_delta = poly_center_parent - center,
        poly_center_local = [
            v_dot(poly_center_delta, ex),
            v_dot(poly_center_delta, ey),
            v_dot(poly_center_delta, ez)
        ],
        edge_pts_local = [[-edge_len_actual / 2, 0, 0], [edge_len_actual / 2, 0, 0]]
    )
    [
        edge_idx,
        center,
        ex,
        ey,
        ez,
        edge_len_actual,
        edge_midradius,
        poly_center_local,
        edge_pts_local,
        e,
        adj_faces_idx,
        undef,
        undef,
        undef,
        undef
    ];

/**
 * Function: Build a canonical vertex placement site from vertices already in a local coordinate system.
 * Params: vertex_idx (vertex index), faces (face list), verts_local (poly vertices in parent-local coordinates), poly_center_local_parent (optional poly center in parent-local coords), eps (degeneracy tolerance)
 * Returns: vertex site record matching `ps_vertex_sites(...)`
 */
function _ps_vertex_site_from_local_poly(vertex_idx, faces, verts_local, poly_center_local_parent=undef, eps=1e-12) =
    let(
        edges = _ps_edges_from_faces(faces),
        center = verts_local[vertex_idx],
        poly_center_parent = is_undef(poly_center_local_parent) ? [0, 0, 0] : poly_center_local_parent,
        radial_raw = center - poly_center_parent,
        ez = (norm(radial_raw) <= eps) ? [0, 0, 1] : v_norm(radial_raw),
        neighbors_idx = _ps_vertex_site_neighbors_idx(edges, vertex_idx),
        neighbor0 = len(neighbors_idx) == 0 ? undef : verts_local[neighbors_idx[0]],
        neighbor_dir = is_undef(neighbor0) ? undef : neighbor0 - center,
        proj = is_undef(neighbor_dir) ? undef : neighbor_dir - ez * v_dot(neighbor_dir, ez),
        ex =
            (!is_undef(proj) && norm(proj) > eps)
                ? v_norm(proj)
                : _ps_any_perp(ez),
        ey = v_cross(ez, ex),
        valence = len(neighbors_idx),
        neighbor_pts_local = [
            for (nj = neighbors_idx)
                let(pw = verts_local[nj] - center)
                    [v_dot(pw, ex), v_dot(pw, ey), v_dot(pw, ez)]
        ],
        neighbor_lens = [for (p = neighbor_pts_local) norm(p)],
        edge_len = len(neighbor_lens) == 0 ? 0 : sum(neighbor_lens) / len(neighbor_lens),
        poly_center_delta = poly_center_parent - center,
        poly_center_local = [
            v_dot(poly_center_delta, ex),
            v_dot(poly_center_delta, ey),
            v_dot(poly_center_delta, ez)
        ]
    )
    [
        vertex_idx,
        center,
        ex,
        ey,
        ez,
        edge_len,
        norm(radial_raw),
        poly_center_local,
        valence,
        neighbors_idx,
        neighbor_pts_local,
        undef,
        undef,
        undef,
        undef
    ];

/**
 * Function: Build edge ids from one foreign face intrusion record.
 * Params: record (foreign face intrusion record), faces (face list), verts_local (target face-local poly vertices), eps (tolerance)
 * Returns: global edge ids for boundary edges of the foreign face that intersect local Z=0
 */
function _ps_proxy_edge_ids_from_face_record(record, faces, verts_local, eps=1e-8) =
    let(
        face_idx = ps_intrusion_foreign_idx(record),
        f = faces[face_idx],
        edges = _ps_edges_from_faces(faces),
        ids = [
            for (k = [0:1:len(f)-1])
                let(
                    a = f[k],
                    b = f[(k + 1) % len(f)],
                    edge_idx = _ps_edge_idx_from_verts(edges, a, b)
                )
                if (edge_idx >= 0 && _ps_seg3_hits_z0(verts_local[a], verts_local[b], eps))
                    edge_idx
        ]
    )
    _ps_unique_values(ids);

/**
 * Function: Build vertex ids from one foreign face intrusion record.
 * Params: record (foreign face intrusion record), faces (face list), verts_local (target face-local poly vertices), eps (tolerance)
 * Returns: vertex ids at the ends of boundary edges of the foreign face that intersect local Z=0
 */
function _ps_proxy_vertex_ids_from_face_record(record, faces, verts_local, eps=1e-8) =
    let(
        face_idx = ps_intrusion_foreign_idx(record),
        f = faces[face_idx],
        raw = [
            for (k = [0:1:len(f)-1])
                let(
                    a = f[k],
                    b = f[(k + 1) % len(f)]
                )
                if (_ps_seg3_hits_z0(verts_local[a], verts_local[b], eps))
                    each [a, b]
        ]
    )
    _ps_unique_values(raw);

/**
 * Function: Test whether a candidate record set already contains a foreign element.
 * Params: records (candidate records), kind (foreign kind), idx (foreign index), i (recursion state)
 * Returns: boolean
 */
function _ps_proxy_candidate_record_seen(records, kind, idx, i=0) =
    (i >= len(records)) ? false :
    (ps_intrusion_foreign_kind(records[i]) == kind && ps_intrusion_foreign_idx(records[i]) == idx) ? true :
    _ps_proxy_candidate_record_seen(records, kind, idx, i + 1);

/**
 * Function: Deduplicate proxy candidate records by foreign kind and source index.
 * Params: records (candidate records), i/acc (recursion state)
 * Returns: first record for each `(foreign_kind, foreign_idx)` pair
 */
function _ps_proxy_candidate_records_dedupe(records, i=0, acc=[]) =
    (i >= len(records)) ? acc :
    let(
        record = records[i],
        hit = _ps_proxy_candidate_record_seen(
            acc,
            ps_intrusion_foreign_kind(record),
            ps_intrusion_foreign_idx(record)
        )
    )
    _ps_proxy_candidate_records_dedupe(
        records,
        i + 1,
        hit ? acc : concat(acc, [record])
    );

/**
 * Function: Build one foreign face replay site from an intrusion record.
 * Params: replay_idx (site index), record (foreign intrusion record), poly_faces_idx/poly_verts_local/poly_center_local (target-local poly context), eps (degeneracy tolerance)
 * Returns: replay site record for `place_on_face_foreign_face_replay_sites(...)`
 */
function _ps_face_foreign_face_replay_site(replay_idx, record, poly_faces_idx, poly_verts_local, poly_center_local, eps=1e-12) =
    let(
        foreign_face_idx = ps_intrusion_foreign_idx(record),
        face_site = _ps_face_site_from_local_poly(foreign_face_idx, poly_faces_idx, poly_verts_local, poly_center_local, eps)
    )
    [
        replay_idx,
        record,
        face_site[1],
        face_site[2],
        face_site[3],
        face_site[4],
        foreign_face_idx,
        face_site[10],
        face_site[11],
        face_site[12],
        face_site[9],
        face_site[13][foreign_face_idx],
        ps_intrusion_foreign_kind(record),
        ps_intrusion_segment2d_local(record),
        ps_intrusion_dihedral(record),
        ps_intrusion_confidence(record),
        face_site,
        undef,
        undef
    ];

/**
 * Function: Build one foreign edge replay site from a candidate intrusion record.
 * Params: replay_idx (site index), record (foreign edge candidate record), poly_faces_idx/poly_verts_local/poly_center_local (target-local poly context), eps (degeneracy tolerance)
 * Returns: replay site record for `place_on_face_foreign_proxy_sites(...)`
 */
function _ps_face_foreign_edge_replay_site(replay_idx, record, poly_faces_idx, poly_verts_local, poly_center_local, eps=1e-12) =
    let(
        foreign_edge_idx = ps_intrusion_foreign_idx(record),
        edge_site = _ps_edge_site_from_local_poly(foreign_edge_idx, poly_faces_idx, poly_verts_local, poly_center_local, eps)
    )
    [
        replay_idx,
        record,
        edge_site[1],
        edge_site[2],
        edge_site[3],
        edge_site[4],
        foreign_edge_idx,
        undef,
        undef,
        undef,
        edge_site[7],
        undef,
        ps_intrusion_foreign_kind(record),
        ps_intrusion_segment2d_local(record),
        ps_intrusion_dihedral(record),
        ps_intrusion_confidence(record),
        undef,
        edge_site,
        undef
    ];

/**
 * Function: Build one foreign vertex replay site from a candidate intrusion record.
 * Params: replay_idx (site index), record (foreign vertex candidate record), poly_faces_idx/poly_verts_local/poly_center_local (target-local poly context), eps (degeneracy tolerance)
 * Returns: replay site record for `place_on_face_foreign_proxy_sites(...)`
 */
function _ps_face_foreign_vertex_replay_site(replay_idx, record, poly_faces_idx, poly_verts_local, poly_center_local, eps=1e-12) =
    let(
        foreign_vertex_idx = ps_intrusion_foreign_idx(record),
        vertex_site = _ps_vertex_site_from_local_poly(foreign_vertex_idx, poly_faces_idx, poly_verts_local, poly_center_local, eps)
    )
    [
        replay_idx,
        record,
        vertex_site[1],
        vertex_site[2],
        vertex_site[3],
        vertex_site[4],
        foreign_vertex_idx,
        undef,
        undef,
        undef,
        vertex_site[7],
        undef,
        ps_intrusion_foreign_kind(record),
        ps_intrusion_segment2d_local(record),
        ps_intrusion_dihedral(record),
        ps_intrusion_confidence(record),
        undef,
        undef,
        vertex_site
    ];

/**
 * Function: Build target-local replay sites for exact foreign face intrusions.
 * Params: face_pts2d (target face loop), face_idx (target face index), poly_faces_idx/poly_verts_local/poly_center_local (current `place_on_faces(...)` metadata), eps (tolerance), mode (foreign face fill rule), filter_parent (drop parent-edge cuts)
 * Returns: replay site records for intruding foreign faces, with frames expressed in the target face-local coordinate system
 * Limitations/Gotchas: only exact `"face"` intrusion records are converted; use `ps_face_foreign_proxy_replay_sites(...)` for edge/vertex candidates
 */
function ps_face_foreign_face_replay_sites(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, poly_center_local=undef, eps=1e-8, mode="nonzero", filter_parent=true) =
    let(
        records = [
            for (r = ps_face_foreign_intrusion_records(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps, mode, filter_parent))
                if (ps_intrusion_foreign_kind(r) == "face")
                    r
        ]
    )
    [
        for (ri = [0:1:len(records)-1])
            _ps_face_foreign_face_replay_site(ri, records[ri], poly_faces_idx, poly_verts_local, poly_center_local, eps)
    ];

/**
 * Function: Build provenance-driven proxy replay sites for foreign face/edge/vertex sources.
 * Params: face_pts2d (target face loop), face_idx (target face index), poly_faces_idx/poly_verts_local/poly_center_local (current `place_on_faces(...)` metadata), eps (tolerance), mode (foreign face fill rule), filter_parent (drop parent-edge cuts)
 * Returns: replay site records for exact foreign faces plus candidate foreign edges/vertices implicated by those face-plane cuts
 * Limitations/Gotchas: edge and vertex sites are candidate/provenance records, not distance-envelope proximity tests
 */
function ps_face_foreign_proxy_replay_sites(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, poly_center_local=undef, eps=1e-8, mode="nonzero", filter_parent=true) =
    let(
        face_records = [
            for (r = ps_face_foreign_intrusion_records(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps, mode, filter_parent))
                if (ps_intrusion_foreign_kind(r) == "face")
                    r
        ],
        edge_records = [
            for (r = face_records)
                for (edge_idx = _ps_proxy_edge_ids_from_face_record(r, poly_faces_idx, poly_verts_local, eps))
                    ["face_plane_cut_candidate", face_idx, "edge", edge_idx, ps_intrusion_segment2d_local(r), ps_intrusion_dihedral(r), "candidate"]
        ],
        vertex_records = [
            for (r = face_records)
                for (vertex_idx = _ps_proxy_vertex_ids_from_face_record(r, poly_faces_idx, poly_verts_local, eps))
                    ["face_plane_cut_candidate", face_idx, "vertex", vertex_idx, ps_intrusion_segment2d_local(r), ps_intrusion_dihedral(r), "candidate"]
        ],
        records = _ps_proxy_candidate_records_dedupe(concat(face_records, edge_records, vertex_records))
    )
    [
        for (ri = [0:1:len(records)-1])
            let(kind = ps_intrusion_foreign_kind(records[ri]))
                kind == "face"
                    ? _ps_face_foreign_face_replay_site(ri, records[ri], poly_faces_idx, poly_verts_local, poly_center_local, eps)
                    : kind == "edge"
                    ? _ps_face_foreign_edge_replay_site(ri, records[ri], poly_faces_idx, poly_verts_local, poly_center_local, eps)
                    : _ps_face_foreign_vertex_replay_site(ri, records[ri], poly_faces_idx, poly_verts_local, poly_center_local, eps)
    ];

/**
 * Function: Get replay site index.
 * Params: site (foreign replay site)
 * Returns: zero-based replay site index
 */
function ps_replay_site_idx(site) = site[0];

/**
 * Function: Get source intrusion record from a replay site.
 * Params: site (foreign replay site)
 * Returns: foreign intrusion record
 */
function ps_replay_site_intrusion_record(site) = site[1];

/**
 * Function: Get target-local replay frame center.
 * Params: site (foreign replay site)
 * Returns: center in current target face-local coordinates
 */
function ps_replay_site_center_local(site) = site[2];

/**
 * Function: Get target-local replay frame X axis.
 * Params: site (foreign replay site)
 * Returns: unit X axis in current target face-local coordinates
 */
function ps_replay_site_ex_local(site) = site[3];

/**
 * Function: Get target-local replay frame Y axis.
 * Params: site (foreign replay site)
 * Returns: unit Y axis in current target face-local coordinates
 */
function ps_replay_site_ey_local(site) = site[4];

/**
 * Function: Get target-local replay frame Z axis.
 * Params: site (foreign replay site)
 * Returns: unit Z axis in current target face-local coordinates
 */
function ps_replay_site_ez_local(site) = site[5];

/**
 * Function: Get foreign element index from a replay site.
 * Params: site (foreign replay site)
 * Returns: foreign face/edge/vertex index, depending on `ps_replay_site_foreign_kind(site)`
 */
function ps_replay_site_foreign_idx(site) = site[6];

/**
 * Function: Get foreign face 2D points in replay local coordinates.
 * Params: site (foreign replay site)
 * Returns: `pts2d` for face sites in replay local coordinates, or `undef`
 */
function ps_replay_site_face_pts2d(site) = site[7];

/**
 * Function: Get foreign face 3D points in replay local coordinates.
 * Params: site (foreign replay site)
 * Returns: `pts3d` for face sites in replay local coordinates, or `undef`
 */
function ps_replay_site_face_pts3d_local(site) = site[8];

/**
 * Function: Get all poly vertices in replay local coordinates.
 * Params: site (foreign replay site)
 * Returns: vertex list transformed into the replay frame
 */
function ps_replay_site_poly_verts_local(site) = site[9];

/**
 * Function: Get poly center vector in replay local coordinates.
 * Params: site (foreign replay site)
 * Returns: vector from replay origin to poly center in replay local coordinates
 */
function ps_replay_site_poly_center_local(site) = site[10];

/**
 * Function: Get foreign face vertex indices from a replay site.
 * Params: site (foreign replay site)
 * Returns: face vertex index loop for face sites, or `undef`
 */
function ps_replay_site_face_verts_idx(site) = site[11];

/**
 * Function: Get foreign element kind from a replay site.
 * Params: site (foreign replay site)
 * Returns: foreign element kind (`"face"`, `"edge"`, or `"vertex"`)
 */
function ps_replay_site_foreign_kind(site) = site[12];

/**
 * Function: Get target-local intrusion segment from a replay site.
 * Params: site (foreign replay site)
 * Returns: `seg2d` in target face-local coordinates
 */
function ps_replay_site_intrusion_segment2d_local(site) = site[13];

/**
 * Function: Get cut dihedral from a replay site.
 * Params: site (foreign replay site)
 * Returns: face-plane cut dihedral
 */
function ps_replay_site_intrusion_dihedral(site) = site[14];

/**
 * Function: Get confidence/classification from a replay site.
 * Params: site (foreign replay site)
 * Returns: confidence string
 */
function ps_replay_site_intrusion_confidence(site) = site[15];

/**
 * Function: Get canonical face placement site from a replay site.
 * Params: site (foreign replay site)
 * Returns: face site record matching `ps_face_sites(...)`, or `undef`
 */
function ps_replay_site_face_site(site) = site[16];

/**
 * Function: Get canonical edge placement site from a replay site.
 * Params: site (foreign replay site)
 * Returns: edge site record matching `ps_edge_sites(...)`, or `undef`
 */
function ps_replay_site_edge_site(site) = site[17];

/**
 * Function: Get canonical vertex placement site from a replay site.
 * Params: site (foreign replay site)
 * Returns: vertex site record matching `ps_vertex_sites(...)`, or `undef`
 */
function ps_replay_site_vertex_site(site) = site[18];

/**
 * Function: Build face placement site records for `place_on_faces(...)`.
 * Params: poly (poly descriptor), inter_radius (scale input), edge_len (explicit scale override), classify/classify_opts (optional classification context)
 * Returns: list of face site records `[face_idx, center, ex, ey, ez, edge_len, vertex_count, face_midradius, face_radius, poly_center_local, face_pts2d, face_pts3d_local, poly_verts_local, poly_faces_idx, face_planarity_err, face_is_planar, face_family_id, face_family_count, edge_family_count, vertex_family_count, face_neighbors_idx, face_dihedrals]`
 * Limitations: record shape is currently positional; keep the semantics stable even if the internal representation changes later
 */
function ps_face_sites(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef) =
    let(
        exp_edge_len = is_undef(edge_len) ? inter_radius * poly_e_over_ir(poly) : edge_len,
        scale = exp_edge_len,
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        faces0 = ps_orient_all_faces_outward(verts, faces),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [ for (f = faces0) ps_face_normal(verts, f) ],
        cls = _ps_resolve_classify(poly, classify, classify_opts),
        family_counts = is_undef(cls) ? undef : ps_classify_counts(cls),
        face_family_ids = is_undef(cls) ? [] : ps_classify_face_ids(cls, len(faces)),
        edge_family_count = is_undef(family_counts) ? undef : family_counts[1],
        vert_family_count = is_undef(family_counts) ? undef : family_counts[2]
    )
    [
        for (fi = [0:1:len(faces)-1])
            let(
                f = faces[fi],
                center = poly_face_center(poly, fi, scale),
                ex = poly_face_ex(poly, fi, scale),
                ey = poly_face_ey(poly, fi, scale),
                ez = poly_face_ez(poly, fi, scale),
                face_midradius = norm(center),
                rad_vec = [for (vid = f) norm(verts[vid] * scale - center)],
                face_radius = sum(rad_vec) / len(rad_vec),
                poly_center_local_raw = [
                    v_dot(-center, ex),
                    v_dot(-center, ey),
                    v_dot(-center, ez)
                ],
                face_verts_local = [
                    for (vid = f)
                        let(p = verts[vid] * scale - center)
                            [v_dot(p, ex), v_dot(p, ey), v_dot(p, ez)]
                ],
                poly_verts_local_raw = [
                    for (vi = [0:1:len(verts)-1])
                        let(p = verts[vi] * scale - center)
                            [v_dot(p, ex), v_dot(p, ey), v_dot(p, ez)]
                ],
                zvals = [for (p = face_verts_local) p[2]],
                zmean = (len(zvals) == 0) ? 0 : sum(zvals) / len(zvals),
                face_planarity_err = (len(zvals) == 0) ? 0 : max([for (z = zvals) abs(z - zmean)]),
                face_pts3d_local = [for (p = face_verts_local) [p[0], p[1], p[2] - zmean]],
                poly_center_local = [poly_center_local_raw[0], poly_center_local_raw[1], poly_center_local_raw[2] - zmean],
                poly_verts_local = [for (p = poly_verts_local_raw) [p[0], p[1], p[2] - zmean]],
                face_pts2d = ps_xy(face_pts3d_local),
                face_neighbors_idx = _ps_face_site_neighbors_idx(f, fi, faces0, edges, edge_faces),
                face_dihedrals = _ps_face_site_dihedrals(f, fi, faces0, edges, edge_faces, face_n)
            )
            [
                fi,
                center,
                ex,
                ey,
                ez,
                exp_edge_len,
                len(face_pts2d),
                face_midradius,
                face_radius,
                poly_center_local,
                face_pts2d,
                face_pts3d_local,
                poly_verts_local,
                faces,
                face_planarity_err,
                face_planarity_err <= 1e-8,
                is_undef(cls) ? undef : face_family_ids[fi],
                is_undef(family_counts) ? undef : family_counts[0],
                edge_family_count,
                vert_family_count,
                face_neighbors_idx,
                face_dihedrals
            ]
    ];

/**
 * Module: Place children on selected faces of a polyhedron.
 * Params: poly (poly descriptor), inter_radius (scale input), edge_len (explicit scale override), classify/classify_opts (optional classification context), indices (`undef`, scalar face index, or list of face indices)
 * Returns: none; exposes `$ps_face_*` metadata for each selected face
 * Limitations: `indices` filters the placement loop only; `ps_face_sites(...)` still builds the complete site list so element ids and classification metadata remain global
 */
module place_on_faces(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef, indices = undef) {
    sites = ps_face_sites(poly, inter_radius, edge_len, classify, classify_opts);

    for (site = sites) {
        fi = site[0];
        if (_ps_place_idx_selected(fi, indices)) {
            center = site[1];
            ex = site[2];
            ey = site[3];
            ez = site[4];

            // Per-face metadata (local-space friendly) - mean average values where faces are irregular
            $ps_face_idx           = fi;
            $ps_edge_len           = site[5];
            $ps_vertex_count       = site[6];
            $ps_face_midradius     = site[7];
            $ps_face_radius        = site[8];
            $ps_poly_center_local  = site[9];
            $ps_face_pts2d         = site[10];
            $ps_face_pts3d_local   = site[11];
            $ps_poly_verts_local   = site[12];
            $ps_poly_faces_idx     = site[13];
            $ps_face_planarity_err = site[14];
            $ps_face_is_planar     = site[15];
            $ps_face_family_id     = site[16];
            $ps_face_family_count  = site[17];
            $ps_edge_family_count  = site[18];
            $ps_vertex_family_count = site[19];
            $ps_face_neighbors_idx = site[20];
            $ps_face_dihedrals     = site[21];

            multmatrix(ps_frame_matrix(center, ex, ey, ez))
                children();
        }
    }
}

/**
 * Module: Replay exact foreign face intrusion sites inside the current placed target face.
 * Params: mode (foreign face fill rule), eps (tolerance), filter_parent (drop parent-edge cuts), coords (`"element"` or `"parent"`)
 * Returns: none; exposes `$ps_replay_*` metadata and optionally places children in the foreign face replay frame
 * Limitations/Gotchas: requires `place_on_faces(...)`; does not generate or subtract proxy geometry
 */
module place_on_face_foreign_face_replay_sites(mode="nonzero", eps=1e-8, filter_parent=true, coords="element") {
    assert(!is_undef($ps_face_pts2d), "place_on_face_foreign_face_replay_sites: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_idx), "place_on_face_foreign_face_replay_sites: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "place_on_face_foreign_face_replay_sites: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "place_on_face_foreign_face_replay_sites: requires place_on_faces context ($ps_poly_verts_local)");
    assert(coords == "element" || coords == "parent", "place_on_face_foreign_face_replay_sites: coords must be \"element\" or \"parent\"");

    sites = ps_face_foreign_face_replay_sites($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, $ps_poly_center_local, eps, mode, filter_parent);
    for (site = sites) {
        face_site = ps_replay_site_face_site(site);
        $ps_replay_idx = ps_replay_site_idx(site);
        $ps_replay_count = len(sites);
        $ps_replay_intrusion_record = ps_replay_site_intrusion_record(site);
        $ps_replay_kind = "foreign_face";
        $ps_replay_foreign_kind = ps_replay_site_foreign_kind(site);
        $ps_replay_foreign_idx = ps_replay_site_foreign_idx(site);
        $ps_replay_center_local = ps_replay_site_center_local(site);
        $ps_replay_ex_local = ps_replay_site_ex_local(site);
        $ps_replay_ey_local = ps_replay_site_ey_local(site);
        $ps_replay_ez_local = ps_replay_site_ez_local(site);
        $ps_replay_face_pts2d = ps_replay_site_face_pts2d(site);
        $ps_replay_face_pts3d_local = ps_replay_site_face_pts3d_local(site);
        $ps_replay_poly_verts_local = ps_replay_site_poly_verts_local(site);
        $ps_replay_poly_center_local = ps_replay_site_poly_center_local(site);
        $ps_replay_face_verts_idx = ps_replay_site_face_verts_idx(site);
        $ps_replay_intrusion_segment2d_local = ps_replay_site_intrusion_segment2d_local(site);
        $ps_replay_intrusion_dihedral = ps_replay_site_intrusion_dihedral(site);
        $ps_replay_intrusion_confidence = ps_replay_site_intrusion_confidence(site);

        if (coords == "element") {
            $ps_face_idx           = face_site[0];
            $ps_edge_len           = face_site[5];
            $ps_vertex_count       = face_site[6];
            $ps_face_midradius     = face_site[7];
            $ps_face_radius        = face_site[8];
            $ps_poly_center_local  = face_site[9];
            $ps_face_pts2d         = face_site[10];
            $ps_face_pts3d_local   = face_site[11];
            $ps_poly_verts_local   = face_site[12];
            $ps_poly_faces_idx     = face_site[13];
            $ps_face_planarity_err = face_site[14];
            $ps_face_is_planar     = face_site[15];
            $ps_face_family_id     = face_site[16];
            $ps_face_family_count  = face_site[17];
            $ps_edge_family_count  = face_site[18];
            $ps_vertex_family_count = face_site[19];
            $ps_face_neighbors_idx = face_site[20];
            $ps_face_dihedrals     = face_site[21];

            multmatrix(ps_frame_matrix(
                face_site[1],
                face_site[2],
                face_site[3],
                face_site[4]
            ))
                children();
        } else {
            children();
        }
    }
}

/**
 * Module: Replay caller-supplied proxy geometry for foreign sites affecting the current placed face.
 * Params: mode (foreign face fill rule), eps (tolerance), filter_parent (drop parent-edge cuts), coords (`"element"` or `"parent"`), face_child/edge_child/vertex_child (child slots)
 * Returns: none; exposes `$ps_proxy_*` metadata and calls the child slot matching the foreign source kind
 * Limitations/Gotchas: face sites are exact face-plane intrusions; edge/vertex sites are provenance-driven candidates, not distance-envelope proximity tests
 */
module place_on_face_foreign_proxy_sites(
    mode="nonzero",
    eps=1e-8,
    filter_parent=true,
    coords="element",
    face_child=0,
    edge_child=1,
    vertex_child=2
) {
    assert(!is_undef($ps_face_pts2d), "place_on_face_foreign_proxy_sites: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_idx), "place_on_face_foreign_proxy_sites: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "place_on_face_foreign_proxy_sites: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "place_on_face_foreign_proxy_sites: requires place_on_faces context ($ps_poly_verts_local)");
    assert(coords == "element" || coords == "parent", "place_on_face_foreign_proxy_sites: coords must be \"element\" or \"parent\"");
    assert(face_child >= 0 && edge_child >= 0 && vertex_child >= 0, "place_on_face_foreign_proxy_sites: child slot indices must be non-negative");

    sites = ps_face_foreign_proxy_replay_sites($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, $ps_poly_center_local, eps, mode, filter_parent);
    for (site = sites) {
        source_kind = ps_replay_site_foreign_kind(site);
        face_site = ps_replay_site_face_site(site);
        edge_site = ps_replay_site_edge_site(site);
        vertex_site = ps_replay_site_vertex_site(site);
        child_idx =
            source_kind == "face" ? face_child :
            source_kind == "edge" ? edge_child :
            source_kind == "vertex" ? vertex_child :
            undef;

        $ps_proxy_idx = ps_replay_site_idx(site);
        $ps_proxy_count = len(sites);
        $ps_proxy_kind = str("foreign_", source_kind);
        $ps_proxy_source_kind = source_kind;
        $ps_proxy_source_idx = ps_replay_site_foreign_idx(site);
        $ps_proxy_target_face_idx = $ps_face_idx;
        $ps_proxy_child_idx = child_idx;
        $ps_proxy_intrusion_record = ps_replay_site_intrusion_record(site);
        $ps_proxy_intrusion_segment2d_local = ps_replay_site_intrusion_segment2d_local(site);
        $ps_proxy_intrusion_dihedral = ps_replay_site_intrusion_dihedral(site);
        $ps_proxy_intrusion_confidence = ps_replay_site_intrusion_confidence(site);
        $ps_proxy_center_local = ps_replay_site_center_local(site);
        $ps_proxy_ex_local = ps_replay_site_ex_local(site);
        $ps_proxy_ey_local = ps_replay_site_ey_local(site);
        $ps_proxy_ez_local = ps_replay_site_ez_local(site);
        $ps_proxy_face_pts2d = ps_replay_site_face_pts2d(site);
        $ps_proxy_face_pts3d_local = ps_replay_site_face_pts3d_local(site);
        $ps_proxy_face_verts_idx = ps_replay_site_face_verts_idx(site);
        $ps_proxy_edge_pts_local = is_undef(edge_site) ? undef : edge_site[8];
        $ps_proxy_edge_verts_idx = is_undef(edge_site) ? undef : edge_site[9];
        $ps_proxy_edge_adj_faces_idx = is_undef(edge_site) ? undef : edge_site[10];
        $ps_proxy_vertex_valence = is_undef(vertex_site) ? undef : vertex_site[8];
        $ps_proxy_vertex_neighbors_idx = is_undef(vertex_site) ? undef : vertex_site[9];
        $ps_proxy_vertex_neighbor_pts_local = is_undef(vertex_site) ? undef : vertex_site[10];
        $ps_proxy_poly_verts_local = ps_replay_site_poly_verts_local(site);
        $ps_proxy_poly_center_local = ps_replay_site_poly_center_local(site);
        $ps_replay_idx = ps_replay_site_idx(site);
        $ps_replay_count = len(sites);
        $ps_replay_intrusion_record = ps_replay_site_intrusion_record(site);
        $ps_replay_kind = str("foreign_", source_kind);
        $ps_replay_foreign_kind = source_kind;
        $ps_replay_foreign_idx = ps_replay_site_foreign_idx(site);
        $ps_replay_center_local = ps_replay_site_center_local(site);
        $ps_replay_ex_local = ps_replay_site_ex_local(site);
        $ps_replay_ey_local = ps_replay_site_ey_local(site);
        $ps_replay_ez_local = ps_replay_site_ez_local(site);
        $ps_replay_face_pts2d = ps_replay_site_face_pts2d(site);
        $ps_replay_face_pts3d_local = ps_replay_site_face_pts3d_local(site);
        $ps_replay_poly_verts_local = ps_replay_site_poly_verts_local(site);
        $ps_replay_poly_center_local = ps_replay_site_poly_center_local(site);
        $ps_replay_face_verts_idx = ps_replay_site_face_verts_idx(site);
        $ps_replay_edge_pts_local = is_undef(edge_site) ? undef : edge_site[8];
        $ps_replay_edge_verts_idx = is_undef(edge_site) ? undef : edge_site[9];
        $ps_replay_edge_adj_faces_idx = is_undef(edge_site) ? undef : edge_site[10];
        $ps_replay_vertex_valence = is_undef(vertex_site) ? undef : vertex_site[8];
        $ps_replay_vertex_neighbors_idx = is_undef(vertex_site) ? undef : vertex_site[9];
        $ps_replay_vertex_neighbor_pts_local = is_undef(vertex_site) ? undef : vertex_site[10];
        $ps_replay_intrusion_segment2d_local = ps_replay_site_intrusion_segment2d_local(site);
        $ps_replay_intrusion_dihedral = ps_replay_site_intrusion_dihedral(site);
        $ps_replay_intrusion_confidence = ps_replay_site_intrusion_confidence(site);

        if (!is_undef(child_idx) && child_idx < $children) {
            if (coords == "element") {
                if (source_kind == "face") {
                    $ps_face_idx           = face_site[0];
                    $ps_edge_len           = face_site[5];
                    $ps_vertex_count       = face_site[6];
                    $ps_face_midradius     = face_site[7];
                    $ps_face_radius        = face_site[8];
                    $ps_poly_center_local  = face_site[9];
                    $ps_face_pts2d         = face_site[10];
                    $ps_face_pts3d_local   = face_site[11];
                    $ps_poly_verts_local   = face_site[12];
                    $ps_poly_faces_idx     = face_site[13];
                    $ps_face_planarity_err = face_site[14];
                    $ps_face_is_planar     = face_site[15];
                    $ps_face_family_id     = face_site[16];
                    $ps_face_family_count  = face_site[17];
                    $ps_edge_family_count  = face_site[18];
                    $ps_vertex_family_count = face_site[19];
                    $ps_face_neighbors_idx = face_site[20];
                    $ps_face_dihedrals     = face_site[21];

                    multmatrix(ps_frame_matrix(
                        face_site[1],
                        face_site[2],
                        face_site[3],
                        face_site[4]
                    ))
                        children(child_idx);
                } else if (source_kind == "edge") {
                    $ps_edge_idx            = edge_site[0];
                    $ps_edge_len            = edge_site[5];
                    $ps_edge_midradius      = edge_site[6];
                    $ps_poly_center_local   = edge_site[7];
                    $ps_edge_pts_local      = edge_site[8];
                    $ps_edge_verts_idx      = edge_site[9];
                    $ps_edge_adj_faces_idx  = edge_site[10];
                    $ps_edge_family_id      = edge_site[11];
                    $ps_face_family_count   = edge_site[12];
                    $ps_edge_family_count   = edge_site[13];
                    $ps_vertex_family_count = edge_site[14];

                    multmatrix(ps_frame_matrix(
                        edge_site[1],
                        edge_site[2],
                        edge_site[3],
                        edge_site[4]
                    ))
                        children(child_idx);
                } else if (source_kind == "vertex") {
                    $ps_vertex_idx                = vertex_site[0];
                    $ps_vertex_valence            = vertex_site[8];
                    $ps_vertex_neighbors_idx      = vertex_site[9];
                    $ps_vertex_neighbor_pts_local = vertex_site[10];
                    $ps_edge_len                  = vertex_site[5];
                    $ps_vert_radius               = vertex_site[6];
                    $ps_poly_center_local         = vertex_site[7];
                    $ps_vertex_family_id          = vertex_site[11];
                    $ps_face_family_count         = vertex_site[12];
                    $ps_edge_family_count         = vertex_site[13];
                    $ps_vertex_family_count       = vertex_site[14];

                    multmatrix(ps_frame_matrix(
                        vertex_site[1],
                        vertex_site[2],
                        vertex_site[3],
                        vertex_site[4]
                    ))
                        children(child_idx);
                }
            } else {
                children(child_idx);
            }
        }
    }
}

/**
 * Function: Build edge placement site records for `place_on_edges(...)`.
 * Params: poly (poly descriptor), inter_radius (scale input), edge_len (explicit scale override), classify/classify_opts (optional classification context)
 * Returns: list of edge site records `[edge_idx, center, ex, ey, ez, edge_len, edge_midradius, poly_center_local, edge_pts_local, edge_verts_idx, edge_adj_faces_idx, edge_family_id, face_family_count, edge_family_count, vertex_family_count]`
 * Limitations: uses an adjacent-face normal bisector for `+Z` when a usable face pair exists, with radial fallback on boundary or degenerate edges
 */
function ps_edge_sites(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef) =
    let(
        exp_edge_len = is_undef(edge_len) ? inter_radius * poly_e_over_ir(poly) : edge_len,
        scale = exp_edge_len,
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        faces0 = ps_orient_all_faces_outward(verts, faces),
        edges = _ps_edges_from_faces(faces),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [for (f = faces0) ps_face_normal(verts, f)],
        cls = _ps_resolve_classify(poly, classify, classify_opts),
        family_counts = is_undef(cls) ? undef : ps_classify_counts(cls),
        edge_family_ids = is_undef(cls) ? [] : ps_classify_edge_ids(cls, len(edges)),
        face_family_count = is_undef(family_counts) ? undef : family_counts[0],
        edge_family_count = is_undef(family_counts) ? undef : family_counts[1],
        vert_family_count = is_undef(family_counts) ? undef : family_counts[2]
    )
    [
        for (ei = [0:1:len(edges)-1])
            let(
                e = edges[ei],
                v0 = verts[e[0]] * scale,
                v1 = verts[e[1]] * scale,
                center = (v0 + v1) / 2,
                ex = v_norm(v1 - v0),
                adj_faces_idx = edge_faces[ei],
                radial_ref = v_norm(center),
                bisector_raw =
                    (len(adj_faces_idx) < 2)
                        ? radial_ref
                        : face_n[adj_faces_idx[0]] + face_n[adj_faces_idx[1]],
                bisector_signed =
                    (norm(bisector_raw) <= 1e-12)
                        ? radial_ref
                        : ((v_dot(bisector_raw, radial_ref) < 0) ? -bisector_raw : bisector_raw),
                ez_proj = bisector_signed - ex * v_dot(bisector_signed, ex),
                radial_proj = radial_ref - ex * v_dot(radial_ref, ex),
                ez_dir =
                    (norm(ez_proj) <= 1e-12)
                        ? radial_proj
                        : ez_proj,
                ez = v_norm(ez_dir),
                ey = v_norm(v_cross(ez, ex)),
                edge_midradius = norm(center),
                edge_len_actual = norm(v1 - v0),
                poly_center_local = [
                    v_dot(-center, ex),
                    v_dot(-center, ey),
                    v_dot(-center, ez)
                ],
                edge_pts_local = [[-edge_len_actual/2, 0, 0], [edge_len_actual/2, 0, 0]]
            )
            [
                ei,
                center,
                ex,
                ey,
                ez,
                edge_len_actual,
                edge_midradius,
                poly_center_local,
                edge_pts_local,
                e,
                adj_faces_idx,
                is_undef(cls) ? undef : edge_family_ids[ei],
                face_family_count,
                edge_family_count,
                vert_family_count
            ]
    ];

/**
 * Function: Build vertex placement site records for `place_on_vertices(...)`.
 * Params: poly (poly descriptor), inter_radius (scale input), edge_len (explicit scale override), classify/classify_opts (optional classification context)
 * Returns: list of vertex site records `[vertex_idx, center, ex, ey, ez, edge_len, vert_radius, poly_center_local, vertex_valence, vertex_neighbors_idx, vertex_neighbor_pts_local, vertex_family_id, face_family_count, edge_family_count, vertex_family_count]`
 * Limitations: preserves the current radial vertex frame derived from one projected neighbor direction
 */
function ps_vertex_sites(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef) =
    let(
        exp_edge_len = is_undef(edge_len) ? inter_radius * poly_e_over_ir(poly) : edge_len,
        scale = exp_edge_len,
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        edges = _ps_edges_from_faces(faces),
        cls = _ps_resolve_classify(poly, classify, classify_opts),
        family_counts = is_undef(cls) ? undef : ps_classify_counts(cls),
        vert_family_ids = is_undef(cls) ? [] : ps_classify_vert_ids(cls, len(verts)),
        face_family_count = is_undef(family_counts) ? undef : family_counts[0],
        edge_family_count = is_undef(family_counts) ? undef : family_counts[1],
        vert_family_count = is_undef(family_counts) ? undef : family_counts[2]
    )
    [
        for (vi = [0:1:len(verts)-1])
            let(
                v0 = verts[vi] * scale,
                ez = v_norm(v0),
                ni = poly_vertex_neighbor(poly, vi),
                v1 = verts[ni] * scale,
                neighbor_dir = v1 - v0,
                proj = neighbor_dir - ez * v_dot(neighbor_dir, ez),
                proj_len = norm(proj),
                ex = (proj_len == 0) ? [1,0,0] : proj / proj_len,
                ey = v_cross(ez, ex),
                center = v0,
                vert_radius = norm(center),
                neighbors_idx = _ps_vertex_site_neighbors_idx(edges, vi),
                valence = len(neighbors_idx),
                neighbor_pts_local = [
                    for (nj = neighbors_idx)
                        let(pw = verts[nj] * scale - v0)
                            [v_dot(pw, ex), v_dot(pw, ey), v_dot(pw, ez)]
                ]
            )
            [
                vi,
                center,
                ex,
                ey,
                ez,
                exp_edge_len,
                vert_radius,
                [0, 0, -vert_radius],
                valence,
                neighbors_idx,
                neighbor_pts_local,
                is_undef(cls) ? undef : vert_family_ids[vi],
                face_family_count,
                edge_family_count,
                vert_family_count
            ]
    ];

/**
 * Module: Place children on selected vertices of a polyhedron.
 * Params: poly (poly descriptor), inter_radius (scale input), edge_len (explicit scale override), classify/classify_opts (optional classification context), indices (`undef`, scalar vertex index, or list of vertex indices)
 * Returns: none; exposes `$ps_vertex_*` metadata for each selected vertex
 * Limitations: `indices` filters the placement loop only; `ps_vertex_sites(...)` still builds the complete site list so element ids and classification metadata remain global
 */
module place_on_vertices(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef, indices = undef) {
    sites = ps_vertex_sites(poly, inter_radius, edge_len, classify, classify_opts);

    for (site = sites) {
        if (_ps_place_idx_selected(site[0], indices)) {
            // Metadata for children (local-space friendly)
            $ps_vertex_idx                = site[0];
            $ps_vertex_valence            = site[8];
            $ps_vertex_neighbors_idx      = site[9];
            $ps_vertex_neighbor_pts_local = site[10];

            $ps_edge_len                  = site[5];      // (target edge length parameter)
            $ps_vert_radius               = site[6];
            $ps_poly_center_local         = site[7];
            $ps_vertex_family_id          = site[11];
            $ps_face_family_count         = site[12];
            $ps_edge_family_count         = site[13];
            $ps_vertex_family_count       = site[14];

            multmatrix(ps_frame_matrix(site[1], site[2], site[3], site[4]))
                children();
        }
    }
}


/**
 * Module: Place children on selected edges of a polyhedron.
 * Params: poly (poly descriptor), inter_radius (scale input), edge_len (explicit scale override), classify/classify_opts (optional classification context), indices (`undef`, scalar edge index, or list of edge indices)
 * Returns: none; exposes `$ps_edge_*` metadata for each selected edge
 * Limitations: `indices` filters the placement loop only; `ps_edge_sites(...)` still builds the complete site list so element ids and classification metadata remain global
 */
module place_on_edges(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef, indices = undef) {
    sites = ps_edge_sites(poly, inter_radius, edge_len, classify, classify_opts);

    for (site = sites) {
        if (_ps_place_idx_selected(site[0], indices)) {
            // Metadata for children (edge-local)
            $ps_edge_idx            = site[0];
            $ps_edge_len            = site[5];      // actual length of this edge (vs supplied edge_len = scaling factor arg)
            $ps_edge_midradius      = site[6];
            $ps_poly_center_local   = site[7];

            $ps_edge_pts_local      = site[8];
            $ps_edge_verts_idx      = site[9];
            $ps_edge_adj_faces_idx  = site[10];
            $ps_edge_family_id      = site[11];
            $ps_face_family_count   = site[12];
            $ps_edge_family_count   = site[13];
            $ps_vertex_family_count = site[14];

            multmatrix(ps_frame_matrix(site[1], site[2], site[3], site[4]))
                children();
        }
    }
}


module face_debug() {
    // Face index
    color("white") translate([0,0,2])
        text(str($ps_face_idx), size=5, halign="center", valign="center");

    // Local axes
    color("red")   cube([8,1,1], center=false);
    color("green") rotate([0,0,90]) cube([8,1,1], center=false);
    color("blue")  rotate([0,-90,0]) cube([8,1,1], center=false);

    // Radial line to centre
    color("yellow") cylinder(h = -$ps_poly_center_local[2], r = 0.5, center=false);
}
