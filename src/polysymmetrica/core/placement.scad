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
                face_pts2d = [for (p = face_pts3d_local) [p[0], p[1]]],
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

// ---- Generic face-placement driver ----
module place_on_faces(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef) {
    sites = ps_face_sites(poly, inter_radius, edge_len, classify, classify_opts);

    for (site = sites) {
        fi = site[0];
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
        edges = _ps_edges_from_faces(faces0),
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

// ---- Place children on all vertices of a polyhedron ----
module place_on_vertices(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef) {
    sites = ps_vertex_sites(poly, inter_radius, edge_len, classify, classify_opts);

    for (site = sites) {
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


// ---- Place children on all edges of a polyhedron ----
module place_on_edges(poly, inter_radius = 1, edge_len = undef, classify = undef, classify_opts = undef) {
    sites = ps_edge_sites(poly, inter_radius, edge_len, classify, classify_opts);

    for (site = sites) {
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
