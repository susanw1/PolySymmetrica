/*
Proxy-interaction helpers.

The idea is to keep the current analytic layer for:
- visibility / segmentation
- cut provenance
- cut_pair_id / cut_run_id
- feature classification

and add a separate fabrication-oriented path that consumes raw proposal solids.

Important distinction:
- occupancy proxies describe material that really exists in the finished part
- clearance proxies describe empty space intentionally reserved for fit/tolerance

These proxies are local interaction proxies, not already-clipped final
geometry.

The intended fabrication model is:
- `face_shielded(i) = F_raw(i) ∩ Z(i) ∩ ⋂ B(i,e)`
- `V(i) = ⋃ Occ(x)` over non-adjacent intersecting features only
- `face_final(i) = face_shielded(i) - V(i) - C_local(i)`

where:
- adjacent faces are handled by bisector anti-interference shields `B(i,e)`
- foreign cutters are non-adjacent occupied volumes only
- local face/edge/vertex seating is a separate clearance path

This file still contains lower-level experimental helpers, but that model is
the basis future work should follow rather than current fully-landed behavior.
*/

use <placement.scad>
use <funcs.scad>
use <segments.scad>
use <face_regions.scad>

function _ps_proxy_frame_inverse_matrix(center, ex, ey, ez) = [
    [ex[0], ex[1], ex[2], -v_dot(ex, center)],
    [ey[0], ey[1], ey[2], -v_dot(ey, center)],
    [ez[0], ez[1], ez[2], -v_dot(ez, center)],
    [0,     0,     0,      1]
];

function _ps_proxy_resolve_face_indices(poly, face_indices, target_face_idx) =
    is_undef(face_indices)
        ? [for (fi = [0 : 1 : len(poly_faces(poly)) - 1]) if (fi != target_face_idx) fi]
        : [for (fi = face_indices) if (fi != target_face_idx) fi];

function _ps_proxy_resolve_nonadjacent_face_indices(poly, face_indices, target_face_idx, target_face_neighbors_idx) =
    let(
        face_count = len(poly_faces(poly)),
        candidates =
            is_undef(face_indices)
                ? [for (fi = [0 : 1 : face_count - 1]) fi]
                : face_indices,
        excluded = concat(
            [target_face_idx],
            is_undef(target_face_neighbors_idx)
                ? []
                : [for (adj = target_face_neighbors_idx) if (!is_undef(adj)) adj]
        )
    )
    [for (fi = candidates) if (!_ps_list_contains(excluded, fi)) fi];

function _ps_proxy_resolve_edge_indices(poly, edge_indices) =
    is_undef(edge_indices)
        ? [for (ei = [0 : 1 : len(_ps_edges_from_faces(ps_orient_all_faces_outward(poly_verts(poly), poly_faces(poly)))) - 1]) ei]
        : edge_indices;

function _ps_proxy_resolve_vertex_indices(poly, vertex_indices) =
    is_undef(vertex_indices)
        ? [for (vi = [0 : 1 : len(poly_verts(poly)) - 1]) vi]
        : vertex_indices;

function _ps_proxy_target_face_edge_indices(poly, face_idx) =
    let(
        faces0 = ps_orient_all_faces_outward(poly_verts(poly), poly_faces(poly)),
        edges0 = _ps_edges_from_faces(faces0),
        face0 = faces0[face_idx]
    )
    [
        for (k = [0 : 1 : len(face0) - 1])
            ps_find_edge_index(edges0, face0[k], face0[(k + 1) % len(face0)])
    ];

function _ps_proxy_target_face_vertex_indices(poly, face_idx) =
    let(faces0 = ps_orient_all_faces_outward(poly_verts(poly), poly_faces(poly)))
        faces0[face_idx];

function _ps_proxy_face_edge_seg2d(face_pts2d, edge_idx) =
    let(n = len(face_pts2d))
    (is_undef(edge_idx) || edge_idx < 0 || edge_idx >= n)
        ? undef
        : [face_pts2d[edge_idx], face_pts2d[(edge_idx + 1) % n]];

function _ps_proxy_seg_shares_endpoint(seg_a, seg_b, eps=1e-8) =
    is_undef(seg_a) || is_undef(seg_b) ? false :
    _ps_seg_close2(seg_a[0], seg_b[0], eps) ||
    _ps_seg_close2(seg_a[0], seg_b[1], eps) ||
    _ps_seg_close2(seg_a[1], seg_b[0], eps) ||
    _ps_seg_close2(seg_a[1], seg_b[1], eps);

function _ps_proxy_seg_crosses_strict(seg_a, seg_b, eps=1e-8) =
    is_undef(seg_a) || is_undef(seg_b) ? false :
    _ps_proxy_seg_shares_endpoint(seg_a, seg_b, eps) ? false :
    let(
        a = seg_a[0],
        b = seg_a[1],
        c = seg_b[0],
        d = seg_b[1],
        o1 = _ps_orient2(a, b, c),
        o2 = _ps_orient2(a, b, d),
        o3 = _ps_orient2(c, d, a),
        o4 = _ps_orient2(c, d, b)
    )
    (((o1 > eps) && (o2 < -eps)) || ((o1 < -eps) && (o2 > eps))) &&
    (((o3 > eps) && (o4 < -eps)) || ((o3 < -eps) && (o4 > eps)));

function _ps_proxy_source_edges_cross(face_pts2d, edge_idx_a, edge_idx_b, eps=1e-8) =
    let(
        seg_a = _ps_proxy_face_edge_seg2d(face_pts2d, edge_idx_a),
        seg_b = _ps_proxy_face_edge_seg2d(face_pts2d, edge_idx_b)
    )
    _ps_proxy_seg_crosses_strict(seg_a, seg_b, eps);

function _ps_proxy_source_edge_cross_hit(face_pts2d, target_edge_idx, source_edge_idx, eps=1e-8) =
    let(
        seg_a = _ps_proxy_face_edge_seg2d(face_pts2d, target_edge_idx),
        seg_b = _ps_proxy_face_edge_seg2d(face_pts2d, source_edge_idx),
        hit =
            is_undef(seg_a) || is_undef(seg_b)
                ? undef
                : _ps_seg_proper_intersection(seg_a[0], seg_a[1], seg_b[0], seg_b[1], eps)
    )
    hit;

function _ps_proxy_poly_area2(pts2d) =
    let(n = len(pts2d))
    (n < 3) ? 0 :
    sum([
        for (i = [0 : 1 : n - 1])
            let(j = (i + 1) % n)
                pts2d[i][0] * pts2d[j][1] - pts2d[j][0] * pts2d[i][1]
    ]) / 2;

function _ps_proxy_sign_eps(v, eps=1e-8) =
    (v > eps) ? 1 : ((v < -eps) ? -1 : 0);

function _ps_proxy_remove_at(list, idx) =
    [for (i = [0 : 1 : len(list) - 1]) if (i != idx) list[i]];

function _ps_proxy_cross_rec_min_idx(recs, idx=0, best=0) =
    (idx >= len(recs)) ? best :
    _ps_proxy_cross_rec_min_idx(
        recs,
        idx + 1,
        (recs[idx][0] < recs[best][0]) ? idx : best
    );

function _ps_proxy_sort_cross_recs(recs, acc=[]) =
    (len(recs) == 0) ? acc :
    let(mi = _ps_proxy_cross_rec_min_idx(recs))
        _ps_proxy_sort_cross_recs(
            _ps_proxy_remove_at(recs, mi),
            concat(acc, [recs[mi]])
        );

function _ps_proxy_crossing_records(face_pts2d, target_edge_idx, source_edge_indices=undef, exclude_edge_idx=undef, eps=1e-8) =
    _ps_proxy_sort_cross_recs([
        for (source_edge_idx = [0 : 1 : len(face_pts2d) - 1])
            let(hit = _ps_proxy_source_edge_cross_hit(face_pts2d, target_edge_idx, source_edge_idx, eps))
            if (
                (is_undef(source_edge_indices) || _ps_list_contains(source_edge_indices, source_edge_idx)) &&
                (source_edge_idx != exclude_edge_idx) &&
                _ps_proxy_source_edges_cross(face_pts2d, target_edge_idx, source_edge_idx, eps) &&
                !is_undef(hit)
            )
                [hit[0], source_edge_idx, hit[2]]
    ]);

function _ps_proxy_path_edge_ids_forward(nface, edge_a, edge_b) =
    let(span = (edge_b - edge_a + nface) % nface)
    [for (step = [0 : 1 : span]) (edge_a + step) % nface];

function _ps_proxy_path_edge_ids_backward(nface, edge_a, edge_b) =
    let(span = (edge_a - edge_b + nface) % nface)
    [for (step = [0 : 1 : span]) (edge_a - step + nface) % nface];

function _ps_proxy_path_pts_forward(face_pts2d, edge_a, edge_b, hit_a, hit_b) =
    let(
        nface = len(face_pts2d),
        span = (edge_b - edge_a + nface) % nface,
        mids = [for (step = [1 : 1 : span]) face_pts2d[(edge_a + step) % nface]]
    )
    concat([hit_a], mids, [hit_b]);

function _ps_proxy_path_pts_backward(face_pts2d, edge_a, edge_b, hit_a, hit_b) =
    let(
        nface = len(face_pts2d),
        span = (edge_a - edge_b + nface) % nface,
        mids = [for (step = [0 : 1 : span - 1]) face_pts2d[(edge_a - step + nface) % nface]]
    )
    concat([hit_a], mids, [hit_b]);

function _ps_proxy_target_edge_owner_side_sign(face_pts2d, target_edge_idx, eps=1e-8) =
    let(
        nface = len(face_pts2d),
        seg = _ps_proxy_face_edge_seg2d(face_pts2d, target_edge_idx),
        s_next = _ps_proxy_sign_eps(_ps_orient2(seg[0], seg[1], face_pts2d[(target_edge_idx + 2) % nface]), eps),
        s_prev = _ps_proxy_sign_eps(_ps_orient2(seg[0], seg[1], face_pts2d[(target_edge_idx - 1 + nface) % nface]), eps)
    )
    (s_next != 0) ? s_next : s_prev;

function _ps_proxy_path_side_score(path_pts, target_seg, desired_sign, eps=1e-8) =
    sum([
        for (p = path_pts)
            let(s = _ps_proxy_sign_eps(_ps_orient2(target_seg[0], target_seg[1], p), eps))
                (s == desired_sign) ? 1 : ((s == -desired_sign) ? -1 : 0)
    ]);

function _ps_proxy_lobe_rec_from_pair(face_pts2d, target_edge_idx, rec_a, rec_b, eps=1e-8) =
    let(
        nface = len(face_pts2d),
        target_seg = _ps_proxy_face_edge_seg2d(face_pts2d, target_edge_idx),
        owner_side = _ps_proxy_target_edge_owner_side_sign(face_pts2d, target_edge_idx, eps),
        desired_side = -owner_side,
        edge_a = rec_a[1],
        edge_b = rec_b[1],
        hit_a = rec_a[2],
        hit_b = rec_b[2],
        poly_fwd = _ps_proxy_path_pts_forward(face_pts2d, edge_a, edge_b, hit_a, hit_b),
        poly_bwd = _ps_proxy_path_pts_backward(face_pts2d, edge_a, edge_b, hit_a, hit_b),
        fwd_score = _ps_proxy_path_side_score(poly_fwd, target_seg, desired_side, eps),
        bwd_score = _ps_proxy_path_side_score(poly_bwd, target_seg, desired_side, eps),
        area_fwd = abs(_ps_proxy_poly_area2(poly_fwd)),
        area_bwd = abs(_ps_proxy_poly_area2(poly_bwd)),
        use_fwd =
            (fwd_score > bwd_score) ? true :
            (bwd_score > fwd_score) ? false :
            (area_fwd <= area_bwd),
        poly2d = use_fwd ? poly_fwd : poly_bwd,
        edge_ids = use_fwd
            ? _ps_proxy_path_edge_ids_forward(nface, edge_a, edge_b)
            : _ps_proxy_path_edge_ids_backward(nface, edge_a, edge_b)
    )
    [poly2d, edge_ids, [edge_a, edge_b], [rec_a[0], rec_b[0]]];

function _ps_proxy_target_edge_lobe_recs(face_pts2d, target_edge_idx, source_edge_indices=undef, exclude_edge_idx=undef, eps=1e-8) =
    let(
        recs = _ps_proxy_crossing_records(face_pts2d, target_edge_idx, source_edge_indices, exclude_edge_idx, eps),
        _even = assert((len(recs) % 2) == 0, str("_ps_proxy_target_edge_lobe_recs: odd crossing count for edge ", target_edge_idx))
    )
    (len(recs) < 2)
        ? []
        : [
            for (ri = [0 : 2 : len(recs) - 2])
                _ps_proxy_lobe_rec_from_pair(face_pts2d, target_edge_idx, recs[ri], recs[ri + 1], eps)
        ];

function _ps_proxy_cell_color(mask) = [
    0.3 + 0.6 * abs(sin(17 * mask + 0.1)),
    0.3 + 0.6 * abs(sin(29 * mask + 0.6)),
    0.3 + 0.6 * abs(sin(43 * mask + 1.1)),
    1
];

module _ps_proxy_face_influence(face_bounds) {
    z0 = face_bounds[0];
    z1 = face_bounds[1];
    assert(is_list(face_bounds) && len(face_bounds) >= 2 && z1 >= z0, "face_bounds must be [z0, z1] with z1 >= z0");

    translate([0, 0, z0])
        linear_extrude(height = z1 - z0)
            square([1000, 1000], center = true);
}

module _ps_proxy_face_interferer(face_bounds, mode = "raw") {
    if (mode == "raw")
        intersection() {
            children();
            _ps_proxy_face_influence(face_bounds);
        }
    else if (mode == "sweep_to_bounds") {
        z0 = face_bounds[0];
        z1 = face_bounds[1];
        translate([0, 0, z0])
            linear_extrude(height = z1 - z0)
                projection(cut = false)
                    children();
    } else
        assert(false, str("_ps_proxy_face_interferer: unknown mode ", mode));
}

module _ps_proxy_edge_influence(edge_radius, edge_length) {
    if (is_undef(edge_radius) || is_undef(edge_length))
        children();
    else
        intersection() {
            children();
            translate([-edge_length / 2, -edge_radius, -edge_radius])
                cube([edge_length, 2 * edge_radius, 2 * edge_radius]);
        }
}

module _ps_proxy_vertex_influence(vertex_radius) {
    if (is_undef(vertex_radius))
        children();
    else
        intersection() {
            children();
            sphere(r = vertex_radius);
        }
}

module _ps_proxy_place_faces_in_target_frame(poly, target_edge_len, target_inv, face_bounds, face_indices, face_proxy_mode) {
    if (!is_undef(face_indices) && len(face_indices) > 0)
        multmatrix(target_inv)
            place_on_faces(poly, edge_len = target_edge_len, indices = face_indices)
                _ps_proxy_face_interferer(face_bounds, face_proxy_mode)
                    children();
}

module _ps_proxy_place_edges_in_target_frame(poly, target_edge_len, target_inv, edge_radius, edge_length, edge_indices) {
    multmatrix(target_inv)
        place_on_edges(poly, edge_len = target_edge_len, indices = edge_indices)
            _ps_proxy_edge_influence(edge_radius, edge_length)
                children();
}

module ps_clip_local_boundary_edge_by_own_interference_ctx(edge_radius, edge_length, face_bounds, eps = 1e-8) {
    z0 = face_bounds[0];
    z1 = face_bounds[1];

    difference() {
        _ps_proxy_edge_influence(edge_radius, edge_length)
            children();

        ps_face_boundary_seg_default_interference_cutter_local_ctx(z0, z1, eps);
    }
}

module _ps_proxy_emit_crossing_local_boundary_lobe_bodies_in_current_edge_frame(face_bounds, edge_radius, edge_length, source_edge_indices, lobe_indices, exclude_seg_idx, eps = 1e-8) {
    z0 = face_bounds[0];
    z1 = face_bounds[1];
    target_face_pts2d = $ps_face_pts2d;
    target_source_edge_idx = $ps_face_boundary_seg_source_edge_idx;
    target_inv = _ps_proxy_frame_inverse_matrix(
        $ps_edge_center_world,
        $ps_edge_ex_world,
        $ps_edge_ey_world,
        $ps_edge_ez_world
    );
    lobes = _ps_proxy_target_edge_lobe_recs(target_face_pts2d, target_source_edge_idx, source_edge_indices, exclude_seg_idx, eps);

    for (li = [0 : 1 : len(lobes) - 1])
        if (is_undef(lobe_indices) || _ps_list_contains(lobe_indices, li))
            intersection() {
                multmatrix(target_inv)
                    ps_face_crossing_lobe_body_ctx(lobes[li][0], lobes[li][1], z0, z1, eps);

                ps_clip_local_boundary_edge_by_own_interference_ctx(
                    edge_radius,
                    edge_length,
                    face_bounds,
                    eps
                )
                    children();
            }
}

module _ps_proxy_emit_crossing_local_boundary_lobe_bodies_raw_in_current_edge_frame(face_bounds, source_edge_indices, lobe_indices, exclude_seg_idx, eps = 1e-8) {
    z0 = face_bounds[0];
    z1 = face_bounds[1];
    target_face_pts2d = $ps_face_pts2d;
    target_source_edge_idx = $ps_face_boundary_seg_source_edge_idx;
    target_inv = _ps_proxy_frame_inverse_matrix(
        $ps_edge_center_world,
        $ps_edge_ex_world,
        $ps_edge_ey_world,
        $ps_edge_ez_world
    );
    lobes = _ps_proxy_target_edge_lobe_recs(target_face_pts2d, target_source_edge_idx, source_edge_indices, exclude_seg_idx, eps);

    multmatrix(target_inv)
        for (li = [0 : 1 : len(lobes) - 1])
            if (is_undef(lobe_indices) || _ps_list_contains(lobe_indices, li))
                ps_face_crossing_lobe_body_ctx(lobes[li][0], lobes[li][1], z0, z1, eps);
}

module ps_clip_local_boundary_edge_by_crossing_edges_ctx(edge_radius, edge_length, face_bounds, source_edge_indices = undef, lobe_indices = undef, eps = 1e-8) {
    assert(!is_undef($ps_face_boundary_seg_idx), "ps_clip_local_boundary_edge_by_crossing_edges_ctx must be used inside place_on_face_filled_boundary_source_edges(...)");
    assert(!is_undef($ps_edge_center_world), "ps_clip_local_boundary_edge_by_crossing_edges_ctx must be used inside place_on_face_filled_boundary_source_edges(...)");

    difference() {
        ps_clip_local_boundary_edge_by_own_interference_ctx(
            edge_radius,
            edge_length,
            face_bounds,
            eps
        )
            children();

        _ps_proxy_emit_crossing_local_boundary_lobe_bodies_in_current_edge_frame(
            face_bounds,
            edge_radius,
            edge_length,
            source_edge_indices,
            lobe_indices,
            $ps_face_boundary_seg_source_edge_idx,
            eps
        )
            children();
    }
}

module _ps_proxy_place_local_boundary_edges_in_target_face(edge_radius, edge_length, source_edge_indices, face_bounds, eps = 1e-8) {
    place_on_face_filled_boundary_source_edges(source_edge_indices = source_edge_indices, eps = eps)
        ps_clip_local_boundary_edge_by_crossing_edges_ctx(
            edge_radius,
            edge_length,
            face_bounds,
            source_edge_indices,
            undef,
            eps
        )
            children();
}

module _ps_proxy_place_vertices_in_target_frame(poly, target_edge_len, target_inv, vertex_radius, vertex_indices) {
    multmatrix(target_inv)
        place_on_vertices(poly, edge_len = target_edge_len, indices = vertex_indices)
            _ps_proxy_vertex_influence(vertex_radius)
                children();
}

module _ps_proxy_emit_cutter_union_in_target_frame(
    poly,
    target_edge_len,
    target_inv,
    face_bounds,
    face_proxy_mode,
    edge_radius,
    edge_length,
    vertex_radius,
    resolved_face_indices,
    resolved_edge_indices,
    resolved_vertex_indices
) {
    union() {
        if ($children > 0)
            _ps_proxy_place_faces_in_target_frame(
                poly,
                target_edge_len,
                target_inv,
                face_bounds,
                resolved_face_indices,
                face_proxy_mode
            )
                children(0);

        if ($children > 1)
            _ps_proxy_place_edges_in_target_frame(
                poly,
                target_edge_len,
                target_inv,
                edge_radius,
                edge_length,
                resolved_edge_indices
            )
                children(1);

        if ($children > 2)
            _ps_proxy_place_vertices_in_target_frame(
                poly,
                target_edge_len,
                target_inv,
                vertex_radius,
                resolved_vertex_indices
            )
                children(2);
    }
}

module _ps_proxy_emit_face_cutter_by_index(
    poly,
    target_edge_len,
    target_inv,
    face_bounds,
    face_proxy_mode,
    edge_radius,
    edge_length,
    vertex_radius,
    resolved_face_indices,
    resolved_edge_indices,
    resolved_vertex_indices,
    cutter_idx
) {
    face_count = len(resolved_face_indices);
    edge_count = len(resolved_edge_indices);

    if (cutter_idx < face_count) {
        if ($children > 1)
            _ps_proxy_place_faces_in_target_frame(
                poly,
                target_edge_len,
                target_inv,
                face_bounds,
                [resolved_face_indices[cutter_idx]],
                face_proxy_mode
            )
                children(1);
    } else if (cutter_idx < face_count + edge_count) {
        if ($children > 2)
            _ps_proxy_place_edges_in_target_frame(
                poly,
                target_edge_len,
                target_inv,
                edge_radius,
                edge_length,
                [resolved_edge_indices[cutter_idx - face_count]]
            )
                children(2);
    } else {
        if ($children > 3)
            _ps_proxy_place_vertices_in_target_frame(
                poly,
                target_edge_len,
                target_inv,
                vertex_radius,
                [resolved_vertex_indices[cutter_idx - face_count - edge_count]]
            )
                children(3);
    }
}

module _ps_proxy_partition_face_cells(
    poly,
    target_edge_len,
    target_inv,
    face_bounds,
    face_proxy_mode,
    edge_radius,
    edge_length,
    vertex_radius,
    resolved_face_indices,
    resolved_edge_indices,
    resolved_vertex_indices,
    cutter_count,
    idx = 0,
    mask = 0,
    color_cells = false
) {
    if (idx >= cutter_count)
        let($ps_proxy_cell_mask = mask)
            if (color_cells)
                color(_ps_proxy_cell_color(mask))
                    children(0);
            else
                children(0);
    else {
        difference() {
            _ps_proxy_partition_face_cells(
                poly,
                target_edge_len,
                target_inv,
                face_bounds,
                face_proxy_mode,
                edge_radius,
                edge_length,
                vertex_radius,
                resolved_face_indices,
                resolved_edge_indices,
                resolved_vertex_indices,
                cutter_count,
                idx + 1,
                mask,
                color_cells
            ) {
                children(0);
                if ($children > 1) children(1);
                if ($children > 2) children(2);
                if ($children > 3) children(3);
            }
            _ps_proxy_emit_face_cutter_by_index(
                poly,
                target_edge_len,
                target_inv,
                face_bounds,
                face_proxy_mode,
                edge_radius,
                edge_length,
                vertex_radius,
                resolved_face_indices,
                resolved_edge_indices,
                resolved_vertex_indices,
                idx
            ) {
                children(0);
                if ($children > 1) children(1);
                if ($children > 2) children(2);
                if ($children > 3) children(3);
            }
        }
        intersection() {
            _ps_proxy_partition_face_cells(
                poly,
                target_edge_len,
                target_inv,
                face_bounds,
                face_proxy_mode,
                edge_radius,
                edge_length,
                vertex_radius,
                resolved_face_indices,
                resolved_edge_indices,
                resolved_vertex_indices,
                cutter_count,
                idx + 1,
                mask + pow(2, idx),
                color_cells
            ) {
                children(0);
                if ($children > 1) children(1);
                if ($children > 2) children(2);
                if ($children > 3) children(3);
            }
            _ps_proxy_emit_face_cutter_by_index(
                poly,
                target_edge_len,
                target_inv,
                face_bounds,
                face_proxy_mode,
                edge_radius,
                edge_length,
                vertex_radius,
                resolved_face_indices,
                resolved_edge_indices,
                resolved_vertex_indices,
                idx
            ) {
                children(0);
                if ($children > 1) children(1);
                if ($children > 2) children(2);
                if ($children > 3) children(3);
            }
        }
    }
}

/*
Clip one face-local occupancy proxy by raw non-adjacent neighboring face /
edge / vertex occupancy proxies, after applying the local anti-interference
shield to that target face.

Children:
- child 0: face occupancy proxy
  This convenience wrapper also reuses child 0 as the foreign face cutter
  proxy. Use `ps_clip_interference_target_by_feature_proxies(...)` if target
  and foreign face proxies differ.
- child 1: edge occupancy proxy
- child 2: vertex occupancy proxy
*/
module ps_clip_face_by_feature_proxies(
    poly,
    face_idx,
    inter_radius = 1,
    edge_len = undef,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    ps_clip_interference_target_by_feature_proxies(
        poly,
        face_idx,
        inter_radius = inter_radius,
        edge_len = edge_len,
        face_bounds = face_bounds,
        face_proxy_mode = face_proxy_mode,
        edge_radius = edge_radius,
        edge_length = edge_length,
        vertex_radius = vertex_radius,
        include_faces = include_faces,
        include_edges = include_edges,
        include_vertices = include_vertices,
        face_indices = face_indices,
        edge_indices = edge_indices,
        vertex_indices = vertex_indices,
        filter = filter,
        eps = eps
    ) {
        if ($children > 0) children(0);
        if ($children > 0) children(0);
        if ($children > 1) children(1);
        if ($children > 2) children(2);
    }
}

/*
Clip one arbitrary face-local target proxy by raw non-adjacent neighboring
face / edge / vertex occupancy proxies, after applying the local
anti-interference shield to that target face.

Children:
- child 0: target proxy to keep and clip
- child 1: face occupancy cutter proxy
- child 2: edge occupancy cutter proxy
- child 3: vertex occupancy cutter proxy
*/
module ps_clip_interference_target_by_feature_proxies(
    poly,
    face_idx,
    inter_radius = 1,
    edge_len = undef,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    place_on_faces(poly, inter_radius = inter_radius, edge_len = edge_len, indices = [face_idx])
        ps_clip_interference_target_by_feature_proxies_ctx(
            poly = poly,
            face_bounds = face_bounds,
            face_proxy_mode = face_proxy_mode,
            edge_radius = edge_radius,
            edge_length = edge_length,
            vertex_radius = vertex_radius,
            include_faces = include_faces,
            include_edges = include_edges,
            include_vertices = include_vertices,
            face_indices = face_indices,
            edge_indices = edge_indices,
            vertex_indices = vertex_indices,
            filter = filter,
            eps = eps
        )
        {
            if ($children > 0) children(0);
            if ($children > 1) children(1);
            if ($children > 2) children(2);
            if ($children > 3) children(3);
        }
}

/*
Context wrapper variant intended for use inside `place_on_faces(...)`.
*/
module ps_clip_interference_target_by_feature_proxies_ctx(
    poly,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    target_edge_len = $ps_edge_len;
    target_inv = _ps_proxy_frame_inverse_matrix(
        $ps_face_center_world,
        $ps_face_ex_world,
        $ps_face_ey_world,
        $ps_face_ez_world
    );
    resolved_face_indices = _ps_proxy_resolve_nonadjacent_face_indices(
        poly,
        face_indices,
        $ps_face_idx,
        $ps_face_neighbors_idx
    );
    resolved_edge_indices = _ps_proxy_resolve_edge_indices(poly, edge_indices);
    resolved_vertex_indices = _ps_proxy_resolve_vertex_indices(poly, vertex_indices);
    assert(!is_undef(target_edge_len), "ps_clip_interference_target_by_feature_proxies_ctx must be used inside place_on_faces(...)");
    assert(is_undef(filter), "ps_clip_interference_target_by_feature_proxies_ctx: filter is reserved and not yet implemented");

    difference() {
        intersection() {
            children(0);
            ps_face_interference_volume_ctx(face_bounds[0], face_bounds[1], eps);
        }

        _ps_proxy_emit_cutter_union_in_target_frame(
            poly,
            target_edge_len,
            target_inv,
            face_bounds,
            face_proxy_mode,
            include_edges ? edge_radius : undef,
            include_edges ? edge_length : undef,
            include_vertices ? vertex_radius : undef,
            include_faces ? resolved_face_indices : [],
            include_edges ? resolved_edge_indices : [],
            include_vertices ? resolved_vertex_indices : []
        ) {
            if ($children > 1) children(1);
            if ($children > 2) children(2);
            if ($children > 3) children(3);
        }
    }
}

/*
Clip one arbitrary face-local target proxy by raw non-adjacent neighboring
face / edge / vertex occupancy proxies.

Children:
- child 0: target proxy to keep and clip
- child 1: face occupancy cutter proxy
- child 2: edge occupancy cutter proxy
- child 3: vertex occupancy cutter proxy
*/
module ps_clip_target_by_feature_proxies(
    poly,
    face_idx,
    inter_radius = 1,
    edge_len = undef,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    place_on_faces(poly, inter_radius = inter_radius, edge_len = edge_len, indices = [face_idx])
        ps_clip_target_by_feature_proxies_ctx(
            poly = poly,
            face_bounds = face_bounds,
            face_proxy_mode = face_proxy_mode,
            edge_radius = edge_radius,
            edge_length = edge_length,
            vertex_radius = vertex_radius,
            include_faces = include_faces,
            include_edges = include_edges,
            include_vertices = include_vertices,
            face_indices = face_indices,
            edge_indices = edge_indices,
            vertex_indices = vertex_indices,
            filter = filter,
            eps = eps
        )
        {
            if ($children > 0) children(0);
            if ($children > 1) children(1);
            if ($children > 2) children(2);
            if ($children > 3) children(3);
        }
}

/*
Context wrapper variant intended for use inside `place_on_faces(...)`.
*/
module ps_clip_target_by_feature_proxies_ctx(
    poly,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    target_edge_len = $ps_edge_len;
    target_inv = _ps_proxy_frame_inverse_matrix(
        $ps_face_center_world,
        $ps_face_ex_world,
        $ps_face_ey_world,
        $ps_face_ez_world
    );
    resolved_face_indices = _ps_proxy_resolve_face_indices(poly, face_indices, $ps_face_idx);
    resolved_edge_indices = _ps_proxy_resolve_edge_indices(poly, edge_indices);
    resolved_vertex_indices = _ps_proxy_resolve_vertex_indices(poly, vertex_indices);
    assert(!is_undef(target_edge_len), "ps_clip_target_by_feature_proxies_ctx must be used inside place_on_faces(...)");
    assert(is_undef(filter), "ps_clip_target_by_feature_proxies_ctx: filter is reserved and not yet implemented");

    difference() {
        intersection() {
            children(0);
            _ps_proxy_face_influence(face_bounds);
        }

        _ps_proxy_emit_cutter_union_in_target_frame(
            poly,
            target_edge_len,
            target_inv,
            face_bounds,
            face_proxy_mode,
            include_edges ? edge_radius : undef,
            include_edges ? edge_length : undef,
            include_vertices ? vertex_radius : undef,
            include_faces ? resolved_face_indices : [],
            include_edges ? resolved_edge_indices : [],
            include_vertices ? resolved_vertex_indices : []
        ) {
            if ($children > 1) children(1);
            if ($children > 2) children(2);
            if ($children > 3) children(3);
        }
    }
}

/*
Clip local edge / vertex clearance geometry for one face by foreign occupancy.

Children:
- child 0: local edge clearance proxy (placed on the target face's true filled-boundary edges)
- child 1: local vertex clearance proxy (placed on the target face vertices)
- child 2: foreign face occupancy cutter proxy
- child 3: foreign edge occupancy cutter proxy
- child 4: foreign vertex occupancy cutter proxy
*/
module ps_clip_face_local_clearance_by_feature_proxies(
    poly,
    face_idx,
    inter_radius = 1,
    edge_len = undef,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_local_edges = true,
    include_local_vertices = false,
    include_cutter_faces = true,
    include_cutter_edges = true,
    include_cutter_vertices = true,
    local_edge_indices = undef,
    local_vertex_indices = undef,
    cutter_face_indices = undef,
    cutter_edge_indices = undef,
    cutter_vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    place_on_faces(poly, inter_radius = inter_radius, edge_len = edge_len, indices = [face_idx])
        ps_clip_face_local_clearance_by_feature_proxies_ctx(
            poly = poly,
            face_bounds = face_bounds,
            face_proxy_mode = face_proxy_mode,
            edge_radius = edge_radius,
            edge_length = edge_length,
            vertex_radius = vertex_radius,
            include_local_edges = include_local_edges,
            include_local_vertices = include_local_vertices,
            include_cutter_faces = include_cutter_faces,
            include_cutter_edges = include_cutter_edges,
            include_cutter_vertices = include_cutter_vertices,
            local_edge_indices = local_edge_indices,
            local_vertex_indices = local_vertex_indices,
            cutter_face_indices = cutter_face_indices,
            cutter_edge_indices = cutter_edge_indices,
            cutter_vertex_indices = cutter_vertex_indices,
            filter = filter,
            eps = eps
        ) {
            if ($children > 0) children(0);
            if ($children > 1) children(1);
            if ($children > 2) children(2);
            if ($children > 3) children(3);
            if ($children > 4) children(4);
        }
}

module ps_clip_face_local_clearance_by_feature_proxies_ctx(
    poly,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_local_edges = true,
    include_local_vertices = false,
    include_cutter_faces = true,
    include_cutter_edges = true,
    include_cutter_vertices = true,
    local_edge_indices = undef,
    local_vertex_indices = undef,
    cutter_face_indices = undef,
    cutter_edge_indices = undef,
    cutter_vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    target_edge_len = $ps_edge_len;
    target_inv = _ps_proxy_frame_inverse_matrix(
        $ps_face_center_world,
        $ps_face_ex_world,
        $ps_face_ey_world,
        $ps_face_ez_world
    );
    resolved_local_edge_indices =
        include_local_edges && $children > 0
            ? (is_undef(local_edge_indices) ? _ps_proxy_target_face_edge_indices(poly, $ps_face_idx) : local_edge_indices)
            : [];
    resolved_local_vertex_indices =
        include_local_vertices && $children > 1
            ? (is_undef(local_vertex_indices) ? _ps_proxy_target_face_vertex_indices(poly, $ps_face_idx) : local_vertex_indices)
            : [];
    resolved_cutter_face_indices =
        include_cutter_faces && $children > 2
            ? _ps_proxy_resolve_face_indices(poly, cutter_face_indices, $ps_face_idx)
            : [];
    resolved_cutter_edge_indices =
        include_cutter_edges && $children > 3
            ? _ps_proxy_resolve_edge_indices(poly, cutter_edge_indices)
            : [];
    resolved_cutter_vertex_indices =
        include_cutter_vertices && $children > 4
            ? _ps_proxy_resolve_vertex_indices(poly, cutter_vertex_indices)
            : [];
    assert(!is_undef(target_edge_len), "ps_clip_face_local_clearance_by_feature_proxies_ctx must be used inside place_on_faces(...)");
    assert(is_undef(filter), "ps_clip_face_local_clearance_by_feature_proxies_ctx: filter is reserved and not yet implemented");

    difference() {
        union() {
            if ($children > 0)
                _ps_proxy_place_local_boundary_edges_in_target_face(
                    edge_radius,
                    edge_length,
                    resolved_local_edge_indices,
                    face_bounds,
                    eps
                )
                    children(0);

            if ($children > 1)
                _ps_proxy_place_vertices_in_target_frame(
                    poly,
                    target_edge_len,
                    target_inv,
                    vertex_radius,
                    resolved_local_vertex_indices
                )
                    children(1);
        }

        _ps_proxy_emit_cutter_union_in_target_frame(
            poly,
            target_edge_len,
            target_inv,
            face_bounds,
            face_proxy_mode,
            include_cutter_edges ? edge_radius : undef,
            include_cutter_edges ? edge_length : undef,
            include_cutter_vertices ? vertex_radius : undef,
            resolved_cutter_face_indices,
            resolved_cutter_edge_indices,
            resolved_cutter_vertex_indices
        ) {
            if ($children > 2) children(2);
            if ($children > 3) children(3);
            if ($children > 4) children(4);
        }
    }
}

/*
Build a printable face by clipping shielded face occupancy against foreign
occupancy, then subtracting the clipped local clearance.

Children:
- child 0: face occupancy proxy
- child 1: local edge clearance proxy
- child 2: local vertex clearance proxy
- child 3: foreign face occupancy cutter proxy
- child 4: foreign edge occupancy cutter proxy
- child 5: foreign vertex occupancy cutter proxy
*/
module ps_carve_face_by_feature_proxies(
    poly,
    face_idx,
    inter_radius = 1,
    edge_len = undef,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_local_edges = true,
    include_local_vertices = false,
    include_cutter_faces = true,
    include_cutter_edges = true,
    include_cutter_vertices = true,
    cutter_face_indices = undef,
    cutter_edge_indices = undef,
    cutter_vertex_indices = undef,
    local_edge_indices = undef,
    local_vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    difference() {
        ps_clip_interference_target_by_feature_proxies(
            poly,
            face_idx,
            inter_radius = inter_radius,
            edge_len = edge_len,
            face_bounds = face_bounds,
            face_proxy_mode = face_proxy_mode,
            edge_radius = edge_radius,
            edge_length = edge_length,
            vertex_radius = vertex_radius,
            include_faces = include_cutter_faces,
            include_edges = include_cutter_edges,
            include_vertices = include_cutter_vertices,
            face_indices = cutter_face_indices,
            edge_indices = cutter_edge_indices,
            vertex_indices = cutter_vertex_indices,
            filter = filter,
            eps = eps
        ) {
            if ($children > 0) children(0);
            if ($children > 3) children(3);
            if ($children > 4) children(4);
            if ($children > 5) children(5);
        }

        ps_clip_face_local_clearance_by_feature_proxies(
            poly,
            face_idx,
            inter_radius = inter_radius,
            edge_len = edge_len,
            face_bounds = face_bounds,
            face_proxy_mode = face_proxy_mode,
            edge_radius = edge_radius,
            edge_length = edge_length,
            vertex_radius = vertex_radius,
            include_local_edges = include_local_edges,
            include_local_vertices = include_local_vertices,
            include_cutter_faces = include_cutter_faces,
            include_cutter_edges = include_cutter_edges,
            include_cutter_vertices = include_cutter_vertices,
            local_edge_indices = local_edge_indices,
            local_vertex_indices = local_vertex_indices,
            cutter_face_indices = cutter_face_indices,
            cutter_edge_indices = cutter_edge_indices,
            cutter_vertex_indices = cutter_vertex_indices,
            filter = filter,
            eps = eps
        ) {
            if ($children > 1) children(1);
            if ($children > 2) children(2);
            if ($children > 3) children(3);
            if ($children > 4) children(4);
            if ($children > 5) children(5);
        }
    }
}

/*
Partition one face-local proxy into cells created by raw neighboring face /
edge / vertex proxies. The leaf cells are emitted as separate solids.

Children:
- child 0: face occupancy proxy
- child 1: edge occupancy proxy
- child 2: vertex occupancy proxy
*/
module ps_partition_face_by_feature_proxies(
    poly,
    face_idx,
    inter_radius = 1,
    edge_len = undef,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    color_cells = false,
    max_cutters = 10,
    filter = undef,
    eps = 1e-8
) {
    place_on_faces(poly, inter_radius = inter_radius, edge_len = edge_len, indices = [face_idx])
        ps_partition_face_by_feature_proxies_ctx(
            poly = poly,
            face_bounds = face_bounds,
            face_proxy_mode = face_proxy_mode,
            edge_radius = edge_radius,
            edge_length = edge_length,
            vertex_radius = vertex_radius,
            include_faces = include_faces,
            include_edges = include_edges,
            include_vertices = include_vertices,
            face_indices = face_indices,
            edge_indices = edge_indices,
            vertex_indices = vertex_indices,
            color_cells = color_cells,
            max_cutters = max_cutters,
            filter = filter,
            eps = eps
        )
        {
            if ($children > 0) children(0);
            if ($children > 1) children(1);
            if ($children > 2) children(2);
        }
}

module ps_partition_face_by_feature_proxies_ctx(
    poly,
    face_bounds = [-1, 1],
    face_proxy_mode = "raw",
    edge_radius = undef,
    edge_length = undef,
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    color_cells = false,
    max_cutters = 10,
    filter = undef,
    eps = 1e-8
) {
    target_edge_len = $ps_edge_len;
    target_inv = _ps_proxy_frame_inverse_matrix(
        $ps_face_center_world,
        $ps_face_ex_world,
        $ps_face_ey_world,
        $ps_face_ez_world
    );
    resolved_face_indices = include_faces && $children > 0 ? _ps_proxy_resolve_face_indices(poly, face_indices, $ps_face_idx) : [];
    resolved_edge_indices = include_edges && $children > 1 ? _ps_proxy_resolve_edge_indices(poly, edge_indices) : [];
    resolved_vertex_indices = include_vertices && $children > 2 ? _ps_proxy_resolve_vertex_indices(poly, vertex_indices) : [];
    cutter_count = len(resolved_face_indices) + len(resolved_edge_indices) + len(resolved_vertex_indices);
    assert(!is_undef(target_edge_len), "ps_partition_face_by_feature_proxies_ctx must be used inside place_on_faces(...)");
    assert(is_undef(filter), "ps_partition_face_by_feature_proxies_ctx: filter is reserved and not yet implemented");
    assert(cutter_count <= max_cutters, str("ps_partition_face_by_feature_proxies_ctx: cutter_count=", cutter_count, " exceeds max_cutters=", max_cutters, "; pass explicit indices or raise the limit"));

    _ps_proxy_partition_face_cells(
        poly,
        target_edge_len,
        target_inv,
        face_bounds,
        face_proxy_mode,
        edge_radius,
        edge_length,
        vertex_radius,
        resolved_face_indices,
        resolved_edge_indices,
        resolved_vertex_indices,
        cutter_count,
        color_cells = color_cells
    ) {
        intersection() {
            children(0);
            _ps_proxy_face_influence(face_bounds);
        }
        if ($children > 0) children(0);
        if ($children > 1) children(1);
        if ($children > 2) children(2);
    }
}

/*
Future edge-space analogue. Intended to clip an edge-local proxy against raw
neighboring face / edge / vertex proxies.
*/
module ps_clip_edge_by_feature_proxies(
    poly,
    edge_idx,
    inter_radius = 1,
    edge_len = undef,
    edge_radius,
    edge_length,
    face_bounds = [-1, 1],
    vertex_radius = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    assert(false, "ps_clip_edge_by_feature_proxies: API sketch only; not implemented");
}

/*
Future vertex-space analogue. Intended to clip a vertex-local proxy against raw
neighboring face / edge / vertex proxies.
*/
module ps_clip_vertex_by_feature_proxies(
    poly,
    vertex_idx,
    inter_radius = 1,
    edge_len = undef,
    vertex_radius,
    face_bounds = [-1, 1],
    edge_radius = undef,
    edge_length = undef,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = undef,
    edge_indices = undef,
    vertex_indices = undef,
    filter = undef,
    eps = 1e-8
) {
    assert(false, "ps_clip_vertex_by_feature_proxies: API sketch only; not implemented");
}
