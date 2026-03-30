/*
Proxy-interaction API sketch.

The idea is to keep the current analytic layer for:
- visibility / segmentation
- cut provenance
- cut_pair_id / cut_run_id
- feature classification

and add a separate fabrication-oriented path that consumes raw proposal solids:
- child 0: face proxy
- child 1: edge proxy
- child 2: vertex proxy

These proxies are local interaction proxies, not already-clipped final
geometry. The face path below is intentionally non-recursive: neighboring raw
feature proxies are instantiated in their own exact placement frames, clipped to
simple local influence bounds, transformed into the target face-local frame, and
subtracted from the target face proxy.
*/

use <placement.scad>
use <funcs.scad>

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

module _ps_proxy_place_vertices_in_target_frame(poly, target_edge_len, target_inv, vertex_radius, vertex_indices) {
    multmatrix(target_inv)
        place_on_vertices(poly, edge_len = target_edge_len, indices = vertex_indices)
            _ps_proxy_vertex_influence(vertex_radius)
                children();
}

/*
Clip one face-local proxy by raw neighboring face / edge / vertex proxies.

Children:
- child 0: face proxy
- child 1: edge proxy
- child 2: vertex proxy
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
    place_on_faces(poly, inter_radius = inter_radius, edge_len = edge_len, indices = [face_idx])
        ps_clip_face_by_feature_proxies_ctx(
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
        }
}

/*
Context wrapper variant intended for use inside `place_on_faces(...)`.
*/
module ps_clip_face_by_feature_proxies_ctx(
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

    assert(!is_undef(target_edge_len), "ps_clip_face_by_feature_proxies_ctx must be used inside place_on_faces(...)");
    assert(is_undef(filter), "ps_clip_face_by_feature_proxies_ctx: filter is reserved and not yet implemented");

    difference() {
        intersection() {
            children(0);
            _ps_proxy_face_influence(face_bounds);
        }

        union() {
            if (include_faces && $children > 0)
                _ps_proxy_place_faces_in_target_frame(poly, target_edge_len, target_inv, face_bounds, resolved_face_indices, face_proxy_mode)
                    children(0);

            if (include_edges && $children > 1)
                _ps_proxy_place_edges_in_target_frame(poly, target_edge_len, target_inv, edge_radius, edge_length, edge_indices)
                    children(1);

            if (include_vertices && $children > 2)
                _ps_proxy_place_vertices_in_target_frame(poly, target_edge_len, target_inv, vertex_radius, vertex_indices)
                    children(2);
        }
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
