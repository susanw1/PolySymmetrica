use <../../core/prisms.scad>
use <../../core/funcs.scad>
use <../../core/proxy_interaction.scad>
use <../../core/segments.scad>

SC = 1;
IR = 20 * SC;

p = poly_antiprism(n = 7, p = 3, angle = 15);

// Keep the default view to one ordinary face. Add face 0 back when testing
// star-face punch-throughs, but do not make that the default baseline.
SHOW_FACES = [0]; //[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];

// Default to one known intersecting face so the carved result is visually
// obvious. `undef` means "all other faces".
CUTTER_FACE_INDICES = [15];
faces = poly_faces(p);
edges = _ps_edges_from_faces(faces);
CUTTER_VERTEX_INDICES = [];

// There are two distinct concepts here:
// - local edge-seat clearance, aligned to the target face boundary
// - foreign penetrating occupancy from other faces/edges/vertices
//
// This example now shows only the real composed path:
// clipped face occupancy minus clipped local clearance.
SHOW_CUTTER_FACES = true;
SHOW_LOCAL_EDGE_CLEARANCE = true;
SHOW_CUTTER_VERTICES = false;
SHOW_CARVED_RESULT = true;

function target_boundary_edge_indices(face_idx) =
    let(target_face = faces[face_idx])
        [
            for (k = [0 : 1 : len(target_face) - 1])
                ps_find_edge_index(edges, target_face[k], target_face[(k + 1) % len(target_face)])
        ];

EDGE_T = 3.5 * SC;
FACE_T = 1.6 * SC;
BASE_Z = -FACE_T / 4;
PILLOW_THK = 0.4 * SC;
VERTEX_PROXY_R = 2.2 * SC;
EDGE_STRIP_W = 1.2 * SC;
EDGE_STRIP_H = 100 * SC;
EDGE_INFLUENCE_LEN = undef;

module face_proxy() {
    let(
        face_loops = [
            for (seg = ps_face_segments([for (p = $ps_face_pts2d) [p[0], p[1], 0]], "nonzero", 1e-4))
                seg[0]
        ]
    )
    translate([0, 0, BASE_Z])
        linear_extrude(height = FACE_T + PILLOW_THK)
            union() {
                for (loop = face_loops)
                    ps_polygon(points = loop);
            }
}

module edge_proxy_strip() {
    #cube([$ps_edge_len, EDGE_STRIP_W, EDGE_STRIP_H], center = true);
}

module vertex_proxy() {
    sphere(r = VERTEX_PROXY_R, $fn = 24);
}

for (i = [0 : 1 : len(SHOW_FACES) - 1]) {
    fi = SHOW_FACES[i];
    color(["deepskyblue", "limegreen", "gold", "tomato", "plum"][i % 5])
        ps_carve_face_by_feature_proxies(
            p,
            fi,
            inter_radius = IR,
            face_bounds = [BASE_Z, BASE_Z + FACE_T + PILLOW_THK],
            face_proxy_mode = "sweep_to_bounds",
            edge_radius = EDGE_T / 2,
            edge_length = EDGE_INFLUENCE_LEN,
            vertex_radius = VERTEX_PROXY_R,
            include_local_edges = SHOW_LOCAL_EDGE_CLEARANCE,
            include_local_vertices = false,
            include_cutter_faces = SHOW_CUTTER_FACES,
            include_cutter_edges = false,
            include_cutter_vertices = SHOW_CUTTER_VERTICES,
            cutter_face_indices = CUTTER_FACE_INDICES,
            cutter_edge_indices = [],
            cutter_vertex_indices = CUTTER_VERTEX_INDICES,
            local_edge_indices = target_boundary_edge_indices(fi),
            local_vertex_indices = []
        ) {
            face_proxy();
            edge_proxy_strip();
            vertex_proxy();
            face_proxy();
            edge_proxy_strip();
            vertex_proxy();
        }
}
