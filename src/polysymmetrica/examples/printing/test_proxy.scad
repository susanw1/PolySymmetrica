use <../../core/prisms.scad>
use <../../core/funcs.scad>
use <../../core/proxy_interaction.scad>
use <../../core/render.scad>
use <../../core/segments.scad>

SC = 1;
IR = 20 * SC;

p = poly_antiprism(n = 7, p = 3, angle = 15);

SHOW_FACES = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];

// Optional exact neighboring face proxies to instantiate around each target.
// `undef` means "all other faces".
FACE_PROXY_INDICES = undef;
faces = poly_faces(p);
edges = _ps_edges_from_faces(faces);
VERTEX_PROXY_INDICES = [];

SHOW_FACE_PROXIES = true;
SHOW_EDGE_PROXIES = true;
SHOW_VERTEX_PROXIES = len(VERTEX_PROXY_INDICES) > 0;

function target_edge_indices(face_idx) =
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
EDGE_STRIP_H = FACE_T + PILLOW_THK + 1.0 * SC;
EDGE_STRIP_PAD = 6 * SC;

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
    cube([$ps_edge_len + EDGE_STRIP_PAD, EDGE_STRIP_W, EDGE_STRIP_H], center = true);
}

module vertex_proxy() {
    sphere(r = VERTEX_PROXY_R, $fn = 24);
}

for (i = [0 : 1 : len(SHOW_FACES) - 1]) {
    fi = SHOW_FACES[i];
    color(["deepskyblue", "limegreen", "gold", "tomato", "plum"][i % 5])
    ps_clip_face_by_feature_proxies(
        p,
        fi,
        inter_radius = IR,
        face_bounds = [BASE_Z, BASE_Z + FACE_T + PILLOW_THK],
        face_proxy_mode = "sweep_to_bounds",
        edge_radius = EDGE_T / 2,
        edge_length = IR,
        vertex_radius = VERTEX_PROXY_R,
        include_faces = SHOW_FACE_PROXIES,
        include_edges = SHOW_EDGE_PROXIES,
        include_vertices = SHOW_VERTEX_PROXIES,
        face_indices = FACE_PROXY_INDICES,
        edge_indices = target_edge_indices(fi),
        vertex_indices = VERTEX_PROXY_INDICES
    ) {
        face_proxy();
        edge_proxy_strip();
        vertex_proxy();
    }
}
