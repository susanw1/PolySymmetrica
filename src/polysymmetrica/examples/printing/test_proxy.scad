use <../../core/prisms.scad>
use <../../core/funcs.scad>
use <../../core/proxy_interaction.scad>
use <../../core/render.scad>
use <../../core/segments.scad>

SC = 1;
IR = 20 * SC;

p = poly_antiprism(n = 7, p = 3, angle = 15);

// Keep the default view to one ordinary face. Add face 0 back when testing
// star-face punch-throughs, but do not make that the default baseline.
SHOW_FACES = [0,2,4]; //[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];

// `undef` means "all other faces".
FACE_PROXY_INDICES = undef;
faces = poly_faces(p);
edges = _ps_edges_from_faces(faces);
VERTEX_PROXY_INDICES = [];

// Baseline: show the edge-space subtraction alone. Turn face proxies back on
// once the seat itself is visually understood.
SHOW_FACE_PROXIES = false;
SHOW_EDGE_PROXIES = true;
SHOW_VERTEX_PROXIES = false;

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
EDGE_STRIP_PAD = 0;
EDGE_INFLUENCE_LEN = undef;

function edge_strip_height(zMax, strip_w, dihedral, eps = 1e-9) =
    let(
        alpha = (180 - dihedral) / 2,
        y_proj = abs(sin(alpha)),
        z_proj = abs(cos(alpha)),
        rem = max(0, 2 * zMax - strip_w * y_proj)
    )
    (z_proj <= eps)
        ? 2 * zMax
        : rem / z_proj;

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
    strip_h = edge_strip_height(FACE_T + PILLOW_THK, EDGE_STRIP_W, $ps_dihedral);
    #cube([$ps_edge_len, EDGE_STRIP_W, strip_h], center = true);
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
            edge_length = EDGE_INFLUENCE_LEN,
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
