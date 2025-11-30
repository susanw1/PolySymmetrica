use <funcs.scad>
use <polygon-symm.scad>
use <edge-mount.scad>

// ---- Poly descriptor accessors ----
function poly_verts(poly)      = poly[0];
function poly_faces(poly)      = poly[1];
function poly_unit_edge(poly)  = poly[2];
function poly_e_over_ir(poly) = poly[3];

// ---- Generic face-frame helpers ----
function poly_face_center(poly, fi, scale) =
    let(f  = poly_faces(poly)[fi],
        vs = poly_verts(poly),
        v0 = vs[f[0]] * scale,
        v1 = vs[f[1]] * scale,
        v2 = vs[f[2]] * scale)
    (v0 + v1 + v2) / 3;

function poly_face_ez(poly, fi, scale) =
    let(f  = poly_faces(poly)[fi],
        vs = poly_verts(poly),
        v0 = vs[f[0]] * scale,
        v1 = vs[f[1]] * scale,
        v2 = vs[f[2]] * scale)
    v_norm(v_cross(v1 - v0, v2 - v0));   // outward normal


function poly_face_ex(poly, fi, scale) =
    let(f      = poly_faces(poly)[fi],
        vs     = poly_verts(poly),
        center = poly_face_center(poly, fi, scale),
        v0     = vs[f[0]] * scale)
    v_norm(v0 - center);   // local +X points to vertex 0

function poly_face_ey(poly, fi, scale) =
    v_cross(
        poly_face_ez(poly, fi, scale),
        poly_face_ex(poly, fi, scale)
    );

function frame_matrix(center, ex, ey, ez) = [
    [ex[0], ey[0], ez[0], center[0]],
    [ex[1], ey[1], ez[1], center[1]],
    [ex[2], ey[2], ez[2], center[2]],
    [0,      0,     0,     1]
];



// ---- Generic face-placement driver ----
module place_on_faces(poly, edge_len) {
    scale = edge_len / poly_unit_edge(poly);

    for (fi = [0 : len(poly_faces(poly))-1]) {
        center = poly_face_center(poly, fi, scale);
        ex     = poly_face_ex(poly, fi, scale);
        ey     = poly_face_ey(poly, fi, scale);
        ez     = poly_face_ez(poly, fi, scale);

        $ph_facet_idx = fi;
        $ph_edge_len  = edge_len;

        multmatrix(frame_matrix(center, ex, ey, ez))
            children();
    }
}

// ---- Generic inter-radius driver ----
module place_on_faces_ir(poly, inter_radius) {
    edge_len = inter_radius * poly_e_over_ir(poly);
    place_on_faces(poly, edge_len) children();
}


// ---- Canonical Octahedron (edge length = sqrt(2)) ----
octahedron = [
    // verts (index 0)
    [
        [ 1, 0, 0],  // 0
        [-1, 0, 0],  // 1
        [ 0, 1, 0],  // 2
        [ 0,-1, 0],  // 3
        [ 0, 0, 1],  // 4 (top)
        [ 0, 0,-1]   // 5 (bottom)
    ],

    // faces (index 1) – triangles by vertex index
    [
        [0,2,4], [2,1,4], [1,3,4], [3,0,4],
        [2,0,5], [1,2,5], [3,1,5], [0,3,5]
    ],

    // unit_edge (index 2)
    sqrt(2),

    // e_over_ir (index 3): edge_len / inter_radius
    2
];


module octa_faces_sym(edge_len) {
    place_on_faces(octahedron, edge_len) children();
}

module octa_faces_sym_ir(inter_radius) {
    place_on_faces_ir(octahedron, inter_radius) children();
}

//// edge-based
//octa_faces_sym(40) {
//    regular_polygon_frame(3, $ph_edge_len);
//}

// inter-radius-based (Cundy/Rollett style)
octa_faces_sym_ir(30) {
    regular_polygon_frame(3, $ph_edge_len);
}
