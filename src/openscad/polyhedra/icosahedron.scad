use <poly-core.scad>

phi = (1 + sqrt(5)) / 2;

// ---- Canonical Icosahedron: edge length = 2 ----
icosahedron = [

    // verts (index 0)
    [
        // (0, ±1, ±φ)
        [ 0,  1,  phi],   // 0
        [ 0, -1,  phi],   // 1
        [ 0,  1, -phi],   // 2
        [ 0, -1, -phi],   // 3

        // (±1, ±φ, 0)
        [ 1,  phi,  0],   // 4
        [-1,  phi,  0],   // 5
        [ 1, -phi,  0],   // 6
        [-1, -phi,  0],   // 7

        // (±φ, 0, ±1)
        [ phi,  0,  1],   // 8
        [-phi,  0,  1],   // 9
        [ phi,  0, -1],   // 10
        [-phi,  0, -1]    // 11
    ],

    // faces (index 1) – oriented so normals point OUTWARDS
    [
        [0, 1,  8],
        [0, 9,  1],
        [0, 4,  5],
        [0, 8,  4],
        [0, 5,  9],

        [1, 7,  6],
        [1, 6,  8],
        [1, 9,  7],

        [2,10,  3],
        [2, 3, 11],
        [2, 5,  4],
        [2, 4, 10],
        [2,11,  5],

        [3, 6,  7],
        [3,10,  6],
        [3, 7, 11],

        [4, 8, 10],
        [5,11,  9],
        [6,10,  8],
        [7, 9, 11]
    ],

    // unit_edge (index 2)
    2,

    // e_over_ir (index 3): edge_len / inter_radius
    2 / phi
];


// Icosahedral faces (edge-based)
module icosa_faces_sym(edge_len) {
    place_on_faces(icosahedron, edge_len) children();
}


// Icosahedral faces (inter-radius-based)
module icosa_faces_sym_ir(inter_radius) {
    place_on_faces_ir(icosahedron, inter_radius) children();
}


// Icosahedral vertices (edge-based)
module icosa_vertices_sym(edge_len) {
    place_on_vertices(icosahedron, edge_len) children();
}

// Icosahedral vertices (inter-radius-based)
module icosa_vertices_sym_ir(inter_radius) {
    place_on_vertices_ir(icosahedron, inter_radius) children();
}

//////
// TEST DEMOS
//////

icosa_faces_sym_ir(40) {
    cylinder($fn = 3, r = $ph_facet_radius );

    color("white") translate([0,0,1])
        text(str($ph_facet_idx), size = 5, halign="center", valign="center");

    // starts at face centre, points to poly centre
    color("green")
        translate([0,0,2-norm($ph_poly_center_local)])
        cylinder(h = norm($ph_poly_center_local), r = 0.5, center = false);
}

icosa_vertices_sym_ir(40) {
    color("blue") cylinder($fn = 5, r = $ph_vert_radius/2);
    color("red") sphere(5);
    cylinder(h = 20, r = 2, center = false);  // along local +Z
}
