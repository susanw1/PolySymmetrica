use <poly-core.scad>

// ---- Canonical Octahedron (edge length = sqrt(2)) ----
function octahedron() = [
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


// Octahedral faces (edge-based)
module octa_faces_sym(edge_len) {
    place_on_faces(octahedron(), edge_len) children();
}

// Octahedral faces (inter-radius-based)
module octa_faces_sym_ir(inter_radius) {
    place_on_faces_ir(octahedron(), inter_radius) children();
}

// Octahedral vertices (edge-based)
module octa_vertices_sym(edge_len) {
    place_on_vertices(octahedron(), edge_len) children();
}

// Octahedral vertices (inter-radius-based)
module octa_vertices_sym_ir(inter_radius) {
    place_on_vertices_ir(octahedron(), inter_radius) children();
}

// Octahedral edges (edge-based)
module octa_edges_sym(edge_len) {
    place_on_edges(octahedron(), edge_len) children();
}
// Octahedral edges (inter-radius-based)
module octa_edges_sym_ir(inter_radius) {
    place_on_edges_ir(octahedron(), inter_radius) children();
}


//////
// TEST DEMOS
//////

octa_faces_sym_ir(40) {
    color("yellow") cylinder($fn = 3, r = $ph_facet_radius / 2);
    color("green")
        translate([0,0,2-norm($ph_poly_center_local)])
        cylinder(h = norm($ph_poly_center_local), r = 0.5, center = false);
}

octa_vertices_sym(40) {
    color("red") sphere(5);
    cylinder(h = 20, r = 2, center = false);  // along local +Z
}

octa_vertices_sym(40) {
    // Draw an arrow along +X to show incident-edge direction
    color("cyan")
    cube([8,1,1], center=false);
}


octa_faces_sym(40) {
    face_debug();
}

octa_edges_sym_ir(40) {
    color("black") cube([$ph_edge_len/2,1,1], center=true);
}
