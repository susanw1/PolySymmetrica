use <poly-core.scad>

// ---- Canonical Tetrahedron (edge length = 2*sqrt(2)) ----
tetrahedron = [
    // verts (index 0)
    [
        [ 1,  1,  1],   // 0
        [-1, -1,  1],   // 1
        [-1,  1, -1],   // 2
        [ 1, -1, -1]    // 3
    ],

    // faces (index 1) – triangles by vertex index, oriented consistently
    [
        [0, 2, 1],
        [0, 1, 3],
        [0, 3, 2],
        [1, 2, 3]
    ],
    
    // unit_edge (index 2): canonical edge length
    2 * sqrt(2),

    // e_over_ir (index 3): edge_len / inter_radius
    2 * sqrt(2)
];

// Tetrahedral faces (edge-based)
module tetra_faces_sym(edge_len) {
    place_on_faces(tetrahedron, edge_len) children();
}

// Tetrahedral faces (inter-radius-based)
module tetra_faces_sym_ir(inter_radius) {
    place_on_faces_ir(tetrahedron, inter_radius) children();
}

// Tetrahedral vertices (edge-based)
module tetra_vertices_sym(edge_len) {
    place_on_vertices(tetrahedron, edge_len) children();
}

// Tetrahedral vertices (inter-radius-based)
module tetra_vertices_sym_ir(inter_radius) {
    place_on_vertices_ir(tetrahedron, inter_radius) children();
}

// Tetrahedral edges (edge-based)
module tetra_edges_sym(edge_len) {
    place_on_edges(tetrahedron, edge_len) children();
}
// Tetrahedral edges (inter-radius-based)
module tetra_edges_sym_ir(inter_radius) {
    place_on_edges_ir(tetrahedron, inter_radius) children();
}

//////
// TEST DEMOS
//////

tetra_faces_sym(40) {
    cylinder($fn = 3, r = $ph_facet_radius);
}

tetra_vertices_sym(40) {
    color("red") sphere(5);
    color("blue") cylinder(h = 10, r = 1, center = false);  // along local +Z
}

