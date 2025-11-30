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

module tetra_faces_sym(edge_len) {
    place_on_faces(tetrahedron, edge_len) children();
}

module tetra_faces_sym_ir(inter_radius) {
    place_on_faces_ir(tetrahedron, inter_radius) children();
}

// test demo
tetra_faces_sym_ir(40) {
    cylinder($fn = 3, r = $ph_edge_len / 2);
}