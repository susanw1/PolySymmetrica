use <poly-core.scad>

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


// test demo
octa_faces_sym_ir(40) {
    cylinder($fn = 3, r = $ph_edge_len / 2);
}