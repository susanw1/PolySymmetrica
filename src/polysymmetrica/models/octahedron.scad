// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <../core/placement.scad>

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


//////
// TEST DEMOS
//////

place_on_faces_ir(octahedron(), 40) {
    color("yellow") cylinder($fn = 3, r = $ps_facet_radius / 4);
    color("green")
        translate([0,0,1-norm($ps_poly_center_local)])
        cylinder(h = norm($ps_poly_center_local), r = 0.5, center = false);
}

place_on_vertices(octahedron(), 40) {
    color("red") sphere(5);
    cylinder(h = 20, r = 2, center = false);  // along local +Z
}

place_on_vertices(octahedron(), 40) {
    // Draw an arrow along +X to show incident-edge direction
    color("cyan")
    cube([8,1,1], center=false);
}


place_on_faces(octahedron(), 40) {
    face_debug();
}

place_on_edges_ir(octahedron(), 40) {
    color("black") cube([$ps_edge_len/2,1,1], center=true);
}
