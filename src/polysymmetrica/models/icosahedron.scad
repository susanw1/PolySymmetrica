// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <../core/funcs.scad>
use <../core/placement.scad>

phi = (1 + sqrt(5)) / 2;

// ---- Canonical Icosahedron: edge length = 2 ----
function icosahedron() = make_poly(
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
);


////////////////////////////////////////////////////////////////////////
// TEST DEMOS
//////

place_on_faces_ir(icosahedron(), 40) {
    cylinder($fn = 3, r = $ps_facet_radius );

    color("white") translate([0,0,1])
        text(str($ps_facet_idx), size = 5, halign="center", valign="center");

    // starts at face centre, points to poly centre
    color("green")
        translate([0,0,2-norm($ps_poly_center_local)])
        cylinder(h = norm($ps_poly_center_local), r = 0.5, center = false);
}

place_on_vertices_ir(icosahedron(), 40) {
    color("blue") cylinder($fn = 5, r = $ps_vert_radius/2);
    color("red") sphere(5);
    cylinder(h = 10, r = 2, center = false);  // along local +Z
}

place_on_faces_ir(icosahedron(), 40) {
    cylinder($fn = 3, r = $ps_facet_radius);
}

place_on_vertices_ir(icosahedron(), 40) {
    color("blue") cylinder($fn = 5, r = $ps_vert_radius/3);
    color("pink") translate([0,0,10])
        text(str($ps_vertex_idx), size = 5, halign="center", valign="center");
}

place_on_edges_ir(icosahedron(), 50) {
    color("black") sphere(r=3);
    color("gray") translate([0,0,3])
        text(str($ps_edge_idx), size = 5, halign="center", valign="center");
}
