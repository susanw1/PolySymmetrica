// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <../core/funcs.scad>
use <../core/placement.scad>

phi = (1 + sqrt(5)) / 2;

// ---- Canonical Icosahedron: edge length = 2 ----
function icosahedron() = let(unit_edge = 2) make_poly(
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
    ] / unit_edge,

    // faces (index 1) – oriented for LHR (OpenSCAD)
    [
        [0,  8,  1],
        [0,  1,  9],
        [0,  5,  4],
        [0,  4,  8],
        [0,  9,  5],

        [1,  6,  7],
        [1,  8,  6],
        [1,  7,  9],

        [2,  3, 10],
        [2, 11,  3],
        [2,  4,  5],
        [2, 10,  4],
        [2,  5, 11],

        [3,  7,  6],
        [3,  6, 10],
        [3, 11,  7],

        [4, 10,  8],
        [5,  9, 11],
        [6,  8, 10],
        [7, 11,  9]
    ],

    // e_over_ir (index 3): edge_len / inter_radius
    2 / phi
);


////////////////////////////////////////////////////////////////////////
// TEST DEMOS
//////

place_on_faces(icosahedron(), 40) {
    cylinder($fn = 3, r = $ps_facet_radius );

    color("white") translate([0,0,1])
        text(str($ps_facet_idx), size = 5, halign="center", valign="center");

    // starts at face centre, points to poly centre
    color("green")
        translate([0,0,2-norm($ps_poly_center_local)])
        cylinder(h = norm($ps_poly_center_local), r = 0.5, center = false);
}

place_on_vertices(icosahedron(), 40) {
    color("blue") cylinder($fn = 5, r = $ps_vert_radius/2);
    color("red") sphere(5);
    cylinder(h = 10, r = 2, center = false);  // along local +Z
}

place_on_faces(icosahedron(), 40) {
    cylinder($fn = 3, r = $ps_facet_radius);
}

place_on_vertices(icosahedron(), 40) {
    color("blue") cylinder($fn = 5, r = $ps_vert_radius/3);
    color("pink") translate([0,0,10])
        text(str($ps_vertex_idx), size = 5, halign="center", valign="center");
}

place_on_edges(icosahedron(), 50) {
    color("black") sphere(r=3);
    color("gray") translate([0,0,3])
        text(str($ps_edge_idx), size = 5, halign="center", valign="center");
}
