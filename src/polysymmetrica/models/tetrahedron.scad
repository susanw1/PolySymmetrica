// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <../core/funcs.scad>
use <../core/placement.scad>

// ---- Canonical Tetrahedron (edge length = 2*sqrt(2)) ----
function tetrahedron() = let(unit_edge = 2 * sqrt(2)) make_poly(
    // verts (index 0)
    [
        [ 1,  1,  1],   // 0
        [-1, -1,  1],   // 1
        [-1,  1, -1],   // 2
        [ 1, -1, -1]    // 3
    ]/ unit_edge,

    // faces (index 1) – triangles by vertex index, oriented for LHR (OpenSCAD)
    [
        [0, 1, 2],
        [0, 3, 1],
        [0, 2, 3],
        [1, 3, 2]
    ],

    // e_over_ir (index 3): edge_len / inter_radius
    2 * sqrt(2)
);


//////
// TEST DEMOS
//////

place_on_faces(tetrahedron(), 40) {
    cylinder($fn = 3, r = $ps_face_radius);
}

place_on_vertices(tetrahedron(), 40) {
    color("red") sphere(5);
    color("blue") cylinder(h = 10, r = 1, center = false);  // along local +Z
}
