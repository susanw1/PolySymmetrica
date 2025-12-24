// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>
use <edge-mount.scad>

use <../../core/duals.scad>


// triangular-faced regular polyhedra
translate([-100, 0, -100])
color("grey")
place_on_faces(tetrahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}

translate([0, 0, -100])
color("grey")
place_on_faces(octahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}


translate([100, 0, -100])
color("grey")
place_on_faces(icosahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}

// duals

translate([0, 100, -100])
color("grey")
place_on_faces(poly_dual(octahedron()), 30) {
    regular_polygon_frame(4, $ps_edge_len);
}

translate([100, 100, -100])
color("grey")
place_on_faces(poly_dual(icosahedron()), 30) {
    regular_polygon_frame(5, $ps_edge_len);
}
