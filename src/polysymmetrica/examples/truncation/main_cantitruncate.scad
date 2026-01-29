use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>
use <util_demo.scad>

/**
Demonstrates cantitruncation (two-parameter). Try varying t (face shift) and c (edge/vertex expansion).
*/

for (p = with_index([
            tetrahedron(), hexahedron(), octahedron(), dodecahedron(), icosahedron()
        ]),
     t = with_index([0.1, 0.2, 0.3], -1),
     c = with_index([0.0, 0.1, 0.2], -1)) {
    translate([120 * c[0], 120 * p[0], -120 * t[0]])
        demo(poly_cantitruncate(p[1], t[1], c[1]));
}

// Single example
//demo(poly_cantitruncate(hexahedron(), 0.2, 0.1));
