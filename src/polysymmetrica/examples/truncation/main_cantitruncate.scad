use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>
use <util_demo.scad>

/**
Demonstrates cantitruncation (two-parameter). Try varying t (face shift) and c (edge/vertex expansion).
*/

for (p = with_index([
            undef, // adds text to show 't' value
            tetrahedron(), hexahedron(), octahedron(), dodecahedron(), icosahedron()
        ]),
     t = with_index([0.1, 0.25, 0.4], 0),
     c = with_index([0.01, 0.1, 0.2], 0)) {
    translate([120 * c[0], 120 * p[0], -120 * t[0]])
        if (is_undef(p[1])) 
            text(str("t=", t[1]));
        else
            demo(poly_cantitruncate(p[1], t[1], c[1]));
}

// Single example
//demo(poly_cantitruncate(hexahedron(), 0.2, 0.1));
//demo(poly_cantitruncate(tetrahedron(), 0.1, 0.2));

// Uniform solver example (slow)
