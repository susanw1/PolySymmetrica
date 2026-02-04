use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/platonics_all.scad>

use <util_demo.scad>

/**
Demonstrates extended truncations with extreme values of 't'. A value of t=1/3 gives normal regular truncation 
for triangular faces (eg truncated tetrahedron), and t=1/2 produces rectification (eg cuboctahedron).

But values outside the range 0 < t < 0.5 cause the new vertex faces to over-extend.
* t < 0: anti-truncation
* 0.5 < t <= 1: hyper-truncation
* t > 1: quasi-truncation
*/
for (   p = with_index([
                undef,
                tetrahedron(), hexahedron(), octahedron(), dodecahedron(), icosahedron()
            ], -1),
        t = with_index([ -0.3, 0.01, 0.2, 0.45, 0.55, 0.7, 0.9, 1.0, 1.5, 2.5 ], -1)) {
    translate([100 * t[0], 100 * p[0], 0]) 
        if (is_undef(p[1])) 
            linear_extrude(1) text(str("t=", t[1]));
        else
            demo(poly_truncate(p[1], t[1]));
}        

