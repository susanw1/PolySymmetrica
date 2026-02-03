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
for (   p = [ [-1, undef],
            [0, tetrahedron()] , [1, hexahedron()], [2, octahedron()], [3, dodecahedron()], [4, icosahedron()]], 
        t = [[-1, -0.3], [0, 0.01], [1, 0.2], [2, 0.45], [3, 0.55], 
            [4, 0.7], [5, 0.9], [6, 1.0], [7, 1.5], [8, 2.5]]) {
    translate([100 * t[0], 100 * p[0], 0]) 
        if (is_undef(p[1])) 
            text(str("t=", t[1]));
        else
            demo(poly_truncate(p[1], t[1]));
}        


