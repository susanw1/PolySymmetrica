use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/platonics_all.scad>

use <util_demo.scad>

/**
Demonstrates chamfering.
*/
for (   p = with_index([
            undef, // adds text to show 't' value
            tetrahedron(), hexahedron(), octahedron(), dodecahedron(), icosahedron(),
            poly_truncate(tetrahedron()), poly_dual(poly_rectify(icosahedron()))
            ]),
        t = with_index([ -0.4, -0.1, 0, 0.1, 0.2, 0.4, 1.01, 1.5, 2 ], -1)) {
    translate([100 * t[0], 100 * p[0], 0]) 
        if (is_undef(p[1])) 
            text(str("t=", t[1]));
        else
            demo(poly_chamfer(p[1], t[1]));
}        
