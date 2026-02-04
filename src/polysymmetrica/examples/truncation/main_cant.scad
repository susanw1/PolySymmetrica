use <../../core/truncation.scad>

use <../../models/platonics_all.scad>

use <util_demo.scad>


/**
Demonstrates extended cantellations with extreme values of 't'.
*/
for (   p = with_index([ undef, 
                tetrahedron(), hexahedron(), octahedron(), dodecahedron(), icosahedron() 
            ], -1),
        t = with_index([ -1, -0.7, -0.2,  0.001, 0.2, 0.4, 0.6, 1, 2.5], -3)) {

    translate([100 * t[0], 100 * p[0], 0])
        if (is_undef(p[1]))
            linear_extrude(1) text(str("t=", t[1]));
        else
            demo(poly_cantellate(p[1], t[1]));
}
