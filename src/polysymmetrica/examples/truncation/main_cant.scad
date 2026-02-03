use <../../core/truncation.scad>
use <../../core/cantellation.scad>

use <../../models/platonics_all.scad>

use <util_demo.scad>


/**
Demonstrates extended cantellations with extreme values of 't'.
*/
for (   p = [ [-1, undef],
            [0, tetrahedron()] , [1, hexahedron()], [2, octahedron()], [3, dodecahedron()], [4, icosahedron()]],
        t = [[-3, -1],[ -2, -0.7], [-1, -0.2], [0, 0.001], [1, 0.2], [2, 0.4], [3, 0.6], [4, 1], [5, 2.5]]) {

    translate([100 * t[0], 100 * p[0], 0])
        if (is_undef(p[1]))
            text(str("t=", t[1]));
        else
            demo(poly_cantellate(p[1], t[1]));
}
