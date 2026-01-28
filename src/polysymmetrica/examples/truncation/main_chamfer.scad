use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <util_demo.scad>

/**
Demonstrates chamfering.
*/
*for (   p = [ [-1, undef],
            [0, tetrahedron()]], 
        t = [[-1, -0.3], [0, 0.01], [1, 0.2], [2, 0.45], [3, 0.55], 
            [4, 0.7], [5, 0.9], [6, 1.0], [7, 1.5], [8, 2.5]]) {
    translate([100 * t[0], 100 * p[0], 0]) 
        if (is_undef(p[1])) 
            text(str("t=", t[1]));
        else
            demo(poly_chamfer(p[1], t[1]));
}        

p1 = poly_dual(poly_rectify(dodecahedron()));
p = poly_chamfer(p1, 0.05);
demo(p);
//poly_describe(p1);
//poly_render(p, 20);
