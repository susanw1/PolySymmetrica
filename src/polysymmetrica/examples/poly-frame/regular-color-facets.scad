use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>

use <../../core/duals.scad>


COLORS = [ "black", "red", "green", "yellow", "blue", "magenta", "cyan", "white", "orange" ];
module color_map(L, i) { color(COLORS[L[i]]) children(); }

// triangular-faced regular polyhedra 

translate([-100, 0, 100])
place_on_faces_ir(tetrahedron(), 30) {
    color_map([1, 4, 2, 3], $ps_facet_idx)
        circle(r = $ps_facet_radius, $fn = 3);
}

translate([0, 0, 100])
place_on_faces_ir(octahedron(), 30) {
    color_map([1, 4, 2, 3, 2, 3, 1, 4], $ps_facet_idx)
        circle(r = $ps_facet_radius, $fn = 3);
}


translate([100, 0, 100])
place_on_faces_ir(icosahedron(), 30) {
    color_map([1, 4, 2, 3, 8, 3, 8, 2, 3, 2,    1, 8, 4, 4, 1, 8, 4, 3, 2, 1 ], $ps_facet_idx)
        circle(r = $ps_facet_radius, $fn = 3);
}

// duals

translate([0, 100, 100])
place_on_faces_ir(poly_dual(octahedron()), 30) {
    color_map([1, 8, 4, 3, 7, 2 ], $ps_facet_idx)
        circle(r = $ps_facet_radius, $fn = 4);
}

translate([100, 100, 100])
place_on_faces_ir(poly_dual(icosahedron()), 30) {
    color_map([ 1, 4, 3, 4, 4, 2, 3, 2, 2, 3, 1, 1 ], $ps_facet_idx)
        circle(r = $ps_facet_radius, $fn = 5);
}
