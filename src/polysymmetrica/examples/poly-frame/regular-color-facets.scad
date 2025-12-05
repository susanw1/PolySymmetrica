use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>

use <../../core/duals.scad>


// triangular-faced regular polyhedra 

translate([-100, 0, 100])
place_on_faces_ir(tetrahedron(), 30) {
    color(["red", "blue", "green", "yellow"][$ps_facet_idx])
        circle(r = $ps_facet_radius, $fn = 3);
}

translate([0, 0, 100])
place_on_faces_ir(octahedron(), 30) {
    color(["red", "blue", "green", "yellow", "green", "yellow", "red", "blue"][$ps_facet_idx])
        circle(r = $ps_facet_radius, $fn = 3);
}


translate([100, 0, 100])
place_on_faces_ir(icosahedron(), 30) {
    color(["red", "blue", "green", "yellow", "orange", "yellow", "orange", "green","yellow","green",
        "red", "orange", "blue", "blue", "red", "orange", "blue", "yellow", "green", "red"][$ps_facet_idx])
        circle(r = $ps_facet_radius, $fn = 3);
}

// duals

translate([0, 100, 100])
place_on_faces_ir(poly_dual(octahedron()), 30) {
    color(["red", "orange", "blue", "yellow", "white", "green"][$ps_facet_idx])
        circle(r = $ps_facet_radius, $fn = 4);
}

translate([100, 100, 100])
place_on_faces_ir(poly_dual(icosahedron()), 30) {
    color(["red", "blue", "yellow", "blue", "blue", "green", "yellow", "green", "green", "yellow", "red", "red"][$ps_facet_idx])
        circle(r = $ps_facet_radius, $fn = 5);
}
