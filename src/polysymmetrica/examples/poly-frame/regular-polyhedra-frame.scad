use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>
use <edge-mount.scad>

use <../../core/duals.scad>


// triangular-faced regular polyhedra 
translate([-100, 0, -100]) 
color("grey")
place_on_faces_ir(tetrahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}

translate([0, 0, -100])
color("grey")
place_on_faces_ir(octahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}


translate([100, 0, -100])
color("grey")
place_on_faces_ir(icosahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}

// duals

translate([0, 100, -100])
color("grey")
place_on_faces_ir(poly_dual(octahedron()), 30) {
    regular_polygon_frame(4, $ps_edge_len);
}

translate([100, 100, -100])
color("grey")
place_on_faces_ir(poly_dual(icosahedron()), 30) {
    regular_polygon_frame(5, $ps_edge_len);
}
