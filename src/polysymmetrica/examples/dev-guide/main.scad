use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>

use <../../core/duals.scad>

//// examples in README / dev-guide

translate([-100,0,0])
place_on_faces_ir(octahedron(), 30)
    circle(r = $ps_facet_radius, $fn = 3);
    

place_on_edges_ir(icosahedron(), 40)
    cube([$ps_edge_len, 5, 1], center = true);


translate([200,0,0])
place_on_faces_ir(poly_dual(icosahedron()), 40)
        translate([0,0, -$ps_face_midradius + 10]) 
            cylinder(r1 = 0, r2 = $ps_facet_radius, h = $ps_face_midradius, $fn=5);

translate([400,0,0]) {
    color("green")
    place_on_faces_ir(octahedron(), 30)
        circle(r = $ps_facet_radius, $fn = 3);

    color("gold", alpha = 0.2)
    place_on_faces_ir(poly_dual(octahedron()), 40)    // note scaling for vertex/face alignment
        circle(r = $ps_facet_radius, $fn = 4);
}