use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>

use <../../core/duals.scad>
use <../../core/truncation.scad>

//// examples in README / dev-guide

// 5.1 Octahedron face mounts
translate([-100,0,0])
place_on_faces_ir(octahedron(), 30)
    circle(r = $ps_facet_radius, $fn = $ps_vertex_count);
    
// 5.2 Edge frames on an icosahedron
place_on_edges_ir(icosahedron(), 40)
    cube([$ps_edge_len, 5, 1], center = true);


// 5.3 Dodecahedron from dual of icosa
translate([200,0,0])
place_on_faces_ir(poly_dual(icosahedron()), 40)
    translate([0,0, -$ps_face_midradius + 20]) 
        cylinder(r1 = 0, r2 = $ps_facet_radius, h = $ps_face_midradius, $fn = $ps_vertex_count);

// 5.4 Truncated icosa
translate([350,0,0])
place_on_faces_ir(poly_truncate(icosahedron()), 40)
    color($ps_vertex_count == 5? "blue" : "orange")
        cylinder(r = $ps_facet_radius, $fn = $ps_vertex_count, h = 0.2);


// 5.5 Catalan from dual of truncated octahedron
translate([500,0,0])
place_on_faces_ir(poly_dual(poly_truncate(icosahedron())), 40)
    linear_extrude(height = 0.2) polygon(points = $ps_face_pts2d);

// 5.6 Visual debugging: original + dual overlay
translate([700,0,0]) {
    color("red", alpha = 0.5)
    place_on_faces_ir(octahedron(), 30)
        cylinder(r = $ps_facet_radius, $fn = $ps_vertex_count, h = 0.2);

    color("gold", alpha = 0.3)
    place_on_faces_ir(poly_dual(octahedron()), 30)    // note scaling for vertex/face alignment
        cylinder(r = $ps_facet_radius, $fn = $ps_vertex_count, h = 0.2);
}