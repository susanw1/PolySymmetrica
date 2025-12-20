use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>

use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>

use <../../core/render.scad>

IR = 100;

//p = poly_truncate(octahedron(), 0.3);

//p1 = poly_truncate(icosahedron(), 0.45);
//p = poly_truncate(p1, 0.3);

//p = poly_truncate(poly_truncate(poly_truncate(octahedron(), 1/3), 1/3), 1/3);
//p = poly_dual(poly_truncate(icosahedron(), 0.4999));
//p = tetrahedron();
//p = icosahedron();


place_on_faces_ir(p, IR) {
    let (col = $ps_vertex_count == 4? "red" : "blue") {
        color(col) {
            linear_extrude(height = 1) polygon(points = $ps_face_pts2d);
//            translate([0,0,5]) text(str($ps_facet_idx));
        }
    }
}

place_on_edges_ir(p, IR) {
    color("silver") cube([$ps_edge_len, 2, 2], center = true);
//    translate($ps_poly_center_local) cube([1,1,$ps_edge_midradius], center = false);
//    translate([0,0,5]) text(str($ps_edge_idx));
}

//render_poly(p, 30);

place_on_vertices_ir(p, IR) {
    color("gold") sphere(2, $fn=30);
//    translate([0,0,5]) text(str($ps_vertex_idx));
}