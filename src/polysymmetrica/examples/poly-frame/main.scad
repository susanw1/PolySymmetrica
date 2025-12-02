//use <edge-mount.scad>
//use <../../models/tetrahedron.scad>
//use <../../models/octahedron.scad>
//use <../../models/icosahedron.scad>
//
//
////translate([-100, 0, 0])
////tetra_faces_sym_ir(30) {
////    regular_polygon_frame(3, $ps_edge_len);
////}
////
////translate([0, 0, 0])
////octa_faces_sym_ir(30) {
////    regular_polygon_frame(3, $ps_edge_len);
////}
////
////
////translate([100, 0, 0])
////icosa_faces_sym_ir(30) {
////    regular_polygon_frame(3, $ps_edge_len);
////}
//
//difference() {
//    N = 3;
//    tetra_faces_sym_ir(30) {
//        regular_polygon_frame(N, $ps_edge_len);
//    }
//    tetra_faces_sym_ir(30) {
//        // main facet base
//        translate([0, 0, 1]) linear_extrude(height = 5) {
//            regular_polygon_2d(N, $ps_edge_len);
//        }
//    }
//}
//


use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>
use <edge-mount.scad>


//place_on_faces_ir(octahedron(), 30)
//    circle(r = $ps_facet_radius, $fn = 3);
    
place_on_edges_ir(icosahedron(), 40)
//    cube([$ps_edge_len, 1, 1], center = true);
    translate([-$ps_edge_len/2, 0, 0]) edge_mount($ps_edge_len);    