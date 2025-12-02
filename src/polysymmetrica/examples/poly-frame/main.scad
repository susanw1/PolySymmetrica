use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>
use <edge-mount.scad>

use <../../core/duals.scad>


// triangular-faced regular polyhedra 
translate([-100, 0, 0])
place_on_faces_ir(tetrahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}

translate([0, 0, 0])
place_on_faces_ir(octahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}


translate([100, 0, 0])
place_on_faces_ir(icosahedron(), 30) {
    regular_polygon_frame(3, $ps_edge_len);
}

// duals

translate([0, 100, 0])
place_on_faces_ir(poly_dual(octahedron()), 30) {
    regular_polygon_frame(4, $ps_edge_len);
}

translate([100, 100, 0])
place_on_faces_ir(poly_dual(icosahedron()), 30) {
    regular_polygon_frame(5, $ps_edge_len);
}


//difference() {
//    N = 3;
//    place_on_faces_ir(tetrahedron(), 30) {
//        regular_polygon_frame(N, $ps_edge_len);
//    }
//    place_on_faces_ir(tetrahedron(), 30) {
//        // main facet base
//        translate([0, 0, 1]) linear_extrude(height = 5) {
//            regular_polygon_2d(N, $ps_edge_len);
//        }
//    }
//}



//// examples in dev-guide

//use <../../core/placement.scad>
//use <../../core/duals.scad>
//use <../../models/tetrahedron.scad>
//use <../../models/octahedron.scad>
//use <../../models/icosahedron.scad>
//use <edge-mount.scad>
//
//
////place_on_faces_ir(octahedron(), 30)
////    circle(r = $ps_facet_radius, $fn = 3);
//    
//place_on_edges_ir(poly_dual(icosahedron()), 40)
//    cube([$ps_edge_len, 1, 1], center = true);
////    translate([-$ps_edge_len/2, 0, 0]) edge_mount($ps_edge_len);    