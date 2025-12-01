use <edge-mount.scad>
use <tetrahedron.scad>
use <octahedron.scad>
use <icosahedron.scad>


//translate([-100, 0, 0])
//tetra_faces_sym_ir(30) {
//    regular_polygon_frame(3, $ph_edge_len);
//}
//
//translate([0, 0, 0])
//octa_faces_sym_ir(30) {
//    regular_polygon_frame(3, $ph_edge_len);
//}
//
//
//translate([100, 0, 0])
//icosa_faces_sym_ir(30) {
//    regular_polygon_frame(3, $ph_edge_len);
//}

difference() {
    N = 3;
    icosa_faces_sym_ir(30) {
        regular_polygon_frame(N, $ph_edge_len);
    }
    icosa_faces_sym_ir(30) {
        // main facet base
        translate([0, 0, 1]) linear_extrude(height = 5) {
            regular_polygon_2d(N, $ph_edge_len);
        }
    }
}

