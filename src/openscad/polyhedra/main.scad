use <edge-mount.scad>
use <tetrahedron.scad>
use <octahedron.scad>
use <icosahedron.scad>


translate([-100, 0, 0])
tetra_faces_sym_ir(30) {
    regular_polygon_frame(3, $ph_edge_len);
}

translate([0, 0, 0])
octa_faces_sym_ir(30) {
    regular_polygon_frame(3, $ph_edge_len);
}


translate([100, 0, 0])
icosa_faces_sym_ir(30) {
    regular_polygon_frame(3, $ph_edge_len);
}