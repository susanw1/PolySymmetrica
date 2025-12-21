use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>

use <../../models/regular_all.scad>


IR = 30;

COLORS = [ "", "", "", "yellow", "red", "green", "blue", "gray", "red", "white", "red" ];

module demo(p) {
    place_on_faces_ir(p, IR) {
        let (col = COLORS[$ps_vertex_count]) {
            color(col) {
                linear_extrude(height = 1) polygon(points = $ps_face_pts2d);
//                translate([0,0,3]) text(str($ps_facet_idx), halign="center",valign="center", size=4);
            }
        }
    }

    color("silver") 
    place_on_edges_ir(p, IR) {
        cube([$ps_edge_len, 2, 2], center = true);
        translate($ps_poly_center_local) cube([1,1,$ps_edge_midradius], center = false);
//        translate([0,0,2]) text(str($ps_edge_idx), halign="center",valign="center", size=3);
    }

    color("gold") 
    place_on_vertices_ir(p, IR) {
        sphere(2, $fn=30);
//        translate([0,0,2]) text(str($ps_vertex_idx), halign="center",valign="center", size=2);
    }
}

LAYER1 = -100;

//translate([-100, 0, LAYER1]) demo(poly_truncate(tetrahedron()));
//translate([0, 0, LAYER1]) demo(poly_truncate(octahedron()));
//translate([100, 0, LAYER1]) demo(poly_truncate(icosahedron()));
//
//translate([0, 100, LAYER1]) demo(poly_truncate(hexahedron()));
//translate([100, 100, LAYER1]) demo(poly_truncate(dodecahedron()));

LAYER2 = -100;

//translate([-100, 0, LAYER2]) demo(poly_dual(poly_truncate(tetrahedron())));
//translate([0, 0, LAYER2]) demo(poly_dual(poly_truncate(octahedron())));
//translate([100, 0, LAYER2]) demo(poly_dual(poly_truncate(icosahedron())));
////
//translate([0, 100, LAYER2]) demo(poly_dual(poly_truncate(hexahedron())));
//translate([100, 100, LAYER2]) demo(poly_dual(poly_truncate(dodecahedron())));





use <../../core/funcs.scad>



trunc_tet = poly_truncate(tetrahedron());
//color("yellow", 0.4) polyhedron(poly_verts(trunc_tet), poly_faces(trunc_tet));
triakis_tet = poly_dual(trunc_tet);
//color("blue", 0.4) polyhedron(poly_verts(triakis_tet)*1.7, poly_faces(triakis_tet));

color("blue", 0.4) 
place_on_faces_ir(trunc_tet, IR) {
    circle(r = $ps_facet_radius, $fn = $ps_vertex_count);
}
color("yellow", 1) 
place_on_faces_ir(triakis_tet, IR) {
    polygon(points = $ps_face_pts2d);
}


