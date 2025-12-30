use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>


IR = 30;
T = 0.01;

COLORS = [ "", "", "", "yellow", "red", "green", "blue", "gray", "red", "white", "red" ];

module demo(p, ir = IR) {
    place_on_faces(p, ir) {
        let (col = COLORS[$ps_vertex_count]) {
            color(col) {
                linear_extrude(height=T) polygon(points = $ps_face_pts2d);
//                translate([0,0,3]) text(str($ps_facet_idx), halign="center",valign="center", size=4);
            }
        }
    }

    color("silver")
    place_on_edges(p, ir) {
        cube([$ps_edge_len, 1, 1], center = true);
//        translate($ps_poly_center_local) cube([1,1,$ps_edge_midradius], center = false);
//        translate([0,0,2]) text(str($ps_edge_idx), halign="center",valign="center", size=3);
    }

    color("gold")
    place_on_vertices(p, ir) {
        sphere(1.5, $fn=30);
//        translate([0,0,2]) text(str($ps_vertex_idx), halign="center",valign="center", size=2);
    }
}

module combo(p, scale_f = function(p,d) scale_dual_edge_cross(p,d, 0)) {
    d = poly_dual(p);
//    k = max([for (f = poly_faces(p)) len(f)]);
    m = scale_f(p, d);        
    // m = scale_dual_face_radius(p, d, face_k = k, dual_face_k = undef);
    echo ("m", m);
    color("blue", 1)
    place_on_faces(p, IR) {
        translate([0,0,-T]) linear_extrude(height=T) polygon(points = $ps_face_pts2d);
    }
    color("yellow", 1)
    place_on_faces(d, IR * m) {
        translate([0,0,-T]) linear_extrude(height=T) polygon(points = $ps_face_pts2d);
    }
}


LAYER1 = -100;

translate([-100, -100, LAYER1]) demo(poly_truncate(tetrahedron()));
translate([0, -100, LAYER1]) demo(poly_truncate(octahedron()));
translate([100, -100, LAYER1]) demo(poly_truncate(icosahedron()));

translate([-100, 0, LAYER1]) demo(poly_rectify(tetrahedron()));
translate([0, 0, LAYER1]) demo(poly_rectify(octahedron()));
translate([100, 0, LAYER1]) demo(poly_rectify(icosahedron()));

translate([0, 100, LAYER1]) demo(poly_truncate(hexahedron()));
translate([100, 100, LAYER1]) demo(poly_truncate(dodecahedron()));

LAYER2 = 100;

translate([-100, -100, LAYER2]) demo(poly_dual(poly_truncate(tetrahedron())));
translate([0, -100, LAYER2]) demo(poly_dual(poly_truncate(octahedron())));
translate([100, -100, LAYER2]) demo(poly_dual(poly_truncate(icosahedron())));

translate([-100, 0, LAYER2]) demo(poly_dual(poly_rectify(tetrahedron())));
translate([0, 0, LAYER2]) demo(poly_dual(poly_rectify(octahedron())));
translate([100, 0, LAYER2]) demo(poly_dual(poly_rectify(icosahedron())));

translate([0, 100, LAYER2]) demo(poly_dual(poly_truncate(hexahedron())));
translate([100, 100, LAYER2]) demo(poly_dual(poly_truncate(dodecahedron())));

LAYER3 = 0;

translate([-100, -100, LAYER3]) combo(poly_truncate(tetrahedron()));
translate([0, -100, LAYER3]) combo(poly_truncate(octahedron()));
translate([100, -100, LAYER3]) combo(poly_truncate(icosahedron()));

translate([-100, 0, LAYER3]) combo(poly_rectify(tetrahedron()));
translate([0, 0, LAYER3]) combo(poly_rectify(octahedron()), function(p,d) scale_dual_edge_cross(p,d,0,0));
translate([100, 0, LAYER3]) combo(poly_rectify(icosahedron()), function(p,d) scale_dual_edge_cross(p,d,0,0)); // 1.0119

translate([0, 100, LAYER3]) combo(poly_truncate(hexahedron()));
translate([100, 100, LAYER3]) combo(poly_truncate(dodecahedron()));

