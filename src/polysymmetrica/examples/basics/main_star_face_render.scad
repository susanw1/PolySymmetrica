use <../../core/prisms.scad>
use <../../core/placement.scad>
use <../../core/segments.scad>

IR = 30;
TXT_Z = -34;
TXT_H = 1;
TXT_S = 6;

p = poly_antiprism(5, p=2, angle=0);

module _label(s) {
    translate([0, -55, TXT_Z])
        rotate([0, 0, -90])
            linear_extrude(height = TXT_H)
                text(s, size = TXT_S, halign = "center", valign = "center");
}

function _seg_mid2d(seg) =
    [ (seg[0][0] + seg[1][0]) / 2, (seg[0][1] + seg[1][1]) / 2 ];

module _star_face_polygon_fill(poly, ir=IR, thk=0.4) {
    place_on_faces(poly, ir) {
        if ($ps_face_idx == 0)
            color("lightblue")
                linear_extrude(height = thk / 2, center = true)
                    ps_polygon(points = $ps_face_pts2d, mode = "nonzero");
    }
    color("silver")
    place_on_edges(poly, ir) cube([$ps_edge_len, 0.8, 0.8], center = true);
    color("gold")
    place_on_vertices(poly, ir) sphere(1.2, $fn = 12);
}

module _star_face_geom_cut_stencil(poly, ir=IR, thk=0.4) {
    place_on_faces(poly, ir) {
        if ($ps_vertex_count == 3) {
            render()
            color("lightblue")
                difference() {
                    linear_extrude(height = thk, center = true)
                        ps_polygon(points = $ps_face_pts2d, mode = "nonzero");
                    face_cut_stencil(face_thk = thk, kerf = 1.6, extend = 0.6, z_pad = 0.2, eps = 1e-8, filter_parent = true);
                }
            place_on_face_geom_cut_segments(mode = "nonzero", eps = 1e-8, filter_parent = true) {
                color("black")
                    translate([_seg_mid2d($ps_face_cut_segment2d_local)[0], _seg_mid2d($ps_face_cut_segment2d_local)[1], thk / 2 + 0.01])
                        linear_extrude(height = 0.08)
                            text(str($ps_face_cut_idx), size = 1.4, halign = "center", valign = "center");
            }
        } else {
            color("gainsboro")
                linear_extrude(height = thk / 2, center = true)
                    ps_polygon(points = $ps_face_pts2d, mode = "nonzero");
        }
    }
    color("silver")
    place_on_edges(poly, ir) cube([$ps_edge_len, 0.8, 0.8], center = true);
    color("gold")
    place_on_vertices(poly, ir) sphere(1.2, $fn = 12);
}

translate([-70, 0, 0]) {
    _star_face_polygon_fill(p, IR);
    _label("ps_polygon fill");
}

translate([90, 0, 0]) {
    _star_face_geom_cut_stencil(p, IR);
    _label("geom-cut stencil");
}
