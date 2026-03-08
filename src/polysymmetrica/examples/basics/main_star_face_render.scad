use <../../core/prisms.scad>
use <../../core/render.scad>
use <../../core/placement.scad>
use <../../core/segments.scad>
use <../../core/funcs.scad>
use <../truncation/util_demo.scad>

// Focused repro/inspection scene for star-faced rendering behavior.
// Compare the same poly rendered in four ways:
//   1) poly_render(...) -> OpenSCAD polyhedron(points, faces)
//   2) demo(...)         -> util_demo face fill (triangle fan)
//   3) face-local polygon() plates (even-odd behavior from 2D polygon filling)
//   4) face segmentation plates via place_on_face_segments(...)

IR = 30;
TXT_Z = -34;
TXT_H = 1;
TXT_S = 6;

p = poly_antiprism(5, p=2, angle=0);

module _label(s) {
    translate([0, -55, TXT_Z])
        rotate([0, 0, -90])
            linear_extrude(height=TXT_H)
                text(s, size=TXT_S, halign="center", valign="center");
}

function _seg_color(i) =
    (i % 6 == 0) ? "lightblue" :
    (i % 6 == 1) ? "paleturquoise" :
    (i % 6 == 2) ? "powderblue" :
    (i % 6 == 3) ? "skyblue" :
    (i % 6 == 4) ? "aliceblue" :
    "cadetblue";

function _seg_centroid2d(pts) =
    (len(pts) == 0) ? [0,0] : v_scale(v_sum(pts), 1 / len(pts));

function _seg_mid2d(seg) =
    [ (seg[0][0] + seg[1][0]) / 2, (seg[0][1] + seg[1][1]) / 2 ];

module _face_polygon_fill(poly, ir=IR, thk=0.02) {
    place_on_faces(poly, ir) {
        place_on_face_segments(mode="nonzero") {
            color("lightblue")
                linear_extrude(height=thk)
                    polygon(points = $ps_seg_pts2d);
        }
    }
    color("silver")
    place_on_edges(poly, ir) cube([$ps_edge_len, 0.8, 0.8], center=true);
    color("gold")
    place_on_vertices(poly, ir) sphere(1.2, $fn=12);
}

module _face_segment_fill(poly, ir=IR, thk=0.02) {
    place_on_faces(poly, ir) {
        place_on_face_segments(mode="evenodd") {
            color(_seg_color($ps_seg_idx), 1)
                linear_extrude(height=thk)
                    polygon(points = $ps_seg_pts2d);
            color("white")
                translate([_seg_centroid2d($ps_seg_pts2d)[0], _seg_centroid2d($ps_seg_pts2d)[1], thk + 0.01])
                    linear_extrude(height=0.08)
                        text(str($ps_face_idx, ":", $ps_seg_idx), size=2.2, halign="center", valign="center");
            color("white")
                translate([_seg_centroid2d($ps_seg_pts2d)[0], _seg_centroid2d($ps_seg_pts2d)[1], -0.09])
                    mirror([0, 0, 1])
                        linear_extrude(height=0.08)
                            text(str($ps_face_idx, ":", $ps_seg_idx), size=2.2, halign="center", valign="center");
        }
    }
    color("silver")
    place_on_edges(poly, ir) cube([$ps_edge_len, 0.8, 0.8], center=true);
    color("gold")
    place_on_vertices(poly, ir) sphere(1.2, $fn=12);
}

module _face_geom_cut_stencil(poly, ir=IR, thk=0.4) {
    place_on_faces(poly, ir) {
        // Focus on side panels where penetrating geometry matters most.
        if ($ps_vertex_count == 3) { // && $ps_face_idx == 2) {
            render() color("lightblue")
                difference() {
                    linear_extrude(height=thk, center=true)
                        polygon(points = $ps_face_pts2d);
                    face_cut_stencil(face_thk=thk, kerf=1.6, extend=0.6, z_pad=0.2, eps=1e-8, filter_parent=true);
                }
            place_on_face_geom_cut_segments(eps=1e-8, filter_parent=true) {
                color("black")
                    translate([_seg_mid2d($ps_face_cut_segment2d_local)[0], _seg_mid2d($ps_face_cut_segment2d_local)[1], thk/2 + 0.01])
                        linear_extrude(height=0.08)
                            text(str($ps_face_cut_idx), size=1.4, halign="center", valign="center");
            }
        } else {
            place_on_face_segments(mode="nonzero") {
                color("gainsboro")
                    linear_extrude(height=thk, center=true)
                        polygon(points = $ps_seg_pts2d);
            }
        }
    }
    color("silver")
    place_on_edges(poly, ir) cube([$ps_edge_len, 0.8, 0.8], center=true);
    color("gold")
    place_on_vertices(poly, ir) sphere(1.2, $fn=12);
}

translate([-100, 0, 0]) {
    poly_render(p, IR);
    _label("poly_render (polyhedron)");
}

translate([0, 0, 0]) {
    demo(p, ir=IR, name="star antiprism demo");
    _label("demo (fan face-fill)");
}

translate([100, 0, 0]) {
    _face_polygon_fill(p, IR);
    _label("face-local polygon() fill");
}

translate([200, 0, 0]) {
    _face_segment_fill(p, IR);
    _label("face-segment fill");
}

translate([300, 0, 0]) {
    _face_geom_cut_stencil(p, IR);
    _label("geom-cut stencil");
}

// Optional tiny raw repro: same crossing loop in 2D polygon vs 3D polyhedron.
// Some OpenSCAD backends fill these differently for self-intersecting loops.
SHOW_MIN_REPRO = false;
if (SHOW_MIN_REPRO) {
    pts2d = [[0,9],[-5,-5],[8,3],[-8,3],[5,-5]];
    pts3d_top = [for (q = pts2d) [q[0], q[1], 1]];
    pts3d_bot = [for (q = pts2d) [q[0], q[1], -1]];
    pts3d = concat(pts3d_top, pts3d_bot);
    faces3d = concat(
        [[0,1,2,3,4]],
        [[9,8,7,6,5]],
        [for (i = [0:1:4]) [i, (i+1)%5, 5+((i+1)%5), 5+i]]
    );

    translate([-35, 125, 0]) linear_extrude(height=0.5) polygon(points=pts2d);
    color("silver") translate([35, 125, 0]) polyhedron(points=pts3d, faces=faces3d, convexity=10);
}
