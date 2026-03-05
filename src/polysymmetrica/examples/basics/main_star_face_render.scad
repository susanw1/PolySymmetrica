use <../../core/prisms.scad>
use <../../core/render.scad>
use <../../core/placement.scad>
use <../truncation/util_demo.scad>

// Focused repro/inspection scene for star-faced rendering behavior.
// Compare the same poly rendered in three ways:
//   1) poly_render(...) -> OpenSCAD polyhedron(points, faces)
//   2) demo(...)         -> util_demo face fill (triangle fan)
//   3) face-local polygon() plates (even-odd behavior from 2D polygon filling)

IR = 30;
TXT_Z = -34;
TXT_H = 1;
TXT_S = 6;

p = poly_antiprism(5, p=2);

module _label(s) {
    translate([0, -55, TXT_Z])
        rotate([0, 0, -90])
            linear_extrude(height=TXT_H)
                text(s, size=TXT_S, halign="center", valign="center");
}

module _face_polygon_fill(poly, ir=IR, thk=0.02) {
    place_on_faces(poly, ir) {
        color("lightblue")
            linear_extrude(height=thk)
                polygon(points = $ps_face_pts2d);
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

// Optional tiny raw repro: same crossing loop in 2D polygon vs 3D polyhedron.
// Some OpenSCAD backends fill these differently for self-intersecting loops.
SHOW_MIN_REPRO = true;
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

