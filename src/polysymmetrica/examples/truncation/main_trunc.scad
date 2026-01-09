use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <util_demo.scad>

t = undef;

LAYER1 = -100;

translate([-100, -100, LAYER1]) demo(poly_truncate(tetrahedron(), t));
translate([0, -100, LAYER1]) demo(poly_truncate(octahedron(), t));
translate([100, -100, LAYER1]) demo(poly_truncate(icosahedron(), t));

translate([-100, 0, LAYER1]) demo(poly_rectify(tetrahedron()));
translate([0, 0, LAYER1]) demo(poly_rectify(octahedron()));
translate([100, 0, LAYER1]) demo(poly_rectify(icosahedron()));

translate([0, 100, LAYER1]) demo(poly_truncate(hexahedron(), t));
translate([100, 100, LAYER1]) demo(poly_truncate(dodecahedron(), t));

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
translate([0, 0, LAYER3]) combo(poly_rectify(octahedron()));
translate([100, 0, LAYER3]) combo(poly_rectify(icosahedron()));

translate([0, 100, LAYER3]) combo(poly_truncate(hexahedron()));
translate([100, 100, LAYER3]) combo(poly_truncate(dodecahedron()));

LAYER4 = -200;

translate([-100, -100, LAYER4]) demo(poly_truncate(poly_rectify(tetrahedron())));
translate([0, -100, LAYER4]) demo(poly_truncate(poly_rectify(octahedron()), t=undef));
translate([100, -100, LAYER4]) demo(poly_truncate(poly_rectify(icosahedron()), t = undef));

translate([-100, 0, LAYER4]) demo(poly_rectify(poly_rectify(tetrahedron())));
translate([0, 0, LAYER4]) demo(poly_rectify(poly_rectify(octahedron())));
translate([100, 0, LAYER4]) demo(poly_rectify(poly_rectify(icosahedron())));
