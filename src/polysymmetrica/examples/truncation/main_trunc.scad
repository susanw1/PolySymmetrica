use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <util_demo.scad>

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
translate([0, 0, LAYER3]) combo(poly_rectify(octahedron()));
translate([100, 0, LAYER3]) combo(poly_rectify(icosahedron()));

translate([0, 100, LAYER3]) combo(poly_truncate(hexahedron()));
translate([100, 100, LAYER3]) combo(poly_truncate(dodecahedron()));

