use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <util_demo.scad>

LAYER1 = -100;

poly_describe(poly_cantellate(hexahedron(), 0.2, 0, 0));

translate([-100, -100, LAYER1]) demo(poly_cantellate(tetrahedron(), 0.1, 0.1, 0.1));
translate([0, -100, LAYER1]) demo(poly_cantellate(octahedron(), 0, 0, 0));
translate([100, -100, LAYER1]) demo(poly_cantellate(icosahedron(), 0, 0, 0));
//
//translate([-100, 0, LAYER1]) demo(poly_rectify(tetrahedron()));
//translate([0, 0, LAYER1]) demo(poly_rectify(octahedron()));
//translate([100, 0, LAYER1]) demo(poly_rectify(icosahedron()));
//
//translate([0, 100, LAYER1]) demo(poly_cantellate(hexahedron(), 0.2, 0, 0));
//translate([100, 100, LAYER1]) demo(poly_cantellate(dodecahedron(), 0, 0, 0));
//
//LAYER2 = 100;
//
//translate([-100, -100, LAYER2]) demo(poly_dual(poly_cantellate(tetrahedron()), 0, 0, 0));
//translate([0, -100, LAYER2]) demo(poly_dual(poly_cantellate(octahedron()), 0, 0, 0));
//translate([100, -100, LAYER2]) demo(poly_dual(poly_cantellate(icosahedron()), 0, 0, 0));
//
//translate([-100, 0, LAYER2]) demo(poly_dual(poly_rectify(tetrahedron())));
//translate([0, 0, LAYER2]) demo(poly_dual(poly_rectify(octahedron())));
//translate([100, 0, LAYER2]) demo(poly_dual(poly_rectify(icosahedron())));
//
//translate([0, 100, LAYER2]) demo(poly_dual(poly_cantellate(hexahedron()), 0, 0, 0));
//translate([100, 100, LAYER2]) demo(poly_dual(poly_cantellate(dodecahedron()), 0, 0, 0));
//
//LAYER3 = 0;
//
//translate([-100, -100, LAYER3]) combo(poly_cantellate(tetrahedron(), 0, 0, 0));
//translate([0, -100, LAYER3]) combo(poly_cantellate(octahedron(), 0, 0, 0));
//translate([100, -100, LAYER3]) combo(poly_cantellate(icosahedron(), 0, 0, 0));
//
//translate([-100, 0, LAYER3]) combo(poly_rectify(tetrahedron()));
//translate([0, 0, LAYER3]) combo(poly_rectify(octahedron()));
//translate([100, 0, LAYER3]) combo(poly_rectify(icosahedron()));
//
//translate([0, 100, LAYER3]) combo(poly_cantellate(hexahedron(), 0, 0, 0));
//translate([100, 100, LAYER3]) combo(poly_cantellate(dodecahedron(), 0, 0, 0));
//
