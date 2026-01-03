use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <util_demo.scad>

LAYER1 = -100;

DF = 0.2;
DV = 0.0;
DE = 0.0;

poly_describe(poly_cantellate(hexahedron(), DF, DV, DE));

//translate([-100, -100, LAYER1]) demo(poly_cantellate(tetrahedron(), DF, DV, DE));
//translate([0, -100, LAYER1]) demo(poly_cantellate(octahedron(), DF, DV, DE));
//translate([100, -100, LAYER1]) demo(poly_cantellate(icosahedron(), DF, DV, DE));
//
//translate([-100, 0, LAYER1]) demo(poly_rectify(tetrahedron()));
//translate([0, 0, LAYER1]) demo(poly_rectify(octahedron()));
//translate([100, 0, LAYER1]) demo(poly_rectify(icosahedron()));
//
//translate([0, 100, LAYER1]) demo(poly_cantellate(hexahedron(), DF, DV, DE));
//translate([100, 100, LAYER1]) demo(poly_cantellate(dodecahedron(), DF, DV, DE));
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
//translate([-100, -100, LAYER3]) combo(poly_cantellate(tetrahedron(), DF, DV, DE));
//translate([0, -100, LAYER3]) combo(poly_cantellate(octahedron(), DF, DV, DE));
//translate([100, -100, LAYER3]) combo(poly_cantellate(icosahedron(), DF, DV, DE));
//
//translate([-100, 0, LAYER3]) combo(poly_rectify(tetrahedron()));
//translate([0, 0, LAYER3]) combo(poly_rectify(octahedron()));
//translate([100, 0, LAYER3]) combo(poly_rectify(icosahedron()));
//
//translate([0, 100, LAYER3]) combo(poly_cantellate(hexahedron(), DF, DV, DE));
//translate([100, 100, LAYER3]) combo(poly_cantellate(dodecahedron(), DF, DV, DE));
//
