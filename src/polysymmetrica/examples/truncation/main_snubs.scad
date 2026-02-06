use <../../core/truncation.scad>
use <../../core/duals.scad>
use <../../core/render.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

/**
Quick snub demo (snub cube + snub dodecahedron).
*/

spacing = 140;

//translate([-spacing, 0, 0]) demo(poly_snub(tetrahedron(), angle = 5, c = 0.05, t = 5), name="snub_tet");
//translate([-spacing, 0, 0]) demo(poly_snub(tetrahedron(), angle = 5, c = 0.05), name="snub_tet");
//translate([0, 0, 0]) demo(poly_snub(hexahedron(), angle = 20), name="snub_cube");
//translate([spacing, 0, 0]) demo(poly_snub(dodecahedron(), 20), name="snub_dodecahedron");


b1 = poly_rectify(octahedron());
////translate([spacing*2, 0, 0]) demo(poly_cantellate(base), name="cant_cuboctahedron");
//translate([spacing*3, 0, 0]) demo(poly_snub(base, angle = 5, family_k = 4), name="snub_cuboctahedron");

b2 = poly_dual(b1);
translate([spacing*4, 0, 0]) demo(poly_snub(b2, angle = 10, c = 1.3, family_k = 4), name="snub_rh_dod");

b3 = poly_dual(poly_rectify(icosahedron()));
translate([spacing*5, 0, 0]) demo(poly_snub(b3, angle = 20, c = 1.3), name="snub_rh_triac");

