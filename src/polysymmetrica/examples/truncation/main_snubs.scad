use <../../core/truncation.scad>
use <../../core/duals.scad>
use <../../core/render.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

/**
Quick snub demo (snub cube + snub dodecahedron).
*/

spacing = 140;


//translate([-spacing, 0, 0]) demo(poly_snub(tetrahedron(), c = undef, angle = 40), name="snub_tet");
translate([0, 0, 0]) demo(poly_snub(hexahedron()), name="snub_cube");
//translate([spacing, 0, 0]) demo(poly_snub(dodecahedron()), name="snub_dodecahedron");

//b0 = octahedron();
//b1 = poly_rectify(b0);
//translate([spacing*2, 0, 0]) demo(poly_cantellate(base), name="cant_cuboctahedron");
//translate([spacing*3, 0, 0]) demo(poly_snub(b1, angle = undef, c = undef, family_k = 4), name="snub_cuboctahedron");
////
//b2 = poly_dual(b1);
//translate([spacing*4, 0, 0]) demo(poly_snub(b1), name="snub_rh_dod");
////
//b3 = poly_dual(poly_rectify(icosahedron()));
//translate([spacing*5, 0, 0]) demo(poly_snub(b3, angle = 20), name="snub_rh_triac");

//b4 = poly_cantellate(icosahedron());
//poly_render(b4,40);
//translate([0, 0, 0]) demo(poly_snub(b4, angle = 10), ir = 100, name="canttrunc_cuboctahedron");
