use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <util_demo.scad>

LAYER1 = -100;

demo(poly_cantellate(hexahedron(), 0.1));
poly_describe(poly_cantellate(hexahedron(), 0.1));
//demo(hexahedron());
poly_describe(hexahedron());

//poly_describe(poly_cantellate(tetrahedron()));
