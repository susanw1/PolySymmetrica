use <../polysymmetrica/core/funcs.scad>
use <../polysymmetrica/core/placement.scad>
use <../polysymmetrica/core/duals.scad>

use <../polysymmetrica/models/tetrahedron.scad>
use <../polysymmetrica/models/octahedron.scad>
use <../polysymmetrica/models/icosahedron.scad>

use <testing_util.scad>


p1 = octahedron();
p2 = poly_dual(poly_dual(octahedron()));

assert_verts_matches(poly_verts(p1), poly_verts(p2));
assert_facet_matches(p1, p2);

translate([-200,0,0]) placement_tester(poly_dual(tetrahedron()), 60);

translate([0,0,0]) placement_tester(poly_dual(octahedron()), 60);

translate([200,0,0]) placement_tester(poly_dual(icosahedron()), 60);


