use <../polysymmetrica/core/placement.scad>
use <../polysymmetrica/core/duals.scad>

use <../polysymmetrica/models/tetrahedron.scad>
use <../polysymmetrica/models/octahedron.scad>
use <../polysymmetrica/models/icosahedron.scad>

use <testing_util.scad>


translate([-200,0,0]) placement_tester(poly_dual(tetrahedron()), 60);

translate([0,0,0]) placement_tester(poly_dual(octahedron()), 60, 4);

translate([200,0,0]) placement_tester(poly_dual(icosahedron()), 60, 5);


