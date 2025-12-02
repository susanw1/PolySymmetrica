use <../polysymmetrica/core/placement.scad>

use <../polysymmetrica/models/tetrahedron.scad>
use <../polysymmetrica/models/octahedron.scad>
use <../polysymmetrica/models/icosahedron.scad>

use <testing_util.scad>

echo("Tetrahedron:", tetrahedron());
echo("Octahedron:",  octahedron());
echo("Icosahedron:", icosahedron());

translate([-200,0,0]) placement_tester(tetrahedron(), 60);

translate([0,0,0]) placement_tester(octahedron(), 60);

translate([200,0,0]) placement_tester(icosahedron(), 60);
