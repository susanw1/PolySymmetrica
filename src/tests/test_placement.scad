use <../polysymmetrica/core/placement.scad>
use <../polysymmetrica/models/tetrahedron.scad>
use <../polysymmetrica/models/octahedron.scad>
use <../polysymmetrica/models/icosahedron.scad>


echo("Testing face placement:");
translate([-100,0,0]) place_on_faces(tetrahedron(), 25) sphere(2);

echo("Testing edge placement:");
translate([0,0,0]) place_on_edges(octahedron(), 25) sphere(2);

echo("Testing vertex placement:");
translate([100,0,0]) place_on_vertices(icosahedron(), 25) sphere(2);
