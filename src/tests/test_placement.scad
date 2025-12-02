use <../polysymmetrica/core/placement.scad>
use <../polysymmetrica/models/tetrahedron.scad>
use <../polysymmetrica/models/octahedron.scad>
use <../polysymmetrica/models/icosahedron.scad>


echo("Testing face placement:");
place_on_faces_ir(tetrahedron(), 25) sphere(2);

echo("Testing edge placement:");
place_on_edges_ir(octahedron(), 25) sphere(2);

echo("Testing vertex placement:");
place_on_vertices_ir(icosahedron(), 25) sphere(2);
