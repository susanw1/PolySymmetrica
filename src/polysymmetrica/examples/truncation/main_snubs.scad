use <../../core/truncation.scad>
use <../../core/render.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

/**
Quick snub demo (snub cube + snub dodecahedron).
*/

spacing = 140;

translate([0, 0, 0]) demo(poly_snub(hexahedron()), name="snub_cube");
//translate([spacing, 0, 0]) demo(poly_snub(dodecahedron()), name="snub_dodecahedron");
