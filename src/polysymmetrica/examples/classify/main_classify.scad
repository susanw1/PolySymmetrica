use <../../core/truncation.scad>
use <../../core/classify.scad>
use <../../core/render.scad>
use <../../models/platonics_all.scad>
use <../../models/archimedians_all.scad>
use <../../models/catalans_all.scad>
use <../truncation/util_demo.scad>

/**
Classification demo: shows face/edge/vertex family counts for a few shapes.

Quick usage:
// Topology only
// poly_classify(p, detail=0);
// One-pass neighbour refinement
// poly_classify(p, detail=1);
// Iterated refinement with radius=2
// poly_classify(p, detail=2, radius=2);
// Include geometry (avg edge length) in keys
// poly_classify(p, detail=1, include_geom=true);
*/

spacing = 120;

module show_classify(p, name, pos=[0,0,0]) {
    translate(pos) demo(p, name=name);
    show_poly(p, detail = 1);
}

show_classify(hexahedron(), "cube", [spacing, 0, 0]);
show_classify(octahedron(), "octa", [spacing*2, 0, 0]);
show_classify(dodecahedron(), "dodeca", [spacing*3, 0, 0]);
show_classify(icosahedron(), "icosa", [spacing*4, 0, 0]);

show_classify(poly_truncate(octahedron()), "trunc_octa", [0, -spacing, 0]);
show_classify(rhombicuboctahedron(), "rhombicubocta", [spacing, -spacing, 0]);
show_classify(great_rhombicuboctahedron(), "great_rhombicubocta", [spacing*2, -spacing, 0]);
show_classify(deltoidal_icositetrahedron(), "deltoidal_icositetra", [spacing*3, -spacing, 0]);
show_classify(disdyakis_dodecahedron(), "disdyakis_dodeca", [spacing*4, -spacing, 0]);



//show_poly(truncated_tetrahedron());
//show_poly(cuboctahedron(), detail=5);
//show_poly(rhombicuboctahedron(), detail = 3);
//show_poly(great_rhombicuboctahedron(), detail = 3);
//show_poly(poly_truncate(rhombic_triacontahedron()), detail = 3);


//show_classify(hexahedron(), "cube", [spacing, 0, 0]);
//show_classify(octahedron(), "octa", [spacing*2, 0, 0]);
//show_classify(dodecahedron(), "dodeca", [spacing*3, 0, 0]);
//show_classify(icosahedron(), "icosa", [spacing*4, 0, 0]);
//
//show_classify(poly_truncate(octahedron()), "trunc_octa", [0, -spacing, 0]);
//show_classify(rhombicuboctahedron(), "rhombicubocta", [spacing, -spacing, 0]);
//show_classify(great_rhombicuboctahedron(), "great_rhombicubocta", [spacing*2, -spacing, 0]);
//show_classify(deltoidal_icositetrahedron(), "deltoidal_icositetra", [spacing*3, -spacing, 0]);
//show_classify(disdyakis_dodecahedron(), "disdyakis_dodeca", [spacing*4, -spacing, 0]);
