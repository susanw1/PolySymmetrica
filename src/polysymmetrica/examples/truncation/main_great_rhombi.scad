use <../../core/truncation.scad>
use <../../core/render.scad>
use <../../models/archimedians_all.scad>
use <util_demo.scad>

/**
Great rhombi* via cantitruncate uniform solver.
*/

translate([0, 0, 0])
    demo(great_rhombicuboctahedron(4, 2));

translate([80, 0, 0])
    demo(great_rhombicosidodecahedron(4, 2));
