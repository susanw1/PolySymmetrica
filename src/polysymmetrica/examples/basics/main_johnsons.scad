use <../../core/construction.scad>
use <../../core/render.scad>
use <../../models/johnsons_all.scad>
use <../truncation/util_demo.scad>

// Johnson previews and direct construction examples.

johnsons = johnsons_all();

spacing_x = 120;
spacing_y = 250;
label_height = 1;
label_size = 8;
label_y = -100;
label_z = -20;

for (j = with_index(johnsons)) {
    i = j[0];
    name = j[1][0];
    p = j[1][1]();

    translate([spacing_x * i, 0, 0]) demo(p, name = name);
    translate([spacing_x * i, label_y, label_z])
        rotate([0, 0, -90])
            linear_extrude(height = label_height)
                text(name, size = label_size, halign = "left", valign = "center");
}

//constructors = [
//    ["poly_pyramid(4)", function() poly_pyramid(4)],
//    ["poly_pyramid(5)", function() poly_pyramid(5)],
//    ["poly_pyramid(5,2)", function() poly_pyramid(5, 2)],
//    ["poly_cupola(3)", function() poly_cupola(3)],
//    ["poly_cupola(4)", function() poly_cupola(4)],
//    ["poly_cupola(5)", function() poly_cupola(5)],
//    ["poly_rotunda()", function() poly_rotunda()],
//    ["poly_elongate(cupola3)", function() poly_elongate(poly_cupola(3), f=0)],
//    ["poly_gyroelongate(cupola3)", function() poly_gyroelongate(poly_cupola(3), f=0)]
//];
//
//for (c = with_index(constructors)) {
//    i = c[0];
//    name = c[1][0];
//    p = c[1][1]();
//
//    translate([spacing_x * i, spacing_y, 0]) demo(p, name = name);
//    translate([spacing_x * i, spacing_y - label_y, label_z])
//        rotate([0, 0, -90])
//            linear_extrude(height = label_height)
//                text(name, size = label_size, halign = "right", valign = "center");
//}
