use <../../core/prisms.scad>
use <../truncation/util_demo.scad>

// Regular prisms and antiprisms.

ns = [3, 4, 5, 6, 7, 8];
spacing_x = 120;
spacing_y = 140;
label_size = 8;
label_h = 1;
label_z = -30;

for (n_pair = with_index(ns)) {
    i = n_pair[0];
    n = n_pair[1];
    x = spacing_x * i;

    p_prism = poly_prism(n);
    p_antiprism = poly_antiprism(n);

    translate([x, 0, 0]) demo(p_prism, name=str("prism n=", n));
    translate([x, spacing_y, 0]) demo(p_antiprism, name=str("antiprism n=", n));

    translate([x, -50, label_z])
        rotate([0, 0, -90])
            linear_extrude(height=label_h)
                text(str("prism n=", n), size=label_size, halign="left", valign="center");

    translate([x, spacing_y + 50, label_z])
        rotate([0, 0, -90])
            linear_extrude(height=label_h)
                text(str("antiprism n=", n), size=label_size, halign="right", valign="center");
}
