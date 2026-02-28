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

// Parameter playground (explicit height/angle variants).
play_y = 300;

prism_cases = [
    ["prism n=6 (tall)", function() poly_prism(6, height_scale=1.5)],
    ["prism n=8 (flat)", function() poly_prism(8, height_scale=0.7)],
    ["prism n=7 (explicit h)", function() poly_prism(7, height=0.8, height_scale=1.25)]
];

antiprism_cases = [
    ["antiprism n=6 (+10deg)", function() poly_antiprism(6, angle=10)],
    ["antiprism n=6 (-10deg)", function() poly_antiprism(6, angle=-10)],
    ["antiprism n=7 (explicit h)", function() poly_antiprism(7, height=0.9, height_scale=1.3)]
];

for (item = with_index(prism_cases)) {
    translate([item[0] * spacing_x, play_y, 0]) demo(item[1][1](), name=item[1][0]);
}

for (item = with_index(antiprism_cases)) {
    translate([item[0] * spacing_x, play_y + spacing_y, 0]) demo(item[1][1](), name=item[1][0]);
}
