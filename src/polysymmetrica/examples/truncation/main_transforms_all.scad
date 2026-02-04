use <../../core/truncation.scad>
use <../../core/cantellation.scad>

use <../../models/platonics_all.scad>
use <../../models/archimedians_all.scad>

use <util_demo.scad>

/**
Showcase default transforms (truncate, rectify, chamfer, cantellate, cantitruncate)
applied to the 5 platonics plus cuboctahedron and icosidodecahedron.
*/

shapes = concat(
    platonics_all(),
    [
        ["cuboctahedron", cuboctahedron()],
        ["icosidodecahedron", icosidodecahedron()]
    ]
);

transforms = [
    ["identity", function(p) p],
    ["truncate", function(p) poly_truncate(p)],
    ["rectify", function(p) poly_rectify(p)],
    ["chamfer", function(p) poly_chamfer(p)],
    ["cantellate", function(p) poly_cantellate(p)],
    ["cantitruncate", function(p) poly_cantitruncate(p)]
];

spacing_x = 140;
spacing_y = 140;
label_h = 1;
label_size = 8;
label_z = -20;

// Column labels
for (s = with_index(shapes)) {
    translate([spacing_x * s[0], -70, label_z])
        rotate([0, 0, -90])
            linear_extrude(height=label_h)
                text(s[1][0], size=label_size, halign="left", valign="center");
}

// Row labels + demos
for (t = with_index(transforms)) {
    translate([-80, spacing_y * t[0], label_z])
        rotate([0, 0, -90])
            linear_extrude(height=label_h)
                text(t[1][0], size=label_size, halign="left", valign="center");

    for (s = with_index(shapes)) {
        p = t[1][1](s[1][1]);
        name = str(s[1][0], " (", t[1][0], ")");
        translate([spacing_x * s[0], spacing_y * t[0], 0]) demo(p, name=name);
    }
}
