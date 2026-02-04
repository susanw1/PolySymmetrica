use <../../core/render.scad>
use <../../models/archimedians_all.scad>
use <../truncation/util_demo.scad>

// Render all Archimedean solids (snubs are stubs).

archimedians = archimedians_all();

spacing = 120;
label_height = 1;
label_size = 8;
label_y = -60;
label_z = -15;

for (a = with_index(archimedians)) {
    i = a[0];
    name = a[1][0];
    p = a[1][1];
    if (!is_undef(p)) {
        translate([spacing * i, 0, 0]) demo(p, name=name);
        translate([spacing * i, label_y, label_z])
            rotate([0, 0, -90])
                linear_extrude(height=label_height)
                    text(name, size=label_size, halign="left", valign="center");
    } else {
        echo(str("archimedians: ", name, " is TODO"));
        translate([spacing * i, label_y, label_z])
            rotate([0, 0, -90])
                linear_extrude(height=label_height)
                    text(str(name, " (TODO)"), size=label_size, halign="left", valign="center");
    }
}
