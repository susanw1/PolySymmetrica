use <../../core/render.scad>
use <../../models/archimedians_all.scad>
use <../../models/catalans_all.scad>
use <../truncation/util_demo.scad>

// Archimedean + combo + Catalan display

archs = archimedians_all();
cats = catalans_all();

rows = [
    for (a = archs)
        let(
            name = a[0],
            p = a[1],
            d_map = archimedean_to_catalan_name(name),
            d_idx = [for (i = [0:1:len(cats)-1]) if (cats[i][0] == d_map) i][0],
            d_name = is_undef(d_idx) ? str(name, "_dual") : cats[d_idx][0],
            d = is_undef(d_idx) ? undef : cats[d_idx][1]
        )
        [name, p, d_name, d]
];

spacing_x = 140;
spacing_y = 140;
label_height = 1;
label_size = 8;
label_x = -40;
label_z = -40;

for (r = with_index(rows)) {
    i = r[0];
    poly_name = r[1][0];
    p = r[1][1];
    d_name = r[1][2];
    d = r[1][3];

    // Archimedean
    if (!is_undef(p)) {
        translate([spacing_x * i, 0, 0]) demo(p, name = poly_name);
        translate([spacing_x * i + label_x, 0, label_z])
            rotate([0, 0, -90])
                linear_extrude(height=label_height)
                    text(poly_name, size=label_size, halign="center", valign="center");
    } else {
        echo(str("archimedians: ", poly_name, " is TODO"));
    }

    // Combo (archimedean + dual)
    if (!is_undef(p) && !is_undef(d))
        translate([spacing_x * i, spacing_y, 0]) combo(p);

    // Catalan
    if (!is_undef(d)) {
        translate([spacing_x * i, 2 * spacing_y, 0]) demo(d, name = d_name);
        translate([spacing_x * i + label_x, 2 * spacing_y, label_z])
            rotate([0, 0, -90])
                linear_extrude(height=label_height)
                    text(d_name, size=label_size, halign="center", valign="center");
    } else {
        echo(str("catalans: ", d_name, " is TODO"));
    }
}
