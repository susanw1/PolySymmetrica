use <../../core/render.scad>
use <../../models/archimedians_all.scad>
use <../../models/catalans_all.scad>
use <../truncation/util_demo.scad>

// Archimedean + combo + Catalan display

rows = [
    ["truncated_tetrahedron", truncated_tetrahedron(), "triakis_tetrahedron", triakis_tetrahedron()],
    ["truncated_cube", truncated_cube(), "triakis_octahedron", triakis_octahedron()],
    ["truncated_octahedron", truncated_octahedron(), "tetrakis_hexahedron", tetrakis_hexahedron()],
    ["truncated_dodecahedron", truncated_dodecahedron(), "triakis_icosahedron", triakis_icosahedron()],
    ["truncated_icosahedron", truncated_icosahedron(), "pentakis_dodecahedron", pentakis_dodecahedron()],
    ["cuboctahedron", cuboctahedron(), "rhombic_dodecahedron", rhombic_dodecahedron()],
    ["icosidodecahedron", icosidodecahedron(), "rhombic_triacontahedron", rhombic_triacontahedron()],
    ["rhombicuboctahedron", rhombicuboctahedron(), "deltoidal_icositetrahedron", deltoidal_icositetrahedron()],
    ["rhombicosidodecahedron", rhombicosidodecahedron(), "deltoidal_hexecontahedron", deltoidal_hexecontahedron()],
    ["great_rhombicuboctahedron", great_rhombicuboctahedron(), "disdyakis_dodecahedron", disdyakis_dodecahedron()],
    ["great_rhombicosidodecahedron", great_rhombicosidodecahedron(), "disdyakis_triacontahedron", disdyakis_triacontahedron()],
    ["snub_cube", snub_cube(), "pentagonal_icositetrahedron", pentagonal_icositetrahedron()],
    ["snub_dodecahedron", snub_dodecahedron(), "pentagonal_hexecontahedron", pentagonal_hexecontahedron()]
];

spacing_x = 140;
spacing_y = 140;
label_height = 1;
label_size = 8;
label_x = -40;
label_z = -40;

for (i = [0:1:len(rows)-1]) {
    poly_name = rows[i][0];
    p = rows[i][1];
    d_name = rows[i][2];
    d = rows[i][3];

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
