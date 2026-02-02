use <../../core/render.scad>
use <../../models/archimedians_all.scad>
use <../truncation/util_demo.scad>

// Render all Archimedean solids (snubs are stubs).

archimedians = [
    ["truncated_tetrahedron", truncated_tetrahedron()],
    ["truncated_cube", truncated_cube()],
    ["truncated_octahedron", truncated_octahedron()],
    ["truncated_dodecahedron", truncated_dodecahedron()],
    ["truncated_icosahedron", truncated_icosahedron()],
    ["cuboctahedron", cuboctahedron()],
    ["icosidodecahedron", icosidodecahedron()],
    ["rhombicuboctahedron", rhombicuboctahedron()],
    ["rhombicosidodecahedron", rhombicosidodecahedron()],
    ["great_rhombicuboctahedron", great_rhombicuboctahedron()],
    ["great_rhombicosidodecahedron", great_rhombicosidodecahedron()],
    ["snub_cube", snub_cube()],
    ["snub_dodecahedron", snub_dodecahedron()]
];

spacing = 120;

for (i = [0:1:len(archimedians)-1]) {
    name = archimedians[i][0];
    p = archimedians[i][1];
    if (!is_undef(p))
        translate([spacing * i, 0, 0]) demo(p);
    else
        echo(str("archimedians: ", name, " is TODO"));
}
