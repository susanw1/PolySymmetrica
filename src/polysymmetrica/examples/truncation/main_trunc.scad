use <../../core/truncation.scad>
use <../../core/duals.scad>

use <../../models/platonics_all.scad>

use <util_demo.scad>

LAYER1 = -100;

rows1 = [
    [ poly_truncate(tetrahedron()), poly_truncate(octahedron()), poly_truncate(icosahedron()) ],
    [ poly_rectify(tetrahedron()), poly_rectify(octahedron()), poly_rectify(icosahedron()) ],
    [ undef, poly_truncate(hexahedron()), poly_truncate(dodecahedron()) ]
];

for (r = with_index(rows1, -1))
    for (c = with_index(r[1], -1))
        if (!is_undef(c[1]))
            translate([100 * c[0], 100 * r[0], LAYER1]) demo(c[1]);

LAYER2 = 100;

rows2 = [
    [ poly_dual(poly_truncate(tetrahedron())), poly_dual(poly_truncate(octahedron())), poly_dual(poly_truncate(icosahedron())) ],
    [ poly_dual(poly_rectify(tetrahedron())), poly_dual(poly_rectify(octahedron())), poly_dual(poly_rectify(icosahedron())) ],
    [ undef, poly_dual(poly_truncate(hexahedron())), poly_dual(poly_truncate(dodecahedron())) ]
];

for (r = with_index(rows2, -1))
    for (c = with_index(r[1], -1))
        if (!is_undef(c[1]))
            translate([100 * c[0], 100 * r[0], LAYER2]) demo(c[1]);

LAYER3 = 0;

rows3 = [
    [ poly_truncate(tetrahedron()), poly_truncate(octahedron()), poly_truncate(icosahedron()) ],
    [ poly_rectify(tetrahedron()), poly_rectify(octahedron()), poly_rectify(icosahedron()) ],
    [ undef, poly_truncate(hexahedron()), poly_truncate(dodecahedron()) ]
];

for (r = with_index(rows3, -1))
    for (c = with_index(r[1], -1))
        if (!is_undef(c[1]))
            translate([100 * c[0], 100 * r[0], LAYER3]) combo(c[1]);

LAYER4 = -200;

rows4 = [
    [ poly_truncate(poly_rectify(tetrahedron())), poly_truncate(poly_rectify(octahedron())), poly_truncate(poly_rectify(icosahedron())) ],
    [ poly_rectify(poly_rectify(tetrahedron())), poly_rectify(poly_rectify(octahedron())), poly_rectify(poly_rectify(icosahedron())) ]
];

for (r = with_index(rows4, -1))
    for (c = with_index(r[1], -1))
        if (!is_undef(c[1]))
            translate([100 * c[0], 100 * r[0], LAYER4]) demo(c[1]);
