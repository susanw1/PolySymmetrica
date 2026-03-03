use <../../core/attach.scad>
use <../../core/funcs.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

spacing = 120;

items = [
    ["tetra+t (face 0)", poly_attach(tetrahedron(), tetrahedron(), f1=0, f2=0)],
    ["cube+cube (face 0)", poly_attach(hexahedron(), hexahedron(), f1=0, f2=0)],
    ["cube+cube (faces [0,1])", poly_attach(hexahedron(), hexahedron(), f1=[0,1], f2=0)],
    ["cube+cube (rotate_step=1)", poly_attach(hexahedron(), hexahedron(), f1=0, f2=0, rotate_step=1)],
    ["cube+cube (scaled p2, fit_edge)", poly_attach(hexahedron(), [[for (v = poly_verts(hexahedron())) v * 2], poly_faces(hexahedron()), poly_e_over_ir(hexahedron())], f1=0, f2=0, scale_mode="fit_edge")],
    ["octa+icosa (face 0)", poly_attach(octahedron(), icosahedron(), f1=0, f2=0)],
];

for (row = with_index(items, 0)) {
    ix = row[0];
    name = row[1][0];
    p = row[1][1];
    translate([ix * spacing, 0, 0]) {
        demo(p, name=name);
    }
}
