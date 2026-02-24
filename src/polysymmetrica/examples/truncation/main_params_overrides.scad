use <../../core/funcs.scad>
use <../../core/truncation.scad>
use <../../models/platonics_all.scad>
use <util_demo.scad>

/**
Showcase selective params_overrides on tetrahedron and cube.

Pattern used throughout:
- set scalar controls to 0 (or near-0) so unaffected elements stay unchanged
- selectively override only chosen face/vertex ids
*/

spacing_x = 150;
spacing_y = 125;
label_size = 8;
label_h = 1;
label_z = -20;

row_labels = [
    "truncate: some corners",
    "cantellate: top face only",
    "chamfer: top face only",
    "cantitruncate: top family slice",
    "snub: top face/verts only"
];

function _idx_of_max(vals) =
    let(
        m = max(vals),
        idxs = [for (i = [0:1:len(vals)-1]) if (vals[i] == m) i]
    )
    idxs[0];

function _top_face_idx(poly) =
    let(
        faces = poly_faces(poly),
        zs = [for (fi = [0:1:len(faces)-1]) poly_face_center(poly, fi, 1)[2]]
    )
    _idx_of_max(zs);

function _top_face_verts(poly) =
    let(
        faces = poly_faces(poly),
        fi = _top_face_idx(poly)
    )
    faces[fi];

function _some_top_verts(poly, n=2) =
    let(vs = _top_face_verts(poly))
    [for (i = [0:1:min(n-1, len(vs)-1)]) vs[i]];

module show_param_ops(poly, poly_name) {
    top_f = _top_face_idx(poly);
    top_vs = _top_face_verts(poly);
    some_vs = _some_top_verts(poly, (len(top_vs) >= 4) ? 4 : 2);

    // 1) Truncate some corners only.
    p0 = poly_truncate(
        poly,
        t=0.001,
        params_overrides=[
            ["vert", "id", some_vs, ["t", 0.25]]
        ]
    );

    // 2) Cantellate just the top face.
    p1 = poly_cantellate(
        poly,
        df=0,
        params_overrides=[
            ["face", "id", [top_f], ["df", 0.14]]
        ]
    );

    // 3) Chamfer just the top face.
    p2 = poly_chamfer(
        poly,
        t=0,
        params_overrides=[
            ["face", "id", [top_f], ["t", 0.30]]
        ]
    );

    // 4) Cantitruncate top slice only (corner t + one face c).
    p3 = poly_cantitruncate(
        poly,
        t=0,
        c=0,
        params_overrides=[
            ["vert", "id", top_vs, ["t", 0.20]],
            ["face", "id", [top_f], ["c", 0.10]]
        ]
    );

    // 5) Snub top face/verts only.
    p4 = poly_snub(
        poly,
        angle=0,
        c=0,
        df=0,
        params_overrides=[
            ["face", "id", [top_f], ["df", 0.08], ["angle", 16]],
            ["vert", "id", top_vs, ["c", 0.06]]
        ]
    );

    rows = [p0, p1, p2, p3, p4];
    for (r = with_index(rows)) {
        y = spacing_y * r[0];
        translate([0, y, 0]) demo(r[1], name=str(poly_name, " - ", row_labels[r[0]]));
    }
}

// Column labels
shape_labels = ["tetrahedron", "hexahedron"];
for (s = with_index(shape_labels)) {
    translate([spacing_x * s[0], -80, label_z])
        rotate([0, 0, -90])
            linear_extrude(height=label_h)
                text(s[1], size=label_size, halign="left", valign="center");
}

// Row labels
for (r = with_index(row_labels)) {
    translate([-95, spacing_y * r[0], label_z])
        rotate([0, 0, -90])
            linear_extrude(height=label_h)
                text(r[1], size=label_size, halign="left", valign="center");
}

translate([0, 0, 0]) show_param_ops(tetrahedron(), "tetrahedron");
translate([spacing_x, 0, 0]) show_param_ops(hexahedron(), "hexahedron");
