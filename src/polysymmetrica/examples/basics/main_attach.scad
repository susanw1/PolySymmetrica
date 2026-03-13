use <../../core/construction.scad>
use <../../core/classify.scad>
use <../../core/funcs.scad>
use <../../core/truncation.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

spacing = 120;

octa = octahedron();
tetra = tetrahedron();
// Asymmetric triangular attachee so rotate_step/mirror effects are visible.
tetra_skew = let(
    v = poly_verts(tetra),
    f = poly_faces(tetra),
    apex = v[3] + [0.28, -0.12, 0.22]
) [[v[0], v[1], v[2], apex], f, poly_e_over_ir(tetra)];
tetra_skew_scaled = [[for (v = poly_verts(tetra_skew)) v * 1.6], poly_faces(tetra_skew), poly_e_over_ir(tetra_skew)];

trunc_dod = poly_truncate(dodecahedron());
// Face-size selection only; topology detail is sufficient here.
cls_td = poly_classify(trunc_dod, detail=0);
tri_faces_td = ps_classify_face_idxs_by_n(cls_td, 3);

items = [
    ["octa + skew tetra", poly_attach(octa, tetra_skew)],
    ["octa + skew tetra (rotate_step=1)", poly_attach(octa, tetra_skew, rotate_step=1)],
    ["octa + skew tetra (mirror=true)", poly_attach(octa, tetra_skew, mirror=true)],
    ["octa + scaled skew tetra (fit_edge)", poly_attach(octa, tetra_skew_scaled, scale_mode="fit_edge")],
    ["octa + tetra (faces [0,2,4])", poly_attach(octa, tetra, f1=[0,2,4], f2=0)],
    ["trunc-dod + tetra on triangle family", poly_attach(trunc_dod, tetra, f1=tri_faces_td, f2=0)],
];

for (row = with_index(items, 0)) {
    ix = row[0];
    name = row[1][0];
    p = row[1][1];
    translate([ix * spacing, 0, 0]) {
        demo(p, name=name);
    }
}
