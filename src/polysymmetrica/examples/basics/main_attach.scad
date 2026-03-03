use <../../core/attach.scad>
use <../../core/classify.scad>
use <../../core/funcs.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

spacing = 120;

octa = octahedron();
tetra = tetrahedron();
tetra_scaled = [[for (v = poly_verts(tetra)) v * 1.6], poly_faces(tetra), poly_e_over_ir(tetra)];

trunc_dod = poly_truncate(dodecahedron());
// Face-size selection only; topology detail is sufficient here.
cls_td = poly_classify(trunc_dod, detail=0);
tri_faces_td = ps_classify_face_idxs_by_n(cls_td, 3);

items = [
    ["octa + tetra", poly_attach(octa, tetra)],
    ["octa + tetra (rotate_step=1)", poly_attach(octa, tetra, f1=0, f2=0, rotate_step=1)],
    ["octa + scaled tetra (fit_edge)", poly_attach(octa, tetra_scaled, f1=0, f2=0, scale_mode="fit_edge")],
    ["octa + tetra (fit_edge reference)", poly_attach(octa, tetra, f1=0, f2=0, scale_mode="fit_edge")],
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
