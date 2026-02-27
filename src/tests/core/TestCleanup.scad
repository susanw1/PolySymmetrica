use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/models/platonics_all.scad>

EPS = 1e-12;

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

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

function _maxabs3(verts) =
    max([for (p = verts) max([abs(p[0]), abs(p[1]), abs(p[2])])]);

module test_cleanup_regression__truncate_preserves_e_over_ir() {
    p = hexahedron();
    top_vs = _top_face_verts(p);
    some_vs = [for (i = [0:1:min(3, len(top_vs)-1)]) top_vs[i]];

    q_no = poly_truncate(
        p,
        t=0.001,
        params_overrides=[
            ["vert", "id", some_vs, ["t", 0.25]]
        ],
        cleanup=false
    );

    q_yes = poly_truncate(
        p,
        t=0.001,
        params_overrides=[
            ["vert", "id", some_vs, ["t", 0.25]]
        ],
        cleanup=true,
        cleanup_eps=1e-8
    );

    assert_near(
        poly_e_over_ir(q_yes),
        poly_e_over_ir(q_no),
        1e-12,
        "cleanup should preserve operator descriptor scale (e_over_ir)"
    );
    assert_int_eq(len(poly_verts(q_yes)), len(poly_verts(q_no)), "cleanup verts count");
    assert_int_eq(len(poly_faces(q_yes)), len(poly_faces(q_no)), "cleanup faces count");
    assert_near(_maxabs3(poly_verts(q_yes)), _maxabs3(poly_verts(q_no)), 1e-9, "cleanup coordinates");
}

module run_TestCleanup() {
    test_cleanup_regression__truncate_preserves_e_over_ir();
}

