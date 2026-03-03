use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/attach.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/platonics_all.scad>

EPS = 1e-8;

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}

function _poly_scaled(poly, s) =
    [[for (v = poly_verts(poly)) v * s], poly_faces(poly), poly_e_over_ir(poly)];

module test_poly_attach__tetra_to_tetra_counts() {
    p = poly_attach(tetrahedron(), tetrahedron(), f1=0, f2=0);
    assert(poly_valid(p, "closed"), "attach tetra+t tetra should be closed valid");
    assert_int_eq(len(poly_verts(p)), 5, "attach tetra+t verts");
    assert_int_eq(len(poly_faces(p)), 6, "attach tetra+t faces");
    assert_int_eq(len(poly_edges(p)), 9, "attach tetra+t edges");
}

module test_poly_attach__cube_to_cube_counts() {
    p = poly_attach(hexahedron(), hexahedron(), f1=0, f2=0);
    assert(poly_valid(p, "closed"), "attach cube+cube should be closed valid");
    assert_int_eq(len(poly_verts(p)), 12, "attach cube+cube verts");
    assert_int_eq(len(poly_faces(p)), 10, "attach cube+cube faces");
    assert_int_eq(len(poly_edges(p)), 20, "attach cube+cube edges");
}

module test_poly_attach__rotate_step_modulo() {
    p1 = poly_attach(hexahedron(), hexahedron(), f1=0, f2=0, rotate_step=1);
    p2 = poly_attach(hexahedron(), hexahedron(), f1=0, f2=0, rotate_step=5);
    assert_int_eq(len(poly_verts(p1)), len(poly_verts(p2)), "attach rotate modulo verts len");
    assert_int_eq(len(poly_faces(p1)), len(poly_faces(p2)), "attach rotate modulo faces len");
    for (i = [0:1:len(poly_verts(p1))-1]) {
        assert_near(norm(poly_verts(p1)[i] - poly_verts(p2)[i]), 0, 1e-10, str("attach rotate modulo vert#", i));
    }
}

module test_poly_attach__fit_edge_matches_unscaled_reference() {
    ref = poly_attach(hexahedron(), hexahedron(), f1=0, f2=0, scale_mode="none");
    fit = poly_attach(hexahedron(), _poly_scaled(hexahedron(), 2), f1=0, f2=0, scale_mode="fit_edge");
    assert_int_eq(len(poly_verts(ref)), len(poly_verts(fit)), "attach fit-edge verts len");
    assert_int_eq(len(poly_faces(ref)), len(poly_faces(fit)), "attach fit-edge faces len");
    for (i = [0:1:len(poly_verts(ref))-1]) {
        assert_near(norm(poly_verts(ref)[i] - poly_verts(fit)[i]), 0, 1e-10, str("attach fit-edge vert#", i));
    }
}

module run_TestAttach() {
    test_poly_attach__tetra_to_tetra_counts();
    test_poly_attach__cube_to_cube_counts();
    test_poly_attach__rotate_step_modulo();
    test_poly_attach__fit_edge_matches_unscaled_reference();
}

