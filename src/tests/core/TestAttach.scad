use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/construction.scad>
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

function _tetra_skew_apex() =
    let(
        t = tetrahedron(),
        v = poly_verts(t),
        f = poly_faces(t),
        apex = v[3] + [0.25, -0.1, 0.2],
        v2 = [v[0], v[1], v[2], apex]
    )
    [v2, f, poly_e_over_ir(t)];

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

module test_poly_attach__f1_scalar_shorthand_equals_singleton_list() {
    p_scalar = poly_attach(hexahedron(), hexahedron(), f1=0, f2=0);
    p_list = poly_attach(hexahedron(), hexahedron(), f1=[0], f2=0);
    assert_int_eq(len(poly_verts(p_scalar)), len(poly_verts(p_list)), "attach f1 scalar/list verts len");
    assert_int_eq(len(poly_faces(p_scalar)), len(poly_faces(p_list)), "attach f1 scalar/list faces len");
    for (i = [0:1:len(poly_verts(p_scalar))-1]) {
        assert_near(norm(poly_verts(p_scalar)[i] - poly_verts(p_list)[i]), 0, 1e-10, str("attach f1 scalar/list vert#", i));
    }
}

module test_poly_attach__multi_face_list() {
    p = poly_attach(hexahedron(), hexahedron(), f1=[0,1], f2=0);
    assert(poly_valid(p, "closed"), "attach cube+cube two-face list should be closed valid");
    assert_int_eq(len(poly_verts(p)), 16, "attach cube+cube two-face verts");
    assert_int_eq(len(poly_faces(p)), 14, "attach cube+cube two-face faces");
    assert_int_eq(len(poly_edges(p)), 28, "attach cube+cube two-face edges");
}

module test_poly_attach__mirror_flag_changes_result() {
    p2 = _tetra_skew_apex();
    q0 = poly_attach(octahedron(), p2, f1=0, f2=0, mirror=false);
    q1 = poly_attach(octahedron(), p2, f1=0, f2=0, mirror=true);
    assert(poly_valid(q0, "closed"), "attach mirror=false closed");
    assert(poly_valid(q1, "closed"), "attach mirror=true closed");
    // A skew (non-symmetric) p2 should produce a different glued result when mirrored.
    dmax = max([for (i = [0:1:len(poly_verts(q0))-1]) norm(poly_verts(q0)[i] - poly_verts(q1)[i])]);
    assert(dmax > 1e-6, str("attach mirror flag should change geometry dmax=", dmax));
}

module run_TestAttach() {
    test_poly_attach__tetra_to_tetra_counts();
    test_poly_attach__cube_to_cube_counts();
    test_poly_attach__rotate_step_modulo();
    test_poly_attach__fit_edge_matches_unscaled_reference();
    test_poly_attach__f1_scalar_shorthand_equals_singleton_list();
    test_poly_attach__multi_face_list();
    test_poly_attach__mirror_flag_changes_result();
}

run_TestAttach();
