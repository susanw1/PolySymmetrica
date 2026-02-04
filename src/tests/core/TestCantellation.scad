use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../testing_util.scad>

EPS = 1e-7;

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

function _count_faces_of_size(poly, k) =
    sum([ for (f = poly_faces(poly)) (len(f)==k) ? 1 : 0 ]);

module test_poly_cantellate__tetra_counts() {
    p = tetrahedron();
    q = poly_cantellate(p, 0.2);
    expected_faces = len(poly_faces(p)) + len(_ps_edges_from_faces(poly_faces(p))) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantellate faces count");
    assert_int_eq(_count_faces_of_size(q, 3), 8, "cantellate tetra: 8 triangles");
    assert_int_eq(_count_faces_of_size(q, 4), 6, "cantellate tetra: 6 quads");
}

module test_poly_cantellate__cube_counts() {
    p = hexahedron();
    q = poly_cantellate(p, 0.2);
    expected_faces = len(poly_faces(p)) + len(_ps_edges_from_faces(poly_faces(p))) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantellate cube faces count");
    assert_int_eq(_count_faces_of_size(q, 3), 8, "cantellate cube: 8 triangles");
    assert_int_eq(_count_faces_of_size(q, 4), 18, "cantellate cube: 18 quads");
}

module run_TestCantellation() {
    test_poly_cantellate__tetra_counts();
    test_poly_cantellate__cube_counts();
}

run_TestCantellation();
