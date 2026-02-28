use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/core/validate.scad>

EPS = 1e-8;

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}

function _count_faces_of_size(poly, k) =
    sum([for (f = poly_faces(poly)) (len(f) == k) ? 1 : 0]);

function _edge_rel_spread(poly) =
    let(
        verts = poly_verts(poly),
        edges = poly_edges(poly),
        lens = [for (e = edges) norm(verts[e[0]] - verts[e[1]])],
        avg = (len(lens) == 0) ? 0 : sum(lens) / len(lens)
    )
    (len(lens) == 0 || avg == 0) ? 0 : ((max(lens) - min(lens)) / avg);

module test_poly_prism__counts_and_validity() {
    for (n = [3:1:6]) {
        p = poly_prism(n);
        assert_int_eq(len(poly_verts(p)), 2*n, str("prism n=", n, " verts"));
        assert_int_eq(len(poly_faces(p)), n + 2, str("prism n=", n, " faces"));
        assert_int_eq(len(poly_edges(p)), 3*n, str("prism n=", n, " edges"));
        if (n == 4) {
            // Cube case: caps and sides are all quads.
            assert_int_eq(_count_faces_of_size(p, 4), n + 2, "prism n=4 all faces are quads");
        } else {
            assert_int_eq(_count_faces_of_size(p, n), 2, str("prism n=", n, " two n-gons"));
            assert_int_eq(_count_faces_of_size(p, 4), n, str("prism n=", n, " n quads"));
        }
        assert(poly_valid(p, "closed"), str("prism n=", n, " closed valid"));
        assert(_edge_rel_spread(p) < 1e-10, str("prism n=", n, " uniform edges"));
    }
}

module test_poly_antiprism__counts_and_validity() {
    for (n = [3:1:6]) {
        p = poly_antiprism(n);
        assert_int_eq(len(poly_verts(p)), 2*n, str("antiprism n=", n, " verts"));
        assert_int_eq(len(poly_faces(p)), 2*n + 2, str("antiprism n=", n, " faces"));
        assert_int_eq(len(poly_edges(p)), 4*n, str("antiprism n=", n, " edges"));
        if (n == 3) {
            // Tetrahedral case: caps and side faces are all triangles.
            assert_int_eq(_count_faces_of_size(p, 3), 2*n + 2, "antiprism n=3 all faces are triangles");
        } else {
            assert_int_eq(_count_faces_of_size(p, n), 2, str("antiprism n=", n, " two n-gons"));
            assert_int_eq(_count_faces_of_size(p, 3), 2*n, str("antiprism n=", n, " 2n triangles"));
        }
        assert(poly_valid(p, "closed"), str("antiprism n=", n, " closed valid"));
        assert(_edge_rel_spread(p) < 1e-10, str("antiprism n=", n, " uniform edges"));
    }
}

module test_poly_antiprism__angle_and_height_controls() {
    p0 = poly_antiprism(5, edge=1, angle=0);
    p1 = poly_antiprism(5, edge=1, angle=12);
    p2 = poly_antiprism(5, edge=1, angle=0, height=0.7, height_scale=2);
    assert(poly_valid(p0, "closed"), "antiprism angle=0 valid");
    assert(poly_valid(p1, "closed"), "antiprism angle!=0 valid");
    assert(poly_valid(p2, "closed"), "antiprism explicit height valid");
    // angle should change geometry
    assert(abs(poly_verts(p0)[5][0] - poly_verts(p1)[5][0]) > 1e-6, "antiprism angle affects top ring");
    // explicit height*scale should be honored
    zspan = max([for (v = poly_verts(p2)) v[2]]) - min([for (v = poly_verts(p2)) v[2]]);
    assert_near(zspan, 1.4, 1e-9, "antiprism height*height_scale");
}

module test_poly_prism__height_controls() {
    p = poly_prism(6, edge=1, height=0.8, height_scale=1.25);
    assert(poly_valid(p, "closed"), "prism explicit height valid");
    zspan = max([for (v = poly_verts(p)) v[2]]) - min([for (v = poly_verts(p)) v[2]]);
    assert_near(zspan, 1.0, 1e-9, "prism height*height_scale");
}

module test_poly_prism__height_scale_keeps_nominal_edge_normalization() {
    p1 = poly_prism(6, edge=1, height_scale=1);
    p2 = poly_prism(6, edge=1, height_scale=2);
    assert_near(poly_e_over_ir(p1), 1, 1e-10, "prism n=6 baseline e_over_ir");
    assert_near(poly_e_over_ir(p2), 1, 1e-10, "prism n=6 tall e_over_ir");
}

module test_poly_antiprism__height_scale_keeps_nominal_edge_normalization() {
    p1 = poly_antiprism(6, edge=1, angle=0, height_scale=1);
    p2 = poly_antiprism(6, edge=1, angle=0, height_scale=2);
    assert_near(poly_e_over_ir(p1), poly_e_over_ir(p2), 1e-10, "antiprism n=6 e_over_ir invariant vs height_scale");
}

module run_TestPrisms() {
    test_poly_prism__counts_and_validity();
    test_poly_antiprism__counts_and_validity();
    test_poly_antiprism__angle_and_height_controls();
    test_poly_prism__height_controls();
    test_poly_prism__height_scale_keeps_nominal_edge_normalization();
    test_poly_antiprism__height_scale_keeps_nominal_edge_normalization();
}
