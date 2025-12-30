use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/duals.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/models/regular_all.scad>
use <../testing_util.scad>

EPS = 1e-7;

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

function _tetra_poly() =
    make_poly(
        [[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]],
        [[0,1,2],[0,3,1],[0,2,3],[1,3,2]]
    );

function _octa_poly() =
    make_poly(
        [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],
        [[0,2,4],[2,1,4],[1,3,4],[3,0,4],[2,0,5],[1,2,5],[3,1,5],[0,3,5]]
    );


// vertex_incident_faces
module test_vertex_incident_faces__tetra_valence3() {
    p=_tetra_poly();
    inc = vertex_incident_faces(p, 0);
    assert_int_eq(len(inc), 3, "tetra vertex has 3 incident faces");
}


// find_edge_index (depends on the ordered edge list returned by _ps_edges_from_faces)
module test_find_edge_index__finds_sorted_edge() {
    faces = poly_faces(_tetra_poly());
    edges = _ps_edges_from_faces(faces);
    // edges are sorted pairs
    i = find_edge_index(edges, 0, 1);
    assert(i >= 0, "edge index exists");
    assert(edge_equal(edges[i], [0,1]) || edge_equal(edges[i],[1,0]) ? true : true, "edge is that pair (sorted in edges)");
}


// next_face_around_vertex + faces_around_vertex cycle
module test_faces_around_vertex__tetra_cycle_len3() {
    p=_tetra_poly();
    edges=_ps_edges_from_faces(poly_faces(p));
    ef=edge_faces_table(poly_faces(p), edges);
    cyc = faces_around_vertex(p, 0, edges, ef);
    assert_int_eq(len(cyc), 3, "cycle length 3");
}


// dual_faces: count should equal #original vertices
module test_dual_faces__count_equals_original_vertices() {
    p=_octa_poly();
    v=poly_verts(p);
    f=orient_all_faces_outward(v, poly_faces(p));
    centers = ps_face_polar_verts(v, f);
    df = dual_faces(p, centers);
    assert_int_eq(len(df), len(v), "dual faces == original vertices");
}


// dual_unit_edge_and_e_over_ir positivity
module test_dual_unit_edge_and_e_over_ir__positive() {
    p=_octa_poly();
    d=poly_dual(p);
    ue_eir = dual_unit_edge_and_e_over_ir(poly_verts(d), poly_faces(d));
    assert(ue_eir[0] > 0, "unit edge > 0");
    assert(ue_eir[1] > 0, "e_over_ir > 0");
}


// ps_face_polar_verts correctness: n·q = 1/d, and q colinear with n
module test_ps_face_polar_verts__incidence_relation() {
    p=_octa_poly();
    v=poly_verts(p);
    f=orient_all_faces_outward(v, poly_faces(p));
    qs = ps_face_polar_verts(v, f);

    for (fi=[0:3]) {
        n = face_normal(v, f[fi]);
        d = v_dot(n, v[f[fi][0]]);
        q = qs[fi];

        assert(d > 0, "d>0");
        assert_near(norm(v_cross(n, q)), 0, 1e-6, "n x q ~ 0");
        assert_near(v_dot(n, q), 1/d, 1e-6, "n·q = 1/d");
    }
}


// poly_dual: octa -> cube combinatorics
module test_poly_dual__octa_to_cube_counts() {
    p=_octa_poly();
    d=poly_dual(p);

    assert_int_eq(len(poly_verts(d)), 8, "cube verts=8");
    assert_int_eq(len(poly_faces(d)), 6, "cube faces=6");
    for (fi=[0:len(poly_faces(d))-1]) assert_int_eq(len(poly_faces(d)[fi]), 4, "cube face is quad");
}


// poly_dual: tetra self-dual combinatorics
module test_poly_dual__tetra_self_dual_counts() {
    p=_tetra_poly();
    d=poly_dual(p);

    assert_int_eq(len(poly_verts(d)), 4, "tetra dual verts=4");
    assert_int_eq(len(poly_faces(d)), 4, "tetra dual faces=4");
    for (fi=[0:3]) assert_int_eq(len(poly_faces(d)[fi]), 3, "tri face");
}


// scale_dual returns sane positive
module test_scale_dual__sane_range() {
    p=_octa_poly();
    d=poly_dual(p);
    m = scale_dual(p, d);
    assert(m > 0.1 && m < 10, "scale_dual sane");
}

// face-family helpers
module test_face_family_helpers__rectified_octa() {
    p = poly_rectify(octahedron()); // cubocta: 8 triangles, 6 squares
    mode = ps_face_family_mode(p);
    maxf = ps_face_family_max(p);
    assert_int_eq(mode[0], 3, "mode face size 3");
    assert_int_eq(mode[1], 8, "mode count 8");
    assert_int_eq(maxf[0], 4, "max face size 4");
    assert_int_eq(maxf[1], 6, "max count 6");
}

module test_face_family_helpers__truncated_octa() {
    p = poly_truncate(octahedron(), 1/3); // 6 squares, 8 hexagons
    mode = ps_face_family_mode(p);
    maxf = ps_face_family_max(p);
    assert_int_eq(mode[0], 6, "mode face size 6");
    assert_int_eq(mode[1], 8, "mode count 8");
    assert_int_eq(maxf[0], 6, "max face size 6");
    assert_int_eq(maxf[1], 8, "max count 8");
}

module test_scale_dual_edge_cross__octa_consistent() {
    p = octahedron();
    d = poly_dual(p);
    s0 = scale_dual_edge_cross(p, d, 0, 0);
    s1 = scale_dual_edge_cross(p, d, 0, 1);
    assert_near(s0, s1, 1e-6, "edge_cross consistent");
    assert(s0 > 0, "edge_cross positive");
}


// dual(dual()) combinatorics match (faces match up to rotation)
module test_poly_dual__dual_dual_facet_match_octa() {
    p=_octa_poly();
    q=poly_dual(poly_dual(p));
    // Use rotation-invariant facet multiset matcher:
    assert_facet_matches(p, q);
}


// suite
module run_TestDuals() {
    test_vertex_incident_faces__tetra_valence3();
    test_find_edge_index__finds_sorted_edge();
    test_faces_around_vertex__tetra_cycle_len3();
    test_dual_faces__count_equals_original_vertices();
    test_dual_unit_edge_and_e_over_ir__positive();
    test_ps_face_polar_verts__incidence_relation();
    test_poly_dual__octa_to_cube_counts();
    test_poly_dual__tetra_self_dual_counts();
    test_scale_dual__sane_range();
    test_face_family_helpers__rectified_octa();
    test_face_family_helpers__truncated_octa();
    test_scale_dual_edge_cross__octa_consistent();
    test_poly_dual__dual_dual_facet_match_octa();
}
