use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/validate.scad>
use <../testing_util.scad>

EPS = 1e-9;

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}

module assert_vec3_near(v, w, eps=EPS, msg="") {
    assert(len(v)==3 && len(w)==3, str(msg, " not vec3"));
    for (i=[0:2]) assert_near(v[i], w[i], eps, str(msg, " i=", i));
}

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}


// --- poly accessors ---
module test_poly_accessors__basic() {
    p = [[[0,0,0],[1,0,0],[0,1,0],[0,0,1]], [[0,1,2],[0,1,3],[0,2,3],[1,2,3]], 2.0];
    assert_int_eq(len(poly_verts(p)), 4, "poly_verts");
    assert_int_eq(len(poly_faces(p)), 4, "poly_faces");
    assert_near(poly_e_over_ir(p), 2.0, 0, "poly_e_over_ir");
}


// --- ps_indices_in_range / ps_faces_valid ---
module test_all_indices_in_range__true_false() {
    assert(ps_indices_in_range([0,1,2], 3), "in range");
    assert(!ps_indices_in_range([0,1,3], 3), "out of range");
}

module test_all_faces_valid__true_false() {
    verts = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]];
    faces_ok  = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
    faces_bad = [[0,1,2],[0,1,3],[0,2,3],[1,2,9]];
    assert(ps_faces_valid(verts, faces_ok), "faces valid");
    assert(!ps_faces_valid(verts, faces_bad), "faces invalid");
}


// --- make_poly (positive path) ---
module test_make_poly__valid_tetra_computes_e_over_ir() {
    verts = [[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]];
    faces = [[0,1,2],[0,3,1],[0,2,3],[1,3,2]];
    p = make_poly(verts, faces);
    assert_int_eq(len(poly_verts(p)), 4, "verts");
    assert_int_eq(len(poly_faces(p)), 4, "faces");
    assert(poly_e_over_ir(p) > 0, "e_over_ir > 0");
}

// --- poly_fix_winding ---
module test_poly_fix_winding__repairs_edge_dirs() {
    verts = [[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]];
    // Two faces share edge (0,3) in the same direction -> winding mismatch.
    faces = [[0,2,1],[0,3,1],[0,3,2],[1,2,3]];
    p = [verts, faces, 1];
    assert(!poly_validate_winding(p), "winding mismatch");
    fixed = poly_fix_winding(p);
    assert(poly_validate_winding(fixed), "winding repaired");
}

// --- vector helpers ---
module test_v_ops__dot_cross_norm() {
    assert_near(v_dot([1,2,3],[4,5,6]), 32, 0, "v_dot");
    assert_vec3_near(v_cross([1,0,0],[0,1,0]), [0,0,1], EPS, "v_cross");
    assert_near(v_len([3,4,0]), 5, EPS, "v_len");
    assert_vec3_near(v_norm([0,0,0]), [0,0,0], 0, "v_norm zero");
    assert_vec3_near(v_norm([3,0,0]), [1,0,0], EPS, "v_norm axis");
}

// internal helpers: v_ordered
module test_v_ordered() {
    assert(v_ordered(1,2) == [1,2], "sorted 1,2");
    assert(v_ordered(2,1) == [1,2], "sorted 2,1");
}

// --- edge_equal ---
module test_edge_equal__ordered() {
    assert(edge_equal([0,1],[0,1]), "ordered equal");
    assert(!edge_equal([0,1],[1,0]), "reverse not equal");
}


// --- sum / v_sum ---
module test_sum__numbers() {
    assert_near(sum([1,2,3,4]), 10, 0, "sum ints");
    assert_near(sum([0.1,0.2,0.3]), 0.6, 1e-12, "sum floats");
}

module test_v_sum__vec3_list() {
    assert_vec3_near(v_sum([[1,2,3],[4,5,6]]), [5,7,9], EPS, "v_sum");
}

// --- _ps_ordered_pair ---
module test_ps_ordered_pair__basic() {
    assert(_ps_ordered_pair(1,2) == [1,2], "ordered 1,2");
    assert(_ps_ordered_pair(2,1) == [1,2], "ordered 2,1");
}

// --- _ps_reverse ---
module test_ps_reverse__edge_cases() {
    assert(_ps_reverse([]) == [], "reverse empty");
    assert(_ps_reverse([1]) == [1], "reverse singleton");
    assert(_ps_reverse([1,2,3]) == [3,2,1], "reverse list");
}

// --- _ps_list_contains / _ps_index_of ---
module test_ps_list_contains_index__basic() {
    assert(_ps_list_contains([1,2,3], 2), "contains 2");
    assert(!_ps_list_contains([1,2,3], 4), "not contains 4");
    assert_int_eq(_ps_index_of([1,2,3], 3), 2, "index of 3");
    assert_int_eq(_ps_index_of([1,2,3], 4), -1, "index missing");
}

// --- ps_point_eq ---
module test_ps_point_eq__eps() {
    assert(ps_point_eq([0,0,0], [0,0,1e-6], 1e-5), "point eq within eps");
    assert(!ps_point_eq([0,0,0], [0,0,1e-3], 1e-5), "point eq outside eps");
}


module test_find_edge_index__finds() {
    faces = poly_faces(_tetra_poly());
    edges = _ps_edges_from_faces(faces);
    i = ps_find_edge_index(edges, 0, 1);
    assert(i >= 0, "edge exists");
    assert(edges[i] == [0,1], "stored sorted");
}

// --- polygon helpers ---
module test_calc_edge_radius_roundtrip() {
    for (n=[3,4,5,6]) {
        R=2.5;
        e = ps_calc_edge(n, R);
        R2 = ps_calc_radius(n, e);
        assert_near(R2, R, 1e-9, str("roundtrip n=", n));
    }
}

module test_calc_edge__known_values() {
    assert_near(ps_calc_edge(3,1), sqrt(3), 1e-9, "triangle R=1");
    assert_near(ps_calc_edge(4,1), sqrt(2), 1e-9, "square R=1");
}


// --- ps_face_centroid / ps_face_normal ---
module test_face_centroid__triangle() {
    verts = [[0,0,0],[2,0,0],[0,2,0]];
    c = ps_face_centroid(verts, [0,1,2]);
    assert_vec3_near(c, [2/3, 2/3, 0], 1e-12, "centroid");
}

module test_face_normal__orientation() {
    verts = [[0,0,0],[1,0,0],[0,1,0]];
    n1 = ps_face_normal(verts, [0,1,2]);
    n2 = ps_face_normal(verts, [0,2,1]);
    assert_vec3_near(n1, [0,0,-1], 1e-12, "normal ccw");
    assert_vec3_near(n2, [0,0,1], 1e-12, "normal cw");
}


// --- poly_vertex_neighbor ---
module test_poly_vertex_neighbor__returns_incident_vertex() {
    verts = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]];
    faces = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
    p = make_poly(verts, faces, 1);
    n = poly_vertex_neighbor(p, 0);
    assert(n==1 || n==2 || n==3, "neighbor must be incident");
}


// --- _ps_edges_from_faces / ps_face_has_edge / ps_edge_faces_table ---
module test_edges_from_faces__tetra_counts() {
    faces = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
    edges = _ps_edges_from_faces(faces);
    assert_int_eq(len(edges), 6, "tetra edges=6");
}

module test_edges_from_faces__quad_edges() {
    faces = [[0,1,2,3]];
    edges = _ps_edges_from_faces(faces);
    assert_int_eq(len(edges), 4, "quad edges=4");
}

module test_face_has_edge__undirected_adjacent_only() {
    f = [0,1,2,3];
    assert(ps_face_has_edge(f, 1,2), "has (1,2)");
    assert(ps_face_has_edge(f, 2,1), "has (2,1)");
    assert(!ps_face_has_edge(f, 0,2), "no diagonal (0,2)");
}

module test_edge_faces_table__tetra_each_edge_two_faces() {
    faces = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
    edges = _ps_edges_from_faces(faces);
    t = ps_edge_faces_table(faces, edges);
    for (ei=[0:len(edges)-1]) assert_int_eq(len(t[ei]), 2, str("edge ", ei, " has 2 faces"));
}


// --- face-frame helpers ---
module test_poly_face_center__matches_centroid_for_triangle() {
    verts = [[0,0,0],[2,0,0],[0,2,0],[0,0,2]];
    faces = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
    p = make_poly(verts, faces, 1);
    c = poly_face_center(p, 0, 1);
    assert_vec3_near(c, [2/3,2/3,0], 1e-12, "face center");
}

module test_poly_face_axes__orthonormalish() {
    verts = [[0,0,0],[2,0,0],[0,2,0],[0,0,2]];
    faces = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
    p = make_poly(verts, faces, 1);

    ex = poly_face_ex(p, 0, 1);
    ey = poly_face_ey(p, 0, 1);
    ez = poly_face_ez(p, 0, 1);

    assert_near(v_len(ex), 1, 1e-9, "ex unit");
    assert_near(v_len(ez), 1, 1e-9, "ez unit");
    assert_near(v_dot(ex, ez), 0, 1e-8, "ex ⟂ ez");
    // ey depends on cross(ez,ex) so should be orthogonal
    assert_near(v_dot(ey, ez), 0, 1e-8, "ey ⟂ ez");
}


// --- ps_orient_face_outward / ps_orient_all_faces_outward ---
module test_orient_face_outward__makes_centroid_dot_normal_nonnegative() {
    verts = [[1,0,1],[0,1,1],[-1,0,1]];
    f_ccw = [0,1,2];
    f_cw  = [0,2,1];

    g1 = ps_orient_face_outward(verts, f_ccw);
    g2 = ps_orient_face_outward(verts, f_cw);

    // invariant we actually care about:
    c1 = ps_face_centroid(verts, g1);
    n1 = ps_face_normal(verts, g1);
    assert(v_dot(c1, n1) >= -1e-12, "oriented face has centroid·normal >= 0");

    c2 = ps_face_centroid(verts, g2);
    n2 = ps_face_normal(verts, g2);
    assert(v_dot(c2, n2) >= -1e-12, "flipped face has centroid·normal >= 0");
}

module test_orient_all_faces_outward__length_preserved() {
    verts = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]];
    faces = [[0,2,1],[0,1,3],[0,3,2],[1,2,3]];
    out = ps_orient_all_faces_outward(verts, faces);
    assert_int_eq(len(out), len(faces), "face count preserved");
}

module test_ps_sort__numbers() {
    v = [3,1,4,1,5,9,2];
    s = _ps_sort(v);
    assert(s == [1,1,2,3,4,5,9], "sort ints");
}

module test_ps_sort__floats() {
    v = [0.3, 0.1, 0.2];
    s = _ps_sort(v);
    assert_near(s[0], 0.1, 1e-12, "sort floats 0");
    assert_near(s[1], 0.2, 1e-12, "sort floats 1");
    assert_near(s[2], 0.3, 1e-12, "sort floats 2");
}

module test_ps_sort__empty() {
    s = _ps_sort([]);
    assert_int_eq(len(s), 0, "sort empty");
}

// --- _ps_solve3 ---
module test_ps_solve3__identity_and_scale() {
    m1 = [[1,0,0],[0,1,0],[0,0,1]];
    b1 = [1,2,3];
    x1 = _ps_solve3(m1, b1);
    assert_vec3_near(x1, b1, EPS, "solve identity");

    m2 = [[2,0,0],[0,3,0],[0,0,4]];
    b2 = [4,9,8];
    x2 = _ps_solve3(m2, b2);
    assert_vec3_near(x2, [2,3,2], EPS, "solve scale");
}

// --- ps_frame_matrix ---
module test_ps_frame_matrix__basic_layout() {
    m = ps_frame_matrix([1,2,3], [1,0,0], [0,1,0], [0,0,1]);
    assert_vec3_near([m[0][0], m[0][1], m[0][2]], [1,0,0], EPS, "row0 xyz");
    assert_near(m[0][3], 1, EPS, "row0 w");
    assert_vec3_near([m[1][0], m[1][1], m[1][2]], [0,1,0], EPS, "row1 xyz");
    assert_near(m[1][3], 2, EPS, "row1 w");
    assert_vec3_near([m[2][0], m[2][1], m[2][2]], [0,0,1], EPS, "row2 xyz");
    assert_near(m[2][3], 3, EPS, "row2 w");
}


// ---- suite ----
module run_TestFuncs() {
    test_poly_accessors__basic();

    test_all_indices_in_range__true_false();
    test_all_faces_valid__true_false();
    test_make_poly__valid_tetra_computes_e_over_ir();
    test_poly_fix_winding__repairs_edge_dirs();

    test_v_ops__dot_cross_norm();
//    test_v_ordered();
    test_edge_equal__ordered();
    test_sum__numbers();
    test_v_sum__vec3_list();
    test_ps_ordered_pair__basic();
    test_ps_reverse__edge_cases();
    test_ps_list_contains_index__basic();
    test_ps_point_eq__eps();
    test_find_edge_index__finds();

    test_calc_edge_radius_roundtrip();
    test_calc_edge__known_values();

    test_face_centroid__triangle();
    test_face_normal__orientation();

    test_poly_vertex_neighbor__returns_incident_vertex();

    test_edges_from_faces__tetra_counts();
    test_edges_from_faces__quad_edges();
    test_face_has_edge__undirected_adjacent_only();
    test_edge_faces_table__tetra_each_edge_two_faces();

    test_poly_face_center__matches_centroid_for_triangle();
    test_poly_face_axes__orthonormalish();

    test_orient_face_outward__makes_centroid_dot_normal_nonnegative();
    test_orient_all_faces_outward__length_preserved();
    test_ps_sort__numbers();
    test_ps_sort__floats();
    test_ps_sort__empty();
    test_ps_solve3__identity_and_scale();
    test_ps_frame_matrix__basic_layout();
}

run_TestFuncs();
