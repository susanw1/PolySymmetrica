use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/cleanup.scad>
use <../../polysymmetrica/core/params.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/platonics_all.scad>
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

function _count_faces_of_size(poly, k) =
    sum([for (f = poly_faces(poly)) (len(f) == k) ? 1 : 0]);

function _cube_with_cycle_noise() =
    let(
        p = hexahedron(),
        v = poly_verts(p),
        f = poly_faces(p),
        f0 = f[0],
        f0_noisy = [f0[0], f0[1], f0[1], f0[2], f0[3], f0[0]],
        f_new = concat([f0_noisy], [for (i = [1:1:len(f)-1]) f[i]])
    )
    [v, f_new, poly_e_over_ir(p)];

function _cube_with_degenerate_extra_face() =
    let(
        p = hexahedron(),
        v = poly_verts(p),
        f = poly_faces(p)
    )
    [v, concat(f, [[0,1,0]]), poly_e_over_ir(p)];

function _pyramid_with_nonplanar_base() =
    let(
        v = [[1,1,0],[-1,1,0],[-1,-1,0.2],[1,-1,0],[0,0,1]],
        f = [[0,1,2,3],[0,4,1],[1,4,2],[2,4,3],[3,4,0]]
    )
    [v, f, 1];

function _tetra_with_duplicate_vertex() =
    let(
        p = _tetra_poly(),
        v0 = poly_verts(p),
        v = concat(v0, [v0[0]]),
        f = [[4,1,2],[4,3,1],[4,2,3],[1,3,2]]
    )
    [v, f, poly_e_over_ir(p)];


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

// --- prefix offsets / face edge helpers ---
module test__ps_prefix_offsets__basic() {
    counts = [3, 2, 4];
    offs = _ps_prefix_offsets(counts, [0]);
    assert(offs == [0,3,5,9], "prefix offsets");
}

module test__ps_face_offsets__basic() {
    faces = [[0,1,2],[0,1,2,3]];
    offs = _ps_face_offsets(faces);
    assert(offs == [0,3,7], "face offsets");
}

module test__ps_face_edge_offsets__basic() {
    faces = [[0,1,2],[0,1,2,3]];
    offs = _ps_face_edge_offsets(faces);
    assert(offs == [0,6,14], "face edge offsets");
}

module test__ps_face_edge_site__basic() {
    assert(_ps_face_edge_site(10, 2, false) == 14, "edge site base+2k");
    assert(_ps_face_edge_site(10, 2, true) == 15, "edge site base+2k+1");
}

module test__ps_face_edge_sites_for_face_edge__square() {
    faces = [[0,1,2,3]];
    base = 0;
    s = _ps_face_edge_sites_for_face_edge(faces, 0, 1, 2, base);
    assert(s == [2,3], "edge sites for v1->v2");
    s2 = _ps_face_edge_sites_for_face_edge(faces, 0, 2, 1, base);
    assert(s2 == [3,2], "edge sites for reversed edge");
}

module test__ps_project_edge_pts_for_face_edge__square_xy() {
    verts = [[0,0,0],[1,0,0],[1,1,0],[0,1,0]];
    faces = [[0,1,2,3]];
    edges = _ps_edges_from_faces(faces);
    edge_pts = [
        for (ei = [0:1:len(edges)-1])
            let(a = edges[ei][0], b = edges[ei][1], A = verts[a], B = verts[b])
                [A + 0.25 * (B - A), B + 0.25 * (A - B)]
    ];
    n_f = [0,0,1];
    p0 = [0,0,0];
    pair = _ps_project_edge_pts_for_face_edge(verts, edges, edge_pts, n_f, p0, 0, 1);
    // ordered along +X (edge 0->1)
    assert(pair[0][0] < pair[1][0], "project edge pts order");
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

// --- params-overrides helpers ---
module test_params_overrides__lookup_and_count() {
    pbf = [
        ["face", "family", 0, ["df", 0.1], ["angle", 5], ["df", 0.2]],
        ["face", "family", 0, ["df", undef], ["angle", 7]],
        ["face", "family", 1, ["df", 0.3]],
        ["vert", "family", 0, ["c", 0.4], ["de", 0.5]],
        ["vert", "family", 0, ["c", 0.6]]
    ];

    assert_int_eq(ps_params_count_kind(pbf, "face"), 3, "params count face rows");
    assert_int_eq(ps_params_count_kind(pbf, "vert"), 2, "params count vert rows");
    assert_int_eq(ps_params_count_kind(undef, "face"), 0, "params count undef");

    assert_near(ps_params_get(pbf, "face", "df", 0, 0), 0.2, EPS, "params get face0 df");
    assert_near(ps_params_get(pbf, "face", "angle", 0, 0), 7, EPS, "params get face0 angle");
    assert_near(ps_params_get(pbf, "face", "df", 1, 1), 0.3, EPS, "params get face1 df");
    assert_near(ps_params_get(pbf, "vert", "c", 0, 0), 0.6, EPS, "params get vert0 c");
    assert_near(ps_params_get(pbf, "vert", "de", 0, 0), 0.5, EPS, "params get vert0 de");
    assert(is_undef(ps_params_get(pbf, "face", "df", 2, 2)), "params get missing family -> undef");
    assert(is_undef(ps_params_get(pbf, "face", "missing", 0, 0)), "params get missing key -> undef");
}

module test_params_overrides__compile_dense_arrays() {
    pbf = [
        ["face", "family", 0, ["df", 0.2], ["angle", 7]],
        ["face", "family", 2, ["angle", 9]],
        ["vert", "family", 1, ["c", 0.6], ["de", 0.7]]
    ];
    compiled = ps_params_compile_specs(pbf, [
        ["face", "df", 4, [0,1,2,3]],
        ["face", "angle", 4, [0,1,2,3]],
        ["vert", "c", 3, [0,1,2]],
        ["vert", "de", 3, [0,1,2]]
    ]);
    face_df = compiled[0];
    face_angle = compiled[1];
    vert_c = compiled[2];
    vert_de = compiled[3];

    assert_int_eq(len(face_df), 4, "compiled face_df len");
    assert_int_eq(len(face_angle), 4, "compiled face_angle len");
    assert_int_eq(len(vert_c), 3, "compiled vert_c len");
    assert_int_eq(len(vert_de), 3, "compiled vert_de len");

    assert_near(face_df[0], 0.2, EPS, "compiled face_df[0]");
    assert(is_undef(face_df[1]), "compiled face_df[1] undef");
    assert(is_undef(face_df[2]), "compiled face_df[2] undef");
    assert(is_undef(face_df[3]), "compiled face_df[3] undef");

    assert_near(face_angle[0], 7, EPS, "compiled face_angle[0]");
    assert(is_undef(face_angle[1]), "compiled face_angle[1] undef");
    assert_near(face_angle[2], 9, EPS, "compiled face_angle[2]");
    assert(is_undef(face_angle[3]), "compiled face_angle[3] undef");

    assert(is_undef(vert_c[0]), "compiled vert_c[0] undef");
    assert_near(vert_c[1], 0.6, EPS, "compiled vert_c[1]");
    assert(is_undef(vert_c[2]), "compiled vert_c[2] undef");

    assert(is_undef(vert_de[0]), "compiled vert_de[0] undef");
    assert_near(vert_de[1], 0.7, EPS, "compiled vert_de[1]");
    assert(is_undef(vert_de[2]), "compiled vert_de[2] undef");

    assert_near(ps_compiled_param_get(face_angle, 2), 9, EPS, "compiled get in-bounds");
    assert(is_undef(ps_compiled_param_get(face_angle, -1)), "compiled get negative -> undef");
    assert(is_undef(ps_compiled_param_get(face_angle, 4)), "compiled get out-of-bounds -> undef");
    assert(is_undef(ps_compiled_param_get(undef, 0)), "compiled get undef arr -> undef");
}

module test_params_overrides__compile_single_key() {
    pbf = [
        ["face", "family", 0, ["df", 0.2]],
        ["face", "family", 2, ["df", 0.9]]
    ];
    arr = ps_params_compile_key(pbf, "face", "df", 4, [0,1,2,3]);
    assert_int_eq(len(arr), 4, "compile key len");
    assert_near(arr[0], 0.2, EPS, "compile key face0");
    assert(is_undef(arr[1]), "compile key face1 undef");
    assert_near(arr[2], 0.9, EPS, "compile key face2");
    assert(is_undef(arr[3]), "compile key face3 undef");
}

module test_params__selector_precedence_and_compile() {
    p = [
        ["face", "all", ["df", 0.1], ["angle", 5]],
        ["face", "family", 1, ["df", 0.2]],
        ["face", "id", [3, 4], ["df", 0.4]],
        ["face", "id", 4, ["angle", 11]]
    ];

    // id > family > all
    assert_near(ps_params_get(p, "face", "df", 2, 1), 0.2, EPS, "params precedence family");
    assert_near(ps_params_get(p, "face", "df", 3, 1), 0.4, EPS, "params precedence id over family");
    assert_near(ps_params_get(p, "face", "df", 9, 0), 0.1, EPS, "params precedence all fallback");
    assert_near(ps_params_get(p, "face", "angle", 4, 1), 11, EPS, "params id scalar selector");
    assert_near(ps_params_get(p, "face", "angle", 1, 1), 5, EPS, "params all fallback angle");

    fam_ids = [1, 1, 1, 1, 0];
    arr = ps_params_compile_key(p, "face", "df", 5, fam_ids);
    assert_int_eq(len(arr), 5, "compile key element len");
    assert_near(arr[0], 0.2, EPS, "compile key elem0");
    assert_near(arr[1], 0.2, EPS, "compile key elem1");
    assert_near(arr[2], 0.2, EPS, "compile key elem2");
    assert_near(arr[3], 0.4, EPS, "compile key elem3 id override");
    assert_near(arr[4], 0.4, EPS, "compile key elem4 id override");
}

// --- cleanup ---
module test_poly_cleanup__normalizes_face_cycles() {
    p = _cube_with_cycle_noise();
    q = poly_cleanup(p, eps=1e-8, fix_winding=true, drop_degenerate=true);
    f0 = poly_faces(q)[0];
    assert_int_eq(len(f0), 4, "cleanup normalized noisy face to quad");
    assert(poly_valid(q, "closed"), "cleanup cycle normalization should produce closed-valid poly");
}

module test_poly_cleanup__drops_degenerate_faces() {
    p = _cube_with_degenerate_extra_face();
    q = poly_cleanup(p, eps=1e-8, fix_winding=true, drop_degenerate=true);
    assert_int_eq(len(poly_faces(q)), 6, "cleanup dropped degenerate extra face");
    assert(poly_valid(q, "closed"), "cleanup dropped degenerate face should remain closed-valid");
}

module test_poly_cleanup__triangulates_nonplanar_faces() {
    p = _pyramid_with_nonplanar_base();
    q = poly_cleanup(p, eps=1e-8, triangulate_nonplanar=true, fix_winding=true);
    assert_int_eq(len(poly_faces(q)), 6, "cleanup triangulates nonplanar quad base into two triangles");
    assert_int_eq(_count_faces_of_size(q, 4), 0, "cleanup triangulated nonplanar quads");
    assert(poly_valid(q, "struct"), "cleanup triangulated output should be structurally valid");
}

module test_poly_cleanup__merges_and_compacts_vertices() {
    p = _tetra_with_duplicate_vertex();
    q = poly_cleanup(p, eps=1e-8, merge_vertices=true, remove_unreferenced=true, fix_winding=true);
    assert_int_eq(len(poly_verts(q)), 4, "cleanup merged duplicate vertex and compacted unreferenced");
    assert(poly_valid(q, "closed"), "cleanup merge/compact should yield closed-valid tetra");
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
    test_params_overrides__lookup_and_count();
    test_params_overrides__compile_dense_arrays();
    test_params_overrides__compile_single_key();
    test_params__selector_precedence_and_compile();
    test_poly_cleanup__normalizes_face_cycles();
    test_poly_cleanup__drops_degenerate_faces();
    test_poly_cleanup__triangulates_nonplanar_faces();
    test_poly_cleanup__merges_and_compacts_vertices();
}

run_TestFuncs();
