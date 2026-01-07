use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/duals.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/regular_all.scad>
use <../testing_util.scad>

EPS = 1e-7;

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}
module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

function _count_faces_of_size(poly, k) =
    sum([ for (f = poly_faces(poly)) (len(f)==k) ? 1 : 0 ]);


// internal helpers: _ps_edge_point_near, equality/find/unique/remap

module test__ps_edge_point_near__picks_correct_end() {
    p=_tetra_poly();
    verts=poly_verts(p);
    faces=poly_faces(p);
    edges=_ps_edges_from_faces(faces);

    t=0.2;
    edge_pts = [
        for (ei=[0:len(edges)-1])
            let(a=edges[ei][0], b=edges[ei][1], A=verts[a], B=verts[b])
            [A + t*(B-A), B + t*(A-B)]
    ];

    // pick edge (0,1) near 0 should equal first point of that edge entry
    P0 = _ps_edge_point_near(edges, edge_pts, 0, 1, 0);
    ei = find_edge_index(edges, 0, 1);
    assert(norm(P0 - edge_pts[ei][0]) < 1e-12, "near 0 uses [0]");
}

// point eq / find / unique
module test__ps_unique_points__dedups_with_eps() {
    pts = [[0,0,0],[0,0,0.0],[1,0,0],[1,0,0.0000000001]];
    uniq = _ps_unique_points(pts, 1e-6);
    assert_int_eq(len(uniq), 2, "dedup to 2");
}

module test__ps_face_points_to_indices__maps() {
    uniq = [[0,0,0],[1,0,0],[0,1,0]];
    face = [[0,0,0],[0,1,0],[1,0,0]];
    idx = _ps_face_points_to_indices(uniq, face, 1e-9);
    assert(idx == [0,2,1], "maps indices");
}


// truncation: default t, corner t
module test__ps_corner_t__equilateral_is_one_thirdish() {
    verts = [[0,0,0],[1,0,0],[0.5,sqrt(3)/2,0]];
    face = [0,1,2];
    t = _ps_corner_t(verts, face, 0);
    assert_near(t, 1/3, 1e-6, "equilateral corner t");
}

module test__ps_truncate_default_t__tetra_one_thirdish() {
    p=_tetra_poly();
    t=_ps_truncate_default_t(p, tol=1e-2, fallback=0.2);
    assert_near(t, 1/3, 1e-2, "default t tetra ~1/3");
}


// truncation main: combinatorics for t=1/3
module test_poly_truncate__tetra_counts_at_one_third() {
    p=_tetra_poly();
    q=poly_truncate(p, 1/3);

    assert_int_eq(len(poly_verts(q)), 12, "trunc tetra verts=12");
    assert_int_eq(len(poly_faces(q)), 8, "trunc tetra faces=8");

    sizes = [for (f=poly_faces(q)) len(f)];
    assert_int_eq(sum([for (s=sizes) if (s==3) 1]), 4, "4 tri faces");
    assert_int_eq(sum([for (s=sizes) if (s==6) 1]), 4, "4 hex faces");
}


// truncation t=0: counts preserved (geometry changes because of dedup path; still should match input)
module test_poly_truncate__t_zero_counts_preserved() {
    p=_tetra_poly();
    q=poly_truncate(p, 0);
    assert_int_eq(len(poly_verts(q)), len(poly_verts(p)), "t=0 verts");
    assert_int_eq(len(poly_faces(q)), len(poly_faces(p)), "t=0 faces");
    assert_facet_matches(p, q);
}


// truncation + dual: size relations
module test_poly_truncate_then_dual__counts_relations() {
    p=_tetra_poly();
    q=poly_truncate(p, 1/3);
    d=poly_dual(q);

    assert_int_eq(len(poly_verts(d)), len(poly_faces(q)), "dual verts==faces");
    assert_int_eq(len(poly_faces(d)), len(poly_verts(q)), "dual faces==verts");
}

module test_poly_rectify__tetra_counts() {
    p=_tetra_poly();
    q=poly_rectify(p);
    assert_poly_valid(q);

    // rectified tetrahedron is an octahedron
    assert(len(poly_verts(q)) == 6, "rectify tetra verts");
    assert(len(poly_faces(q)) == 8, "rectify tetra faces");
    assert(_count_faces_of_size(q,3) == 8, "rectify tetra: 8 triangles");
}

module test_poly_cantellate__tetra_counts() {
    p = _tetra_poly();
    edges = _ps_edges_from_faces(poly_faces(p));
    q = poly_cantellate(p, 0.2, 0.2);
    assert_poly_valid(q);
    assert_poly_valid_mode(q, "struct");

    expected_faces = len(poly_faces(p)) + len(edges) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantellate faces count");
}

module test_truncate__validity() {
    p = poly_truncate(tetrahedron(), 1/3);
    assert_poly_valid_mode(p, "closed");
}

module test_truncate__tetra_archimedean_counts() {
    p = poly_truncate(tetrahedron(), 1/3);
    assert_poly_valid(p);

    assert(len(poly_verts(p)) == 12, "trunc tetra verts");
    assert(len(poly_faces(p)) == 8,  "trunc tetra faces");
    assert(_count_faces_of_size(p,3) == 4, "trunc tetra: 4 triangles");
    assert(_count_faces_of_size(p,6) == 4, "trunc tetra: 4 hexagons");
}

module test_truncate__octa_archimedean_counts() {
    p = poly_truncate(octahedron(), 1/3);
    assert_poly_valid(p);

    // truncated octahedron: 24 verts, 14 faces = 8 hex + 6 squares
    assert(len(poly_verts(p)) == 24, "trunc octa verts");
    assert(len(poly_faces(p)) == 14, "trunc octa faces");
    assert(_count_faces_of_size(p,4) == 6, "trunc octa: 6 squares");
    assert(_count_faces_of_size(p,6) == 8, "trunc octa: 8 hexagons");
}

module test_truncate__icosa_archimedean_counts() {
    p = poly_truncate(icosahedron(), 1/3);
    assert_poly_valid(p);

    // truncated icosahedron: 60 verts, 32 faces = 12 pent + 20 hex
    assert(len(poly_verts(p)) == 60, "trunc icosa verts");
    assert(len(poly_faces(p)) == 32, "trunc icosa faces");
    assert(_count_faces_of_size(p,5) == 12, "trunc icosa: 12 pentagons");
    assert(_count_faces_of_size(p,6) == 20, "trunc icosa: 20 hexagons");
}

module test_truncate__hexa_archimedean_counts() {
    // If you don't have a hexahedron() model file, use dual(octa)
    c = poly_dual(octahedron());
    p = poly_truncate(c, 1/3);
    assert_poly_valid(p);

    // truncated cube: 24 verts, 14 faces = 8 tri + 6 oct
    assert(len(poly_verts(p)) == 24, "trunc cube verts");
    assert(len(poly_faces(p)) == 14, "trunc cube faces");
    assert(_count_faces_of_size(p,3) == 8, "trunc cube: 8 triangles");
    assert(_count_faces_of_size(p,8) == 6, "trunc cube: 6 octagons");
}

module test_truncate__dodeca_archimedean_counts() {
    // If you don't have dodecahedron() model file, use dual(icosa)
    d = poly_dual(icosahedron());
    p = poly_truncate(d, 1/3);
    assert_poly_valid(p);

    // truncated dodeca: 60 verts, 32 faces = 20 tri + 12 dec
    assert(len(poly_verts(p)) == 60, "trunc dodeca verts");
    assert(len(poly_faces(p)) == 32, "trunc dodeca faces");
    assert(_count_faces_of_size(p,3) == 20, "trunc dodeca: 20 triangles");
    assert(_count_faces_of_size(p,10) == 12, "trunc dodeca: 12 decagons");
}



// suite
module run_TestTruncation() {
    test__ps_edge_point_near__picks_correct_end();
    test__ps_unique_points__dedups_with_eps();
    test__ps_face_points_to_indices__maps();

    test__ps_corner_t__equilateral_is_one_thirdish();
    test__ps_truncate_default_t__tetra_one_thirdish();

    test_poly_truncate__tetra_counts_at_one_third();
    test_poly_truncate__t_zero_counts_preserved();
    test_poly_truncate_then_dual__counts_relations();
    test_poly_rectify__tetra_counts();
    test_poly_cantellate__tetra_counts();
    test_truncate__validity();
    
    test_truncate__tetra_archimedean_counts();
    test_truncate__octa_archimedean_counts();
    test_truncate__icosa_archimedean_counts();
    test_truncate__hexa_archimedean_counts();
    test_truncate__dodeca_archimedean_counts();
}

run_TestTruncation();
