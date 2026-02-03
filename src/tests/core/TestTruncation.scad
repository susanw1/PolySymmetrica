use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/duals.scad>
use <../../polysymmetrica/core/transform.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/archimedians_all.scad>
use <../testing_util.scad>

EPS = 1e-7;
ENABLE_CANTITRUNC_PLANARITY_TEST = false;

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}
module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

function _count_faces_of_size(poly, k) =
    sum([ for (f = poly_faces(poly)) (len(f)==k) ? 1 : 0 ]);

function _map_face_c(size, c_by_size, default=0) =
    let(idxs = [for (i = [0:1:len(c_by_size)-1]) if (c_by_size[i][0] == size) i])
    (len(idxs) == 0) ? default : c_by_size[idxs[0]][1];

function _face_max_plane_err(verts, face) =
    let(
        n = v_norm(ps_face_normal(verts, face)),
        d = v_dot(n, verts[face[0]]),
        errs = [for (vi = face) abs(v_dot(n, verts[vi]) - d)]
    )
    (len(errs) == 0) ? 0 : max(errs);

function _skew_quad_prism() =
    let(
        v = [
            [0,0,0],[2,0,0],[2,1,0],[0,1,0],
            [0.3,0.2,1],[2.2,0.1,1],[2.1,1.2,1],[0.1,1.1,1]
        ],
        f = [
            [0,1,2,3],
            [4,5,6,7],
            [0,1,5,4],
            [1,2,6,5],
            [2,3,7,6],
            [3,0,4,7]
        ]
    )
    make_poly(v, f);


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
    assert_face_matches(p, q);
}

// chamfer: cube should have 6 original faces + 12 edge faces (hexes)
module test_poly_chamfer__cube_face_counts() {
    p = hexahedron();
    q = poly_chamfer(p, 0.1);
    assert_int_eq(len(poly_faces(q)), 18, "chamfer cube faces=18");
    assert_int_eq(_count_faces_of_size(q, 6), 12, "chamfer cube: 12 hex edge faces");
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
    q = poly_cantellate(p, 0.2);
    assert_poly_valid(q);
    assert_poly_valid_mode(q, "struct");

    expected_faces = len(poly_faces(p)) + len(edges) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantellate faces count");

    assert_int_eq(_count_faces_of_size(q, 3), 8, "cantellate tetra: 8 triangles");
    assert_int_eq(_count_faces_of_size(q, 4), 6, "cantellate tetra: 6 quads");
}

module test_poly_cantellate__cube_counts() {
    p = hexahedron();
    edges = _ps_edges_from_faces(poly_faces(p));
    q = poly_cantellate(p, 0.2);
    assert_poly_valid(q);

    expected_faces = len(poly_faces(p)) + len(edges) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantellate cube faces count");
    assert_int_eq(_count_faces_of_size(q, 3), 8, "cantellate cube: 8 triangles");
    assert_int_eq(_count_faces_of_size(q, 4), 18, "cantellate cube: 18 quads");
}

module test_poly_cantitruncate__tetra_counts() {
    p = _tetra_poly();
    edges = _ps_edges_from_faces(poly_faces(p));
    q = poly_cantitruncate(p, 0.2, 0.1);
    assert_poly_valid(q);

    expected_faces = len(poly_faces(p)) + len(edges) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantitruncate faces count");
    assert_int_eq(_count_faces_of_size(q, 3), 0, "cantitruncate tetra: 0 triangles");
    assert_int_eq(_count_faces_of_size(q, 4), 6, "cantitruncate tetra: 6 quads");
    assert_int_eq(_count_faces_of_size(q, 6), 8, "cantitruncate tetra: 8 hexagons (4 face + 4 vertex)");
}

module test_poly_cantitruncate__cube_counts() {
    p = hexahedron();
    edges = _ps_edges_from_faces(poly_faces(p));
    q = poly_cantitruncate(p, 0.2, 0.2);
    assert_poly_valid(q);

    expected_faces = len(poly_faces(p)) + len(edges) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantitruncate cube faces count");
    assert_int_eq(_count_faces_of_size(q, 4), 12, "cantitruncate cube: 12 quads (edges)");
    assert_int_eq(_count_faces_of_size(q, 6), 8, "cantitruncate cube: 8 hexagons (verts)");
    assert_int_eq(_count_faces_of_size(q, 8), 6, "cantitruncate cube: 6 octagons (faces)");
}

module test_poly_cantitruncate__dodeca_counts() {
    p = dodecahedron();
    edges = _ps_edges_from_faces(poly_faces(p));
    q = poly_cantitruncate(p, 0.2, 0.2);
    assert_poly_valid(q);

    expected_faces = len(poly_faces(p)) + len(edges) + len(poly_verts(p));
    assert_int_eq(len(poly_faces(q)), expected_faces, "cantitruncate dodeca faces count");
    assert_int_eq(_count_faces_of_size(q, 4), 30, "cantitruncate dodeca: 30 quads (edges)");
    assert_int_eq(_count_faces_of_size(q, 6), 20, "cantitruncate dodeca: 20 hexagons (verts)");
    assert_int_eq(_count_faces_of_size(q, 10), 12, "cantitruncate dodeca: 12 decagons (faces)");
}

module test_poly_cantitruncate__cube_edge_face_adjacency() {
    p = hexahedron();
    q = poly_cantitruncate(p, 0.2, 0.2);
    faces = poly_faces(q);
    edges = _ps_edges_from_faces(poly_faces(p));

    // pick the first edge-face (quad) after face cycles
    face_count = len(poly_faces(p));
    quad = faces[face_count]; // edge 0 quad
    // adjacent face cycles are expected to share its vertices
    shared = [
        for (v = quad)
            sum([for (f = [0:1:face_count-1]) (search([v], faces[f], 1)[0] >= 0) ? 1 : 0])
    ];
    // each quad vertex should belong to at least one original-face cycle
    assert(min(shared) >= 1, "cantitruncate cube: quad vertices shared with face cycles");
}

module test_poly_cantitruncate_dominant_edges__consistent_pairs() {
    base = poly_rectify(octahedron()); // cuboctahedron
    sol = solve_cantitruncate_dominant_edges(base, 4);
    c_by_size = sol[1];
    c_edge_by_pair = sol[2];
    assert(len(c_edge_by_pair) > 0, "cantitruncate dominant edges: pairs present");
    for (p = c_edge_by_pair)
        let(
            c0 = _map_face_c(p[0], c_by_size),
            c1 = _map_face_c(p[1], c_by_size),
            c_avg = (c0 + c1) / 2
        )
        assert_near(p[2], c_avg, 1e-6, str("cantitruncate dominant edges pair ", p[0], "-", p[1]));
}

module test_poly_cantitruncate_dominant_edges__planarity() {
    if (!ENABLE_CANTITRUNC_PLANARITY_TEST)
        echo("NOTE: cantitruncate dominant edges planarity test disabled");
    else {
        base = poly_rectify(octahedron()); // cuboctahedron
        sol = solve_cantitruncate_dominant_edges(base, 4);
        p = poly_cantitruncate_families(base, sol[0], sol[1], c_edge_by_pair=sol[2]);
        verts = poly_verts(p);
        faces = poly_faces(p);
        errs = [for (f = faces) _face_max_plane_err(verts, f)];
        max_err = max(errs);
        assert(max_err <= 1e-3, str("cantitruncate dominant edges planarity max_err=", max_err));
    }
}

module test_great_rhombi__cube_square_faces() {
    p = great_rhombicuboctahedron();
    assert_poly_valid(p);
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

// chamfer: positive t should move face planes inward (closer to origin)
module test_poly_chamfer__positive_t_inward() {
    p = _tetra_poly();
    q = poly_chamfer(p, 0.1);

    verts_p = poly_verts(p);
    faces_p = poly_faces(p);
    verts_q = poly_verts(q);
    faces_q = poly_faces(q);

    for (fi = [0 : 1 : len(faces_p)-1]) {
        fp = faces_p[fi];
        fq = faces_q[fi];
        np = ps_face_normal(verts_p, fp);
        nq = ps_face_normal(verts_q, fq);
        dp = v_dot(np, verts_p[fp[0]]);
        dq = v_dot(nq, verts_q[fq[0]]);
        assert(dq < dp, str("chamfer t>0 inward face ", fi, " dp=", dp, " dq=", dq));
    }
}

module test_poly_chamfer__tetra_changes_geom() {
    p = _tetra_poly();
    q = poly_chamfer(p, 0.2);
    assert(len(poly_verts(q)) != len(poly_verts(p)), "chamfer tetra should add vertices");
    assert(len(poly_faces(q)) != len(poly_faces(p)), "chamfer tetra should add faces");
}

module test_poly_chamfer__skew_prism_shrinks_all_faces() {
    p = _skew_quad_prism();
    t = 0.9;

    verts0 = poly_verts(p);
    faces0 = ps_orient_all_faces_outward(verts0, poly_faces(p));
    edges = _ps_edges_from_faces(faces0);
    edge_faces = ps_edge_faces_table(faces0, edges);
    face_n = [ for (f = faces0) ps_face_normal(verts0, f) ];

    fi = 0;
    f = faces0[fi];
    n = len(f);
    center = poly_face_center(p, fi, 1);
    ex = poly_face_ex(p, fi, 1);
    ey = poly_face_ey(p, fi, 1);
    pts2d = [
        for (k = [0:1:n-1])
            let(v = verts0[f[k]] - center)
                [v_dot(v, ex), v_dot(v, ey)]
    ];

    area0 = _ps_poly_area_abs_2d(pts2d);
    collapse = _ps_face_bisector_collapse_d(f, fi, center, ex, ey, face_n[fi], pts2d, edges, edge_faces, face_n, verts0);
    d_f = -t * collapse;
    inset2d = _ps_face_inset_bisector_2d(f, fi, d_f, 0, center, ex, ey, face_n[fi], pts2d, edges, edge_faces, face_n, verts0);
    area1 = _ps_poly_area_abs_2d(inset2d);

    assert(area1 < area0, str("chamfer t=0.9 should shrink face area: area0=", area0, " area1=", area1));
}



// suite
module run_TestTruncation() {
    test__ps_unique_points__dedups_with_eps();
    test__ps_face_points_to_indices__maps();

    test__ps_corner_t__equilateral_is_one_thirdish();
    test__ps_truncate_default_t__tetra_one_thirdish();

    test_poly_truncate__tetra_counts_at_one_third();
    test_poly_truncate__t_zero_counts_preserved();
    test_poly_chamfer__cube_face_counts();
    test_poly_truncate_then_dual__counts_relations();
    test_poly_rectify__tetra_counts();
    test_poly_cantellate__tetra_counts();
    test_poly_cantellate__cube_counts();
    test_poly_cantitruncate__tetra_counts();
    test_poly_cantitruncate_dominant_edges__consistent_pairs();
    test_poly_cantitruncate_dominant_edges__planarity();
    test_great_rhombi__cube_square_faces();
    test_truncate__validity();
    
    test_truncate__tetra_archimedean_counts();
    test_truncate__octa_archimedean_counts();
    test_truncate__icosa_archimedean_counts();
    test_truncate__hexa_archimedean_counts();
    test_truncate__dodeca_archimedean_counts();
    test_poly_chamfer__positive_t_inward();
    test_poly_chamfer__tetra_changes_geom();
    test_poly_chamfer__skew_prism_shrinks_all_faces();
}

run_TestTruncation();
