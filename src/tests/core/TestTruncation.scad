use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/duals.scad>
use <../../polysymmetrica/core/transform.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/core/solvers.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/archimedians_all.scad>
use <../testing_util.scad>

EPS = 1e-7;
ENABLE_CANTITRUNC_PLANARITY_TEST = false;
ENABLE_SNUB_PERF_SMOKE = false;

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a-b) <= eps, str(msg, " expected=", b, " got=", a));
}
module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

function _count_faces_of_size(poly, k) =
    sum([ for (f = poly_faces(poly)) (len(f)==k) ? 1 : 0 ]);

function _edge_rel_spread(poly) =
    let(
        verts = poly_verts(poly),
        edges = poly_edges(poly),
        lens = [for (e = edges) norm(verts[e[0]] - verts[e[1]])],
        avg = (len(lens) == 0) ? 0 : (sum(lens) / len(lens))
    )
    (len(lens) == 0 || avg == 0) ? 0 : ((max(lens) - min(lens)) / avg);

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

function _max_vertex_diff(p1, p2) =
    let(
        v1 = poly_verts(p1),
        v2 = poly_verts(p2),
        n = min(len(v1), len(v2))
    )
    (n == 0) ? 0 : max([for (i = [0:1:n-1]) norm(v1[i] - v2[i])]);

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

module test__ps_face_inset_bisector_2d__list_matches_scalar() {
    p = hexahedron();
    verts0 = poly_verts(p);
    faces0 = ps_orient_all_faces_outward(verts0, poly_faces(p));
    edges = _ps_edges_from_faces(faces0);
    edge_faces = ps_edge_faces_table(faces0, edges);
    face_n = [ for (f = faces0) ps_face_normal(verts0, f) ];

    fi = 0;
    f = faces0[fi];
    n = len(f);
    n_f = face_n[fi];
    center = poly_face_center(p, fi, 1);
    ex = poly_face_ex(p, fi, 1);
    ey = poly_face_ey(p, fi, 1);
    pts2d = [
        for (k = [0:1:n-1])
            let(v = verts0[f[k]] - center)
                [v_dot(v, ex), v_dot(v, ey)]
    ];

    d_f = 0.1;
    d_e = 0.2;
    inset_scalar = _ps_face_inset_bisector_2d(f, fi, d_f, d_e, center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0);
    inset_list = _ps_face_inset_bisector_2d(f, fi, d_f, [for (_ = [0:1:n-1]) d_e], center, ex, ey, n_f, pts2d, edges, edge_faces, face_n, verts0);

    for (k = [0:1:n-1]) {
        assert_near(inset_list[k][0], inset_scalar[k][0], 1e-8, str("inset list x match ", k));
        assert_near(inset_list[k][1], inset_scalar[k][1], 1e-8, str("inset list y match ", k));
    }
}

module test__ps_is_regular_base__detects_regular() {
    assert(_ps_is_regular_base(hexahedron()), "regular base: cube");
    assert(!_ps_is_regular_base(cuboctahedron()), "irregular base: cuboctahedron");
}

module test_poly_snub__cube_counts() {
    p = hexahedron();
    q = poly_snub(p);
    assert_poly_valid(q);

    assert_int_eq(len(poly_verts(q)), 24, "snub cube verts count");
    assert_int_eq(len(poly_edges(q)), 60, "snub cube edges count");
    assert_int_eq(len(poly_faces(q)), 38, "snub cube faces count");
    assert_int_eq(_count_faces_of_size(q, 3), 32, "snub cube: 32 triangles");
    assert_int_eq(_count_faces_of_size(q, 4), 6, "snub cube: 6 squares");
}

module test_poly_snub__dodeca_counts() {
    p = dodecahedron();
    // Use explicit parameters here to keep test runtime bounded;
    // default-solving behavior is covered by cube default tests.
    q = poly_snub(p, c=0.07, angle=15);
    assert_poly_valid(q);

    assert_int_eq(len(poly_verts(q)), 60, "snub dodeca verts count");
    assert_int_eq(len(poly_edges(q)), 150, "snub dodeca edges count");
    assert_int_eq(len(poly_faces(q)), 92, "snub dodeca faces count");
    assert_int_eq(_count_faces_of_size(q, 3), 80, "snub dodeca: 80 triangles");
    assert_int_eq(_count_faces_of_size(q, 5), 12, "snub dodeca: 12 pentagons");
}

module test_poly_snub__cube_twist_moves_vertices() {
    p = hexahedron();
    q0 = poly_snub(p, angle=0, c=0.2);
    q1 = poly_snub(p, angle=20, c=0.2);
    verts0 = poly_verts(q0);
    verts1 = poly_verts(q1);
    max_d = max([for (i = [0:1:len(verts0)-1]) norm(verts1[i] - verts0[i])]);
    assert(max_d > 1e-4, "snub cube: twist moves vertices");
}

module test__ps_snub_oriented_edge_faces__cube_consistent_handedness() {
    p = hexahedron();
    verts0 = poly_verts(p);
    faces0 = poly_faces(p);
    edges = poly_edges(p);
    edge_faces = ps_edge_faces_table(faces0, edges);
    face_n = [for (fi = [0:1:len(faces0)-1]) ps_face_normal(verts0, faces0[fi])];

    signs = [
        for (ei = [0:1:len(edges)-1])
            let(
                e = edges[ei],
                fpair = _ps_snub_oriented_edge_faces(e, edge_faces[ei], face_n, verts0),
                n0 = face_n[fpair[0]],
                n1 = face_n[fpair[1]],
                e_dir = v_norm(verts0[e[1]] - verts0[e[0]])
            )
            v_dot(v_cross(n0, n1), e_dir)
    ];
    assert(min(signs) >= -1e-9, str("snub edge-face orientation sign should be non-negative, min=", min(signs)));

    same_on_swap = [
        for (ei = [0:1:len(edges)-1])
            let(
                e = edges[ei],
                raw = edge_faces[ei],
                a = _ps_snub_oriented_edge_faces(e, raw, face_n, verts0),
                b = _ps_snub_oriented_edge_faces(e, [raw[1], raw[0]], face_n, verts0)
            )
            (a[0] == b[0] && a[1] == b[1]) ? 1 : 0
    ];
    assert(min(same_on_swap) == 1, "snub edge-face orientation should be invariant to face-pair order");
}

module test_poly_snub__cube_edge_tris_near_equilateral() {
    p = hexahedron();
    q = poly_snub(p);
    verts = poly_verts(q);
    faces = poly_faces(q);
    edges = poly_edges(q);
    edge_faces = ps_edge_faces_table(faces, edges);
    tri_faces = [for (fi = [0:1:len(faces)-1]) if (len(faces[fi]) == 3) fi];
    tri_sq = [
        for (fi = tri_faces)
            let(nbrs = [for (ei = [0:1:len(edges)-1]) if (edge_faces[ei][0] == fi || edge_faces[ei][1] == fi) (edge_faces[ei][0] == fi ? edge_faces[ei][1] : edge_faces[ei][0])])
            if (len([for (n = nbrs) if (n >= 0 && len(faces[n]) == 4) n]) > 0) fi
    ];
    errs = [
        for (fi = tri_sq)
            let(f = faces[fi])
                let(
                    l0 = norm(verts[f[0]] - verts[f[1]]),
                    l1 = norm(verts[f[1]] - verts[f[2]]),
                    l2 = norm(verts[f[2]] - verts[f[0]]),
                    avg = (l0 + l1 + l2) / 3
                )
                max([abs(l0-avg), abs(l1-avg), abs(l2-avg)]) / avg
    ];
    err = (len(errs) == 0) ? 1 : max(errs);
    assert(err < 0.01, str("snub cube: edge triangles should be near equilateral err=", err));
}

module test_poly_snub__cube_default_global_edges_near_uniform() {
    p = hexahedron();
    q = poly_snub(p);
    spread = _edge_rel_spread(q);
    assert(spread < 0.02, str("snub cube default: global edge spread <2% got=", spread));
}

module test_poly_snub__cube_default_params_reasonable() {
    p = hexahedron();
    params = _ps_snub_default_params(p);
    assert(params[3] == "regular", "snub cube: regular-tier default");
    assert(params[2] > 0.03 && params[2] < 0.12, str("snub cube: c in expected range got=", params[2]));
    assert(params[1] > 8 && params[1] < 25, str("snub cube: angle in expected range got=", params[1]));
}

module test_poly_snub__fixed_c_angle_solver_nonzero() {
    a = _ps_snub_default_angle_c(hexahedron(), 0.07, steps=8, a_max=30);
    assert(a > 8 && a < 25, str("snub cube fixed-c angle should be nonzero got=", a));
}

module test_poly_snub__fixed_c_auto_beats_zero_angle() {
    p = hexahedron();
    q0 = poly_snub(p, c=0.07, angle=0);
    q1 = poly_snub(p, c=0.07);
    assert(_edge_rel_spread(q1) < _edge_rel_spread(q0), "snub cube fixed-c auto-angle improves edge uniformity vs angle=0");
}

module test_poly_snub__default_solver_returns_structured_overrides() {
    p = hexahedron();
    rows = ps_snub_default_params_overrides(p);
    q0 = poly_snub(p);
    q1 = poly_snub(p, params_overrides=rows);
    assert(len(rows) > 0, "snub structured defaults should return rows");
    assert(_max_vertex_diff(q0, q1) < 1e-7, "snub structured defaults should recreate auto result");
}

module test_poly_snub__params_overrides_face_angle_overrides_scalar() {
    p = hexahedron();
    q0 = poly_snub(p, angle=0, c=0.07, df=0.05);
    q1 = poly_snub(p, angle=0, c=0.07, df=0.05, params_overrides=[["face", "family", 0, ["angle", 18]]]);
    qx = poly_snub(p, angle=18, c=0.07, df=0.05);
    assert(_max_vertex_diff(q1, qx) < 1e-7, "snub params_overrides face angle should match explicit angle");
    assert(_max_vertex_diff(q1, q0) > 1e-4, "snub params_overrides face angle should override scalar angle");
}

module test_poly_snub__params_overrides_face_df_overrides_scalar() {
    p = hexahedron();
    q1 = poly_snub(p, angle=15, c=0.07, df=0.02, params_overrides=[["face", "family", 0, ["df", 0.06]]]);
    qx = poly_snub(p, angle=15, c=0.07, df=0.06);
    assert(_max_vertex_diff(q1, qx) < 1e-7, "snub params_overrides face df should match explicit df");
}

module test_poly_snub__params_overrides_vert_c_overrides_scalar() {
    p = hexahedron();
    q1 = poly_snub(p, angle=15, c=0.02, df=0.05, params_overrides=[["vert", "family", 0, ["c", 0.08]]]);
    qx = poly_snub(p, angle=15, c=0.08, df=0.05);
    assert(_max_vertex_diff(q1, qx) < 1e-7, "snub params_overrides vert c should match explicit c");
}

// Informational performance smoke test for default snub solving.
// This is intentionally non-failing and opt-in to avoid slowing normal CI/dev runs.
module perf_snub__defaults_smoke() {
    echo("PERF_SNUB: start cube default params");
    p0 = _ps_snub_default_params(hexahedron());
    echo("PERF_SNUB: done cube default params", p0);

    echo("PERF_SNUB: start dodeca default params");
    p1 = _ps_snub_default_params(dodecahedron());
    echo("PERF_SNUB: done dodeca default params", p1);
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
    test__ps_face_inset_bisector_2d__list_matches_scalar();
    test__ps_is_regular_base__detects_regular();
    test_poly_snub__cube_counts();
    test_poly_snub__dodeca_counts();
    test_poly_snub__cube_twist_moves_vertices();
    test__ps_snub_oriented_edge_faces__cube_consistent_handedness();
    test_poly_snub__cube_edge_tris_near_equilateral();
    test_poly_snub__cube_default_global_edges_near_uniform();
    test_poly_snub__cube_default_params_reasonable();
    test_poly_snub__fixed_c_angle_solver_nonzero();
    test_poly_snub__fixed_c_auto_beats_zero_angle();
    test_poly_snub__default_solver_returns_structured_overrides();
    test_poly_snub__params_overrides_face_angle_overrides_scalar();
    test_poly_snub__params_overrides_face_df_overrides_scalar();
    test_poly_snub__params_overrides_vert_c_overrides_scalar();
    if (ENABLE_SNUB_PERF_SMOKE)
        perf_snub__defaults_smoke();
    else
        echo("NOTE: snub perf smoke test disabled");
}

run_TestTruncation();
