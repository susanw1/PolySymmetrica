use <../../polysymmetrica/core/placement.scad>
use <../../polysymmetrica/core/proxy_interaction.scad>
use <../../polysymmetrica/core/face_regions.scad>
use <../../polysymmetrica/core/classify.scad>
use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/archimedians_all.scad>

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_vec3_near(v, w, eps=1e-9, msg="") {
    assert(norm(v - w) <= eps, str(msg, " expected=", w, " got=", v));
}

function _test_face_cut_ctx(poly, fi) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        face = faces[fi],
        center = poly_face_center(poly, fi, 1),
        ex = poly_face_ex(poly, fi, 1),
        ey = poly_face_ey(poly, fi, 1),
        ez = poly_face_ez(poly, fi, 1),
        face_verts_local = [
            for (vid = face)
                let(p = verts[vid] - center)
                    [v_dot(p, ex), v_dot(p, ey), v_dot(p, ez)]
        ],
        poly_verts_local_raw = [
            for (vi = [0:1:len(verts)-1])
                let(p = verts[vi] - center)
                    [v_dot(p, ex), v_dot(p, ey), v_dot(p, ez)]
        ],
        zvals = [for (p = face_verts_local) p[2]],
        zmean = (len(zvals) == 0) ? 0 : sum(zvals) / len(zvals),
        face_pts2d = [for (p = face_verts_local) [p[0], p[1]]],
        poly_verts_local = [for (p = poly_verts_local_raw) [p[0], p[1], p[2] - zmean]]
    )
    [face_pts2d, poly_verts_local, center, ex, ey];

module test_place_on_faces__family_ids_and_counts_from_classify() {
    p = rhombicuboctahedron();
    faces = poly_faces(p);
    cls = poly_classify(p, 1, 1e-6, 1, false);
    ids = ps_classify_face_ids(cls, len(faces));
    counts = ps_classify_counts(cls);

    place_on_faces(p, classify = cls) {
        assert_int_eq($ps_face_family_id, ids[$ps_face_idx], "face family id");
        assert_int_eq($ps_face_family_count, counts[0], "face family count");
        assert_int_eq($ps_edge_family_count, counts[1], "edge family count from face placement");
        assert_int_eq($ps_vertex_family_count, counts[2], "vertex family count from face placement");
    }
}

module test_place_on_edges__family_ids_and_counts_from_classify() {
    p = rhombicuboctahedron();
    faces = poly_faces(p);
    edges = _ps_edges_from_faces(faces);
    cls = poly_classify(p, 1, 1e-6, 1, false);
    ids = ps_classify_edge_ids(cls, len(edges));
    counts = ps_classify_counts(cls);

    place_on_edges(p, classify = cls) {
        assert_int_eq($ps_edge_family_id, ids[$ps_edge_idx], "edge family id");
        assert_int_eq($ps_face_family_count, counts[0], "face family count from edge placement");
        assert_int_eq($ps_edge_family_count, counts[1], "edge family count");
        assert_int_eq($ps_vertex_family_count, counts[2], "vertex family count from edge placement");
    }
}

module test_place_on_vertices__family_ids_and_counts_from_classify() {
    p = rhombicuboctahedron();
    verts = poly_verts(p);
    cls = poly_classify(p, 1, 1e-6, 1, false);
    ids = ps_classify_vert_ids(cls, len(verts));
    counts = ps_classify_counts(cls);

    place_on_vertices(p, classify = cls) {
        assert_int_eq($ps_vertex_family_id, ids[$ps_vertex_idx], "vertex family id");
        assert_int_eq($ps_face_family_count, counts[0], "face family count from vertex placement");
        assert_int_eq($ps_edge_family_count, counts[1], "edge family count from vertex placement");
        assert_int_eq($ps_vertex_family_count, counts[2], "vertex family count");
    }
}

module test_place_on_faces__auto_classify_matches_precomputed() {
    p = rhombicuboctahedron();
    faces = poly_faces(p);
    opts = [1, 1e-6, 1, false];
    cls = poly_classify(p, opts[0], opts[1], opts[2], opts[3]);
    ids = ps_classify_face_ids(cls, len(faces));

    place_on_faces(p, classify_opts = opts) {
        assert_int_eq($ps_face_family_id, ids[$ps_face_idx], "auto classify face family id");
    }
}

module test_place_on_faces__indices_filters_exact_face() {
    p = hexahedron();
    place_on_faces(p, indices = [2]) {
        assert_int_eq($ps_face_idx, 2, "face index filter should visit only requested face");
    }
}

module test_place_on_all__cube_single_family() {
    p = hexahedron();
    cls = poly_classify(p, 1, 1e-6, 1, false);
    counts = ps_classify_counts(cls);
    assert_int_eq(counts[0], 1, "cube face families");
    assert_int_eq(counts[1], 1, "cube edge families");
    assert_int_eq(counts[2], 1, "cube vertex families");

    place_on_faces(p, classify = cls) {
        assert_int_eq($ps_face_family_id, 0, "cube face family id");
        assert_int_eq($ps_face_family_count, 1, "cube face family count");
    }
    place_on_edges(p, classify = cls) {
        assert_int_eq($ps_edge_family_id, 0, "cube edge family id");
        assert_int_eq($ps_edge_family_count, 1, "cube edge family count");
    }
    place_on_vertices(p, classify = cls) {
        assert_int_eq($ps_vertex_family_id, 0, "cube vertex family id");
        assert_int_eq($ps_vertex_family_count, 1, "cube vertex family count");
    }
}

module test_place_on_edges__no_auto_classify_by_default() {
    p = hexahedron();
    place_on_edges(p) {
        assert(is_undef($ps_edge_family_id), "default placement: edge family id should be undef without classify");
        assert(is_undef($ps_face_family_count), "default placement: face family count should be undef without classify");
        assert(is_undef($ps_edge_family_count), "default placement: edge family count should be undef without classify");
        assert(is_undef($ps_vertex_family_count), "default placement: vertex family count should be undef without classify");
    }
}

module test_place_on_edges__indices_filters_exact_edge() {
    p = hexahedron();
    place_on_edges(p, indices = [3]) {
        assert_int_eq($ps_edge_idx, 3, "edge index filter should visit only requested edge");
    }
}

module test_place_on_edges__ez_world_matches_adjacent_face_bisector() {
    verts = [
        [-1.0, -1.0, -1.0],
        [ 1.0, -1.0, -1.0],
        [ 0.0,  1.0, -1.0],
        [-0.1, -0.6,  1.0],
        [ 2.1, -1.1,  1.0],
        [ 0.8,  0.7,  1.0]
    ];
    faces = [
        [0, 1, 2],
        [5, 4, 3],
        [0, 3, 4, 1],
        [1, 4, 5, 2],
        [2, 5, 3, 0]
    ];
    p = make_poly(verts, faces, 1);
    verts0 = poly_verts(p);
    faces0 = ps_orient_all_faces_outward(verts0, poly_faces(p));
    edges0 = _ps_edges_from_faces(faces0);
    edge_faces0 = ps_edge_faces_table(faces0, edges0);
    face_n0 = [for (f = faces0) ps_face_normal(verts0, f)];
    ei = ps_find_edge_index(edges0, 0, 1);
    adj = edge_faces0[ei];
    radial = v_norm((verts0[0] + verts0[1]) / 2);
    ex = v_norm(verts0[1] - verts0[0]);
    bis0 = face_n0[adj[0]] + face_n0[adj[1]];
    bis1 = (v_dot(bis0, radial) < 0) ? -bis0 : bis0;
    bis_proj = bis1 - ex * v_dot(bis1, ex);
    expected_ez = v_norm(bis_proj);
    dotn = v_dot(face_n0[adj[0]], face_n0[adj[1]]);
    expected_dihedral = 180 - acos(ps_clamp(dotn, -1, 1));

    place_on_edges(p, indices = [ei]) {
        assert_vec3_near($ps_edge_ez_world, expected_ez, 1e-8, "edge ez should bisect adjacent face normals");
        assert(abs($ps_dihedral - expected_dihedral) <= 1e-8, str("edge dihedral wrong expected=", expected_dihedral, " got=", $ps_dihedral));
    }
}

module test_place_on_vertices__indices_filters_exact_vertex() {
    p = hexahedron();
    place_on_vertices(p, indices = [5]) {
        assert_int_eq($ps_vertex_idx, 5, "vertex index filter should visit only requested vertex");
    }
}

module test_place_on_face_segments__star_face_split() {
    p = poly_antiprism(5, 2);
    place_on_faces(p) {
        if ($ps_face_idx == 0) {
            place_on_face_segments(mode="evenodd") {
                assert(!is_undef($ps_seg_idx), "segment idx should be defined");
                assert(!is_undef($ps_seg_count), "segment count should be defined");
                assert($ps_seg_count > 1, "star face should split into multiple segments");
                assert($ps_seg_vertex_count >= 3, "segment should have at least 3 vertices");
                assert(len($ps_seg_pts2d) == $ps_seg_vertex_count, "segment pts2d count");
                assert(len($ps_seg_pts3d_local) == $ps_seg_vertex_count, "segment pts3d count");
                assert(len($ps_seg_parent_face_edge_idx) == $ps_seg_vertex_count, "segment parent-edge mapping");
                assert($ps_face_has_segments == ($ps_seg_count > 1), "segment split flag consistency");
            }
        }
    }
}

module test_place_on_face_segments__default_nonzero_keeps_filled_star() {
    p = poly_antiprism(5, 2);
    place_on_faces(p) {
        if ($ps_face_idx == 0) {
            place_on_face_segments() {
                assert_int_eq($ps_seg_count, 1, "default nonzero face segments should keep filled star as one region");
            }
        }
    }
}

module test_place_on_face_filled_boundary_segments__triangle_tracks_face_edges() {
    p = tetrahedron();
    place_on_faces(p, indices = [0]) {
        assert_int_eq(len(ps_face_filled_boundary_segments($ps_face_pts3d_local, 1e-9)), 3,
            "triangle face should report three filled boundary segments");
        place_on_face_filled_boundary_segments() {
            assert_int_eq($ps_face_boundary_seg_count, 3, "triangle should have three true boundary segments");
            assert(!is_undef($ps_face_boundary_seg_source_edge_idx), "triangle boundary segment should keep source edge idx");
            assert(!is_undef($ps_face_boundary_seg_dihedral), "triangle boundary segment should expose source dihedral");
            assert($ps_face_boundary_seg_len > 0, "triangle boundary segment should have positive length");
            assert(abs(norm($ps_face_boundary_seg_left_normal2d) - 1) < 1e-8, "triangle boundary left normal should be unit");
            assert(abs(norm($ps_face_boundary_seg_inward2d) - 1) < 1e-8, "triangle boundary inward normal should be unit");
            assert(abs(abs(v_dot($ps_face_boundary_seg_left_normal2d, $ps_face_boundary_seg_inward2d)) - 1) < 1e-8,
                "triangle inward direction should align with +/- left normal");
            assert(abs(v_dot($ps_face_boundary_seg2d[1] - $ps_face_boundary_seg2d[0], $ps_face_boundary_seg_inward2d)) < 1e-8,
                "triangle inward should be perpendicular to segment");
        }
    }
}

module test_place_on_face_filled_boundary_segments__pentagram_returns_decagon_segments() {
    p = poly_antiprism(5, 2);
    place_on_faces(p) {
        if ($ps_face_idx == 0) {
            segs = ps_face_filled_boundary_segments($ps_face_pts3d_local, 1e-9);
            src_counts = [for (ei = [0:1:4]) sum([for (s = segs) (s[1] == ei) ? 1 : 0])];

            assert_int_eq(len(segs), 10, "star face should report ten true boundary subsegments");
            for (cnt = src_counts)
                assert_int_eq(cnt, 2, "each star edge should contribute two true boundary subsegments");

            place_on_face_filled_boundary_segments() {
                assert_int_eq($ps_face_boundary_seg_count, 10, "star face should expose ten true boundary subsegments");
                assert($ps_face_boundary_seg_len > 0, "star boundary subsegment should have positive length");
                assert(!is_undef($ps_face_boundary_seg_source_edge_idx), "star boundary subsegment should keep source edge idx");
                assert($ps_face_boundary_seg_source_t0 >= -1e-8 && $ps_face_boundary_seg_source_t0 <= 1 + 1e-8,
                    str("star t0 out of range ", $ps_face_boundary_seg_source_t0));
                assert($ps_face_boundary_seg_source_t1 >= -1e-8 && $ps_face_boundary_seg_source_t1 <= 1 + 1e-8,
                    str("star t1 out of range ", $ps_face_boundary_seg_source_t1));
                assert(abs($ps_face_boundary_seg_source_t1 - $ps_face_boundary_seg_source_t0) > 1e-8,
                    "star boundary subsegment should have non-zero param span");
            }
        }
    }
}

module test_place_on_face_filled_boundary_edges__triangle_maps_to_edge_frame() {
    p = tetrahedron();
    place_on_faces(p, indices = [0]) {
        place_on_face_filled_boundary_edges() {
            assert_int_eq(len($ps_edge_adj_faces_idx), 2, "triangle boundary edge should carry adjacent face pair");
            assert(!is_undef($ps_dihedral), "triangle boundary edge should expose dihedral");
            assert(abs($ps_edge_len - $ps_face_boundary_seg_len) < 1e-8, "boundary edge frame should use subsegment length");
            assert(abs(norm($ps_edge_ex_world) - 1) < 1e-8, "boundary edge ex should be unit");
            assert(abs(norm($ps_edge_ey_world) - 1) < 1e-8, "boundary edge ey should be unit");
            assert(abs(norm($ps_edge_ez_world) - 1) < 1e-8, "boundary edge ez should be unit");
            assert(abs(v_dot($ps_edge_ex_world, $ps_edge_ey_world)) < 1e-8, "boundary edge ex/ey should be orthogonal");
            assert(abs(v_dot($ps_edge_ex_world, $ps_edge_ez_world)) < 1e-8, "boundary edge ex/ez should be orthogonal");
            assert(abs(v_dot($ps_edge_ey_world, $ps_edge_ez_world)) < 1e-8, "boundary edge ey/ez should be orthogonal");
            assert(abs(v_dot($ps_edge_ez_world, [0, 0, 1])) < 1, "boundary edge ez should tilt away from face z when dihedral is finite");
        }
    }
}

module test_place_on_face_filled_boundary_edges__pentagram_uses_true_subsegments() {
    p = poly_antiprism(5, 2);
    place_on_faces(p) {
        if ($ps_face_idx == 0) {
            segs = ps_face_filled_boundary_segments($ps_face_pts3d_local, 1e-9);
            assert_int_eq(len(segs), 10, "star boundary edge frame baseline should have ten subsegments");

            place_on_face_filled_boundary_edges() {
                assert_int_eq($ps_face_boundary_seg_count, 10, "star boundary edge frame should iterate ten subsegments");
                assert(!is_undef($ps_edge_idx), "star boundary edge frame should carry source edge idx");
                assert(!is_undef($ps_dihedral), "star boundary edge frame should carry source dihedral");
                assert(abs($ps_edge_len - $ps_face_boundary_seg_len) < 1e-8, "star boundary edge frame should use subsegment length");
                assert(abs(norm($ps_edge_ez_world) - 1) < 1e-8, "star boundary edge ez should be unit");
                assert(abs(v_dot($ps_edge_ex_world, $ps_edge_ez_world)) < 1e-8, "star boundary edge ex/ez should be orthogonal");
            }
        }
    }
}

module test_ps_face_filled_boundary_segments__triangle_returns_three_edges() {
    pts3 = [[0, 0, 0], [4, 0, 0], [1, 3, 0]];
    segs = ps_face_filled_boundary_segments(pts3, 1e-9);
    src_counts = [for (ei = [0:1:2]) sum([for (s = segs) (s[1] == ei) ? 1 : 0])];

    assert_int_eq(len(segs), 3, "triangle should have three filled boundary segments");
    assert_int_eq(_ps_distinct_count([for (s = segs) s[1]]), 3, "triangle source edge ids should be distinct");
    for (cnt = src_counts)
        assert_int_eq(cnt, 1, "triangle should contribute one boundary segment per source edge");
    for (s = segs) {
        assert(!is_undef(s[1]), "triangle boundary segment should have source edge id");
        assert(!is_undef(s[2]) && !is_undef(s[3]), "triangle boundary segment should carry source params");
        assert(s[2] >= -1e-8 && s[2] <= 1 + 1e-8, str("triangle t0 out of range ", s[2]));
        assert(s[3] >= -1e-8 && s[3] <= 1 + 1e-8, str("triangle t1 out of range ", s[3]));
        assert(abs(s[3] - s[2]) > 1e-8, "triangle boundary segment should have non-zero param span");
    }
}

module test_ps_face_filled_boundary_segments__pentagram_returns_decagon_boundary() {
    pts3 = [[0,9,0], [-5,-5,0], [8,3,0], [-8,3,0], [5,-5,0]];
    cells = ps_face_filled_cells(pts3, 1e-9);
    segs = ps_face_filled_boundary_segments(pts3, 1e-9);
    src_counts = [for (ei = [0:1:4]) sum([for (s = segs) (s[1] == ei) ? 1 : 0])];
    tri_count = sum([for (c = cells) (len(c[0]) == 3) ? 1 : 0]);
    pent_count = sum([for (c = cells) (len(c[0]) == 5) ? 1 : 0]);

    assert_int_eq(len(cells), 6, "nonzero pentagram should yield five arm triangles plus one center pentagon");
    assert_int_eq(tri_count, 5, "nonzero pentagram should have five triangular arm cells");
    assert_int_eq(pent_count, 1, "nonzero pentagram should have one central pentagon cell");
    assert_int_eq(len(segs), 10, "nonzero pentagram should have ten true boundary subsegments");
    for (cnt = src_counts)
        assert_int_eq(cnt, 2, "each original pentagram edge should contribute two boundary subsegments");
    for (s = segs) {
        assert(!is_undef(s[1]), "pentagram boundary segment should have source edge id");
        assert(!is_undef(s[2]) && !is_undef(s[3]), "pentagram boundary segment should carry source params");
        assert(s[2] >= -1e-8 && s[2] <= 1 + 1e-8, str("pentagram t0 out of range ", s[2]));
        assert(s[3] >= -1e-8 && s[3] <= 1 + 1e-8, str("pentagram t1 out of range ", s[3]));
        assert(abs(s[3] - s[2]) > 1e-8, "pentagram boundary segment should have non-zero param span");
    }
}

module test_ps_face_filled_atoms__triangle_stays_single_boundary_atom() {
    pts3 = [[0, 0, 0], [4, 0, 0], [1, 3, 0]];
    atoms = ps_face_filled_atoms(pts3, 1e-9);

    assert_int_eq(len(atoms), 1, "triangle should remain one convex atom");
    assert_int_eq(len(atoms[0][0]), 3, "triangle atom should have three vertices");
    assert_int_eq(sum([for (k = atoms[0][3]) (k == "boundary") ? 1 : 0]), 3, "triangle atom edges should all be boundary");
    assert_int_eq(sum([for (k = atoms[0][3]) (k == "inner") ? 1 : 0]), 0, "triangle atom should have no inner edges");
}

module test_ps_face_filled_atoms__concave_hex_splits_and_preserves_area() {
    pts3 = [[0,0,0], [4,0,0], [4,1,0], [1,1,0], [1,4,0], [0,4,0]];
    atoms = ps_face_filled_atoms(pts3, 1e-9);
    area_poly = abs(_ps_seg_poly_area2([for (p = pts3) [p[0], p[1]]]));
    area_atoms = _list_sum([for (a = atoms) abs(_ps_seg_poly_area2(a[0]))]);
    inner_edges = sum([for (a = atoms) sum([for (k = a[3]) (k == "inner") ? 1 : 0])]);

    assert(len(atoms) > 1, "concave hex should split into multiple convex atoms");
    assert(abs(area_atoms - area_poly) < 1e-6, str("concave atom area mismatch poly=", area_poly, " atoms=", area_atoms));
    assert(inner_edges > 0, "concave atomization should introduce inner edges");
}

module test_place_on_faces__local_z_origin_consistent_for_face_and_poly_verts() {
    // Build a slightly warped cube so at least some faces are non-planar.
    base = hexahedron();
    verts0 = poly_verts(base);
    faces0 = poly_faces(base);
    verts = [
        for (i = [0:1:len(verts0)-1])
            (i == 0) ? (verts0[i] + [0.35, -0.2, 0.45]) : verts0[i]
    ];
    p = make_poly(verts, faces0);

    place_on_faces(p)
        let(
            face = faces0[$ps_face_idx],
            zmean = sum([for (q = $ps_face_pts3d_local) q[2]]) / len($ps_face_pts3d_local)
        ) {
            assert(abs(zmean) < 1e-9, str("face local z should be mean-centered, fi=", $ps_face_idx, " zmean=", zmean));
            for (k = [0:1:len(face)-1])
                let(vi = face[k])
                    assert_vec3_near(
                        $ps_poly_verts_local[vi],
                        $ps_face_pts3d_local[k],
                        1e-8,
                        str("poly/local vertex mismatch at face=", $ps_face_idx, " k=", k, " vi=", vi)
                    );
        }
}

module test_seg_cycle_probe_point__concave_inside() {
    concave = [[0,0], [4,0], [4,1], [1,1], [1,4], [0,4]];
    probe = _ps_seg_cycle_probe_point(concave, 1e-9);
    assert(_ps_seg_point_in_poly_evenodd(probe, concave, 1e-9), str("probe should be inside concave polygon, probe=", probe));
    assert(!_ps_seg_point_on_poly_boundary(probe, concave, 1e-8), str("probe should not lie on boundary, probe=", probe));
}

function _tri2_area(a, b, c) =
    abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2;

function _list_sum(xs, i=0, acc=0) =
    (i >= len(xs)) ? acc : _list_sum(xs, i + 1, acc + xs[i]);

module test_seg_face_tris3__concave_area_preserved() {
    // Concave simple polygon: fan triangulation would over-cover; ear clipping should preserve area.
    pts3 = [[0,0,0], [4,0,0], [4,1,0], [1,1,0], [1,4,0], [0,4,0]];
    tris = _ps_seg_face_tris3([0,1,2,3,4,5], pts3, 1e-9);
    area_poly = abs(_ps_seg_poly_area2([for (p = pts3) [p[0], p[1]]]));
    area_tris = sum([
        for (t = tris)
            _tri2_area(
                [t[0][0], t[0][1]],
                [t[1][0], t[1][1]],
                [t[2][0], t[2][1]]
            )
    ]);
    assert_int_eq(len(tris), 4, "concave hex should triangulate to 4 triangles");
    assert(abs(area_tris - area_poly) < 1e-6, str("concave triangulation area mismatch poly=", area_poly, " tris=", area_tris));
}

module test_seg_face_tris3__star_area_matches_segments() {
    // Pentagram-style self-intersecting face loop.
    pts3 = [[0,9,0], [-5,-5,0], [8,3,0], [-8,3,0], [5,-5,0]];
    tris = _ps_seg_face_tris3([0,1,2,3,4], pts3, 1e-9, "evenodd");
    segs = ps_face_segments(pts3, "evenodd", 1e-9);
    area_tris = _list_sum([
        for (t = tris)
            _tri2_area(
                [t[0][0], t[0][1]],
                [t[1][0], t[1][1]],
                [t[2][0], t[2][1]]
            )
    ]);
    area_segs = _list_sum([for (s = segs) abs(_ps_seg_poly_area2(s[0]))]);
    assert(len(segs) > 1, "star face should split into multiple even-odd regions");
    assert(len(tris) >= 3, "triangulation should produce at least a few triangles");
    assert(abs(area_tris - area_segs) < 1e-6, str("star triangulation area mismatch tris=", area_tris, " segs=", area_segs));
}

module test_ps_face_visible_segments__cube_face_unchanged() {
    p = hexahedron();
    place_on_faces(p) {
        if ($ps_face_idx == 0) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            area_face = abs(_ps_seg_poly_area2($ps_face_pts2d));
            area_vis = _list_sum([for (s = vis) abs(_ps_seg_poly_area2(s[0]))]);
            assert_int_eq(len(vis), 1, "cube face should keep one visible cell");
            assert(abs(area_vis - area_face) < 1e-6, str("cube visible area mismatch face=", area_face, " vis=", area_vis));
            assert_int_eq(len(vis[0][2]), len(vis[0][0]), "cube visible edge-id count");
            assert_int_eq(len(vis[0][3]), len(vis[0][0]), "cube visible edge-kind count");
            assert(sum([for (k = vis[0][3]) (k == "parent") ? 1 : 0]) == len(vis[0][0]), "cube visible cell edges should all be parent");
        }
    }
}

module test_ps_face_visible_segments__star_antiprism_side_reduced() {
    p = poly_antiprism(5, 2);
    faces = poly_faces(p);
    tri_faces = [for (i = [0:1:len(faces)-1]) if (len(faces[i]) == 3) i];
    target = tri_faces[0];
    place_on_faces(p) {
        if ($ps_face_idx == target) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            area_face = abs(_ps_seg_poly_area2($ps_face_pts2d));
            area_vis = _list_sum([for (s = vis) abs(_ps_seg_poly_area2(s[0]))]);
            assert(len(vis) >= 1, "star antiprism side should keep at least one visible cell");
            assert(area_vis > 1e-6, "star antiprism visible area should stay positive");
            assert(area_vis < area_face - 1e-6, str("star antiprism side should lose hidden area face=", area_face, " vis=", area_vis));
            assert(
                sum([
                    for (s = vis)
                        sum([for (k = s[3]) (k == "cut") ? 1 : 0])
                ]) > 0,
                "star antiprism visible cells should include cut edges"
            );
        }
    }
}

module test_ps_face_visible_segments__star_prism_side_keeps_two_side_quads() {
    p = poly_prism(5, 2);
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            assert_int_eq(len(vis), 2, "star prism side should keep two visible side panels");
            for (s = vis) {
                assert_int_eq(len(s[0]), 4, "each star prism side panel should be a quad");
                assert(abs(_ps_seg_poly_area2(s[0])) > 1e-6, "star prism side panel visible area should stay positive");
                assert(sum([for (k = s[3]) (k == "cut") ? 1 : 0]) >= 1, "star prism side panels should include cut edges");
            }
        }
    }
}

module test_ps_face_visible_segments__cells_preserve_parent_winding() {
    p = poly_antiprism(5, 2);
    faces = poly_faces(p);
    tri_faces = [for (i = [0:1:len(faces)-1]) if (len(faces[i]) == 3) i];
    target = tri_faces[0];
    place_on_faces(p) {
        if ($ps_face_idx == target) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            parent_sign = (_ps_seg_poly_area2($ps_face_pts2d) >= 0) ? 1 : -1;
            for (s = vis)
                assert(
                    _ps_seg_poly_area2(s[0]) * parent_sign > 1e-9,
                    str("visible cell winding mismatch parent=", $ps_face_idx, " area=", _ps_seg_poly_area2(s[0]))
                );
        }
    }
}

module test_ps_face_visible_segments__star_prism_cut_edges_reference_cut_entries() {
    p = poly_prism(5, 2);
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            assert(len(entries) > 0, "star prism side should have cut entries");
            for (s = vis) {
                assert_int_eq(len(s[3]), len(s[4]), "star prism visible edge kind/id arity");
                for (k = [0:1:len(s[3])-1]) {
                    if (s[3][k] == "cut") {
                        cid = s[4][k];
                        assert(!is_undef(cid), "star prism cut edge should carry cut-entry id");
                        assert(cid >= 0 && cid < len(entries), str("star prism cut-entry id out of range cid=", cid, " len=", len(entries)));
                    } else {
                        assert(is_undef(s[4][k]), "star prism parent edge should not carry cut-entry id");
                    }
                }
            }
        }
    }
}

module test_ps_face_visible_segments__star_antiprism_cut_edges_reference_cut_entries() {
    p = poly_antiprism(5, 2);
    faces = poly_faces(p);
    tri_faces = [for (i = [0:1:len(faces)-1]) if (len(faces[i]) == 3) i];
    target = tri_faces[0];
    place_on_faces(p) {
        if ($ps_face_idx == target) {
            entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            assert(len(entries) > 0, "star antiprism side should have cut entries");
            assert(
                sum([
                    for (s = vis)
                        sum([for (k = [0:1:len(s[3])-1]) (s[3][k] == "cut") ? 1 : 0])
                ]) > 0,
                "star antiprism visible cells should include cut edges"
            );
            for (s = vis) {
                assert_int_eq(len(s[3]), len(s[4]), "star antiprism visible edge kind/id arity");
                for (k = [0:1:len(s[3])-1]) {
                    if (s[3][k] == "cut") {
                        cid = s[4][k];
                        assert(!is_undef(cid), "star antiprism cut edge should carry cut-entry id");
                        assert(cid >= 0 && cid < len(entries), str("star antiprism cut-entry id out of range cid=", cid, " len=", len(entries)));
                    } else {
                        assert(is_undef(s[4][k]), "star antiprism parent edge should not carry cut-entry id");
                    }
                }
            }
        }
    }
}

module test_ps_face_visible_segments__cut_pair_ids_defined_in_ctx() {
    p = poly_antiprism(n=7, p=3, angle=15);
    place_on_faces(p) {
        if ($ps_face_idx == 2 || $ps_face_idx == 5) {
            place_on_face_visible_segments() {
                for (k = [0:1:len($ps_vis_seg_edge_kinds)-1]) {
                    if ($ps_vis_seg_edge_kinds[k] == "cut") {
                        assert(
                            !is_undef($ps_vis_seg_cut_pair_ids[k]),
                            str("cut pair id should be defined face=", $ps_face_idx, " cell=", $ps_vis_seg_idx, " edge=", k)
                        );
                        assert(
                            !is_undef($ps_vis_seg_cut_run_ids[k]),
                            str("cut run id should be defined face=", $ps_face_idx, " cell=", $ps_vis_seg_idx, " edge=", k)
                        );
                    } else {
                        assert(is_undef($ps_vis_seg_cut_pair_ids[k]), "parent edge should not carry cut pair id");
                        assert(is_undef($ps_vis_seg_cut_run_ids[k]), "parent edge should not carry cut run id");
                    }
                }
            }
        }
    }
}

module test_ps_face_geom_cut_pair_ids__matching_faces_share_join_id() {
    p = poly_antiprism(n=7, p=3, angle=15);
    faces = poly_faces(p);

    ctx2 = _test_face_cut_ctx(p, 2);
    entries2 = ps_face_geom_cut_entries(ctx2[0], 2, faces, ctx2[1], 1e-8, "nonzero", true);
    pairs2 = ps_face_geom_cut_pair_ids(entries2, 2, ctx2[2], ctx2[3], ctx2[4]);

    ctx5 = _test_face_cut_ctx(p, 5);
    entries5 = ps_face_geom_cut_entries(ctx5[0], 5, faces, ctx5[1], 1e-8, "nonzero", true);
    pairs5 = ps_face_geom_cut_pair_ids(entries5, 5, ctx5[2], ctx5[3], ctx5[4]);

    common = [
        for (pid = pairs2)
            if (!is_undef(pid) && len([for (q = pairs5) if (q == pid) 1]) > 0)
                pid
    ];

    assert(len(common) > 0, "matching cut entries across faces should share at least one join id");
}

module test_ps_face_visible_cell_cut_run_ids__wraparound_same_cutter_merges() {
    kinds = ["cut", "cut", "cut"];
    cids = [4, 7, 4];
    runs = ps_face_visible_cell_cut_run_ids(kinds, cids);
    assert_int_eq(runs[0], runs[2], "wrap-around split cut span should keep same run across cycle seam");
    assert(runs[1] != runs[0], "middle distinct cutter contribution should keep distinct run id");
}

module test_ps_face_visible_segments__split_cut_spans_get_distinct_run_ids() {
    p = poly_antiprism(n=7, p=3, angle=15);
    faces = poly_faces(p);
    ctx2 = _test_face_cut_ctx(p, 2);
    entries2 = ps_face_geom_cut_entries(ctx2[0], 2, faces, ctx2[1], 1e-8, "nonzero", true);
    vis2 = ps_face_visible_segments(ctx2[0], 2, faces, ctx2[1], 1e-8, "nonzero", true);
    cell = vis2[0];
    f11_edges = [
        for (k = [0:1:len(cell[3])-1])
            if (
                cell[3][k] == "cut" &&
                ps_cut_entry_cutter_face_idx(entries2[cell[4][k]]) == 11
            )
                k
    ];

    assert_int_eq(len(f11_edges), 2, str("f2/c0 should see two split f11 cut spans edges=", f11_edges));
    assert(!is_undef(cell[5][f11_edges[0]]), "first split cut span should have run id");
    assert(!is_undef(cell[5][f11_edges[1]]), "second split cut span should have run id");
    assert(cell[5][f11_edges[0]] != cell[5][f11_edges[1]], str("split cut spans from same cutter must keep distinct run ids runs=", cell[5]));
}

module test_ps_face_geom_cut_segments__respects_fill_mode() {
    // Target square in z=0 plane plus a star-shaped cutter face tilted through the plane.
    target = [[-6,-6], [6,-6], [6,6], [-6,6]];
    faces = [
        [0,1,2,3],
        [4,5,6,7,8]
    ];
    star_xy = [for (i = [0:1:4]) [10*cos(-144*i), 10*sin(-144*i)]];
    verts_local = concat(
        [for (p = target) [p[0], p[1], 0]],
        [for (p = star_xy) [p[0], p[1], p[1] / 6]]
    );
    segs_evenodd = ps_face_geom_cut_segments(target, 0, faces, verts_local, 1e-8, "evenodd", true);
    segs_nonzero = ps_face_geom_cut_segments(target, 0, faces, verts_local, 1e-8, "nonzero", true);
    len_evenodd = norm(segs_evenodd[0][1] - segs_evenodd[0][0]);
    len_nonzero = norm(segs_nonzero[0][1] - segs_nonzero[0][0]);
    assert(len(segs_evenodd) > 0, "synthetic star cutter should generate some cut geometry");
    assert(len(segs_nonzero) > 0, "synthetic star cutter should generate some nonzero cut geometry");
    assert(segs_nonzero[0] != segs_evenodd[0], str("nonzero and evenodd cut segments should differ evenodd=", segs_evenodd, " nonzero=", segs_nonzero));
    assert(len_nonzero > len_evenodd + 1e-6, str("nonzero star cutter should span a longer merged cut evenodd=", len_evenodd, " nonzero=", len_nonzero));
}

module test_seg_merge_face_cut_group__preserves_disjoint_spans() {
    group = [
        [[[0, 0], [2, 0]], 7, 110],
        [[[4, 0], [6, 0]], 7, 120]
    ];
    merged = _ps_seg_merge_face_cut_group(group, 1e-8);
    assert_int_eq(len(merged), 2, str("disjoint same-face cut spans must remain separate merged=", merged));
}

module test_seg_merge_face_cut_group__merges_touching_spans() {
    group = [
        [[[0, 0], [2, 0]], 7, 110],
        [[[2, 0], [5, 0]], 7, 120]
    ];
    merged = _ps_seg_merge_face_cut_group(group, 1e-8);
    assert_int_eq(len(merged), 1, str("touching same-face cut spans should merge merged=", merged));
    assert(_ps_seg2_eq(merged[0][0], [[0, 0], [5, 0]], 1e-8), str("merged touching span endpoints wrong merged=", merged));
    assert(merged[0][2] == 120, str("merged touching span dihedral should keep max merged=", merged));
}

module test_ps_clip_face_by_feature_proxies__smoke_selected_indices() {
    p = hexahedron();

    ps_clip_face_by_feature_proxies(
        p,
        0,
        edge_len = 1,
        face_bounds = [-0.2, 0.2],
        face_proxy_mode = "sweep_to_bounds",
        edge_radius = 0.2,
        edge_length = 1.2,
        vertex_radius = 0.18,
        face_indices = [1],
        edge_indices = [0, 1, 2],
        vertex_indices = [0, 1, 2]
    ) {
        translate([0, 0, -0.1])
            cube([0.9, 0.9, 0.2], center = true);
        cube([0.8, 0.25, 0.25], center = true);
        sphere(r = 0.2);
    }

    assert(true, "proxy face smoke");
}

module test_ps_carve_face_by_feature_proxies__smoke_selected_indices() {
    p = hexahedron();

    ps_carve_face_by_feature_proxies(
        p,
        0,
        edge_len = 1,
        face_bounds = [-0.2, 0.2],
        face_proxy_mode = "sweep_to_bounds",
        edge_radius = 0.2,
        edge_length = undef,
        vertex_radius = 0.18,
        include_local_edges = true,
        include_local_vertices = false,
        include_cutter_faces = true,
        include_cutter_edges = false,
        include_cutter_vertices = false,
        cutter_face_indices = [1],
        local_edge_indices = [0, 1]
    ) {
        translate([0, 0, -0.1])
            cube([0.9, 0.9, 0.2], center = true);
        cube([0.8, 0.25, 0.25], center = true);
        sphere(r = 0.2);
        translate([0, 0, -0.1])
            cube([0.9, 0.9, 0.2], center = true);
        cube([0.8, 0.25, 0.25], center = true);
        sphere(r = 0.2);
    }

    assert(true, "proxy face carve smoke");
}

module test_ps_clip_face_local_clearance_by_feature_proxies__star_uses_boundary_edge_frames() {
    p = poly_antiprism(5, 2);

    ps_clip_face_local_clearance_by_feature_proxies(
        p,
        0,
        inter_radius = 20,
        face_bounds = [-0.2, 0.2],
        edge_radius = 0.2,
        edge_length = undef,
        vertex_radius = 0.18,
        include_local_edges = true,
        include_local_vertices = false,
        include_cutter_faces = false,
        include_cutter_edges = false,
        include_cutter_vertices = false,
        local_edge_indices = [0]
    ) {
        assert_int_eq($ps_edge_idx, 0, "star local clearance should honor source edge filter");
        assert(!is_undef($ps_dihedral), "star local clearance boundary frame should expose dihedral");
        assert($ps_edge_len > 0, "star local clearance boundary frame should have positive subsegment length");
        cube([0.5, 0.05, 0.05], center = true);
    }

    assert(true, "star local clearance smoke");
}

module test_ps_face_interference_volume_ctx__smoke() {
    p = poly_antiprism(n = 7, p = 3, angle = 15);

    place_on_faces(p, inter_radius = 20, indices = [2]) {
        ps_face_interference_volume_ctx(-0.4, 1.6);
    }

    place_on_faces(p, inter_radius = 20, indices = [0]) {
        ps_face_interference_volume_ctx(-0.4, 1.6);
    }

    assert(true, "face interference smoke");
}

module test_ps_partition_face_by_feature_proxies__smoke_selected_indices() {
    p = hexahedron();

    ps_partition_face_by_feature_proxies(
        p,
        0,
        edge_len = 1,
        face_bounds = [-0.2, 0.2],
        face_proxy_mode = "sweep_to_bounds",
        edge_radius = 0.2,
        edge_length = undef,
        vertex_radius = 0.18,
        include_faces = false,
        include_edges = true,
        include_vertices = false,
        edge_indices = [0, 1],
        color_cells = false,
        max_cutters = 4
    ) {
        translate([0, 0, -0.1])
            cube([0.9, 0.9, 0.2], center = true);
        cube([0.8, 0.25, 0.25], center = true);
        sphere(r = 0.2);
    }

    assert(true, "proxy face partition smoke");
}

module run_TestPlacement() {
    test_place_on_faces__family_ids_and_counts_from_classify();
    test_place_on_edges__family_ids_and_counts_from_classify();
    test_place_on_vertices__family_ids_and_counts_from_classify();
    test_place_on_faces__auto_classify_matches_precomputed();
    test_place_on_faces__indices_filters_exact_face();
    test_place_on_all__cube_single_family();
    test_place_on_edges__no_auto_classify_by_default();
    test_place_on_edges__indices_filters_exact_edge();
    test_place_on_edges__ez_world_matches_adjacent_face_bisector();
    test_place_on_vertices__indices_filters_exact_vertex();
    test_place_on_face_segments__star_face_split();
    test_place_on_face_segments__default_nonzero_keeps_filled_star();
    test_place_on_face_filled_boundary_segments__triangle_tracks_face_edges();
    test_place_on_face_filled_boundary_segments__pentagram_returns_decagon_segments();
    test_ps_face_filled_boundary_segments__triangle_returns_three_edges();
    test_ps_face_filled_boundary_segments__pentagram_returns_decagon_boundary();
    test_ps_face_filled_atoms__triangle_stays_single_boundary_atom();
    test_ps_face_filled_atoms__concave_hex_splits_and_preserves_area();
    test_place_on_faces__local_z_origin_consistent_for_face_and_poly_verts();
    test_seg_cycle_probe_point__concave_inside();
    test_seg_face_tris3__concave_area_preserved();
    test_seg_face_tris3__star_area_matches_segments();
    test_ps_face_visible_segments__cube_face_unchanged();
    test_ps_face_visible_segments__star_antiprism_side_reduced();
    test_ps_face_visible_segments__star_prism_side_keeps_two_side_quads();
    test_ps_face_visible_segments__cells_preserve_parent_winding();
    test_ps_face_visible_segments__star_prism_cut_edges_reference_cut_entries();
    test_ps_face_visible_segments__star_antiprism_cut_edges_reference_cut_entries();
    test_ps_face_visible_segments__cut_pair_ids_defined_in_ctx();
    test_ps_face_geom_cut_pair_ids__matching_faces_share_join_id();
    test_ps_face_visible_cell_cut_run_ids__wraparound_same_cutter_merges();
    test_ps_face_visible_segments__split_cut_spans_get_distinct_run_ids();
    test_ps_face_geom_cut_segments__respects_fill_mode();
    test_seg_merge_face_cut_group__preserves_disjoint_spans();
    test_seg_merge_face_cut_group__merges_touching_spans();
    test_ps_clip_face_by_feature_proxies__smoke_selected_indices();
    test_ps_carve_face_by_feature_proxies__smoke_selected_indices();
    test_ps_clip_face_local_clearance_by_feature_proxies__star_uses_boundary_edge_frames();
    test_ps_face_interference_volume_ctx__smoke();
    test_ps_partition_face_by_feature_proxies__smoke_selected_indices();
}

run_TestPlacement();
