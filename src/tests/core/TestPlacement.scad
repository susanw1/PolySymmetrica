use <../../polysymmetrica/core/placement.scad>
use <../../polysymmetrica/core/classify.scad>
use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/archimedians_all.scad>
use <../../polysymmetrica/models/tetrahedron.scad>

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_vec3_near(v, w, eps=1e-9, msg="") {
    assert(norm(v - w) <= eps, str(msg, " expected=", w, " got=", v));
}

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

module test_ps_face_sites__cube_records_match_face_structure() {
    p = hexahedron();
    faces = poly_faces(p);
    cls = poly_classify(p, 1, 1e-6, 1, false);
    ids = ps_classify_face_ids(cls, len(faces));
    counts = ps_classify_counts(cls);
    sites = ps_face_sites(p, classify = cls);

    assert_int_eq(len(sites), len(faces), "face site count should match face count");

    for (site = sites) {
        fi = site[0];
        face = faces[fi];

        assert_int_eq(site[6], len(face), "site vertex count should match face arity");
        assert_int_eq(len(site[10]), len(face), "site 2d point count should match face arity");
        assert_int_eq(len(site[11]), len(face), "site 3d local point count should match face arity");
        assert_int_eq(len(site[20]), len(face), "site neighbor count should match face arity");
        assert_int_eq(len(site[21]), len(face), "site dihedral count should match face arity");
        assert(site[15], str("cube face site should be planar fi=", fi));
        assert(abs(site[14]) < 1e-9, str("cube face planarity err should be ~0 fi=", fi, " err=", site[14]));
        assert_int_eq(site[16], ids[fi], "site face family id");
        assert_int_eq(site[17], counts[0], "site face family count");
        assert_int_eq(site[18], counts[1], "site edge family count");
        assert_int_eq(site[19], counts[2], "site vertex family count");
    }
}

module test_ps_edge_sites__cube_records_match_edge_structure() {
    p = hexahedron();
    faces = poly_faces(p);
    edges = _ps_edges_from_faces(faces);
    cls = poly_classify(p, 1, 1e-6, 1, false);
    ids = ps_classify_edge_ids(cls, len(edges));
    counts = ps_classify_counts(cls);
    sites = ps_edge_sites(p, classify = cls);

    assert_int_eq(len(sites), len(edges), "edge site count should match edge count");

    for (site = sites) {
        ei = site[0];

        assert_vec3_near(site[8][0], [-site[5] / 2, 0, 0], 1e-9, "edge local start point");
        assert_vec3_near(site[8][1], [site[5] / 2, 0, 0], 1e-9, "edge local end point");
        assert_int_eq(len(site[9]), 2, "edge vertex pair should have arity 2");
        assert_int_eq(len(site[10]), 2, "cube edges should have two adjacent faces");
        assert(abs(norm(site[1]) - site[6]) < 1e-9, str("edge center/midradius mismatch ei=", ei));
        assert_int_eq(site[11], ids[ei], "site edge family id");
        assert_int_eq(site[12], counts[0], "site face family count");
        assert_int_eq(site[13], counts[1], "site edge family count");
        assert_int_eq(site[14], counts[2], "site vertex family count");
    }
}

module test_ps_edge_sites__cube_uses_adjacent_face_normal_bisector() {
    p = hexahedron();
    verts = poly_verts(p);
    faces = poly_faces(p);
    faces0 = ps_orient_all_faces_outward(verts, faces);
    edges = _ps_edges_from_faces(faces0);
    edge_faces = ps_edge_faces_table(faces0, edges);
    face_n = [for (f = faces0) ps_face_normal(verts, f)];
    sites = ps_edge_sites(p);

    for (site = sites) {
        ei = site[0];
        adj = edge_faces[ei];
        nsum = face_n[adj[0]] + face_n[adj[1]];
        expected = v_norm(nsum);
        assert(
            v_dot(site[4], expected) > 1 - 1e-9,
            str("edge site ez should follow adjacent-face bisector ei=", ei, " got=", site[4], " expected=", expected)
        );
    }
}

module test_ps_edge_sites__preserves_raw_edge_order_for_classify_ids() {
    p = hexahedron();
    verts = poly_verts(p);
    faces_rev = [for (f = poly_faces(p)) [for (i = [len(f)-1:-1:0]) f[i]]];
    p_rev = [verts, faces_rev, poly_e_over_ir(p)];
    raw_edges = _ps_edges_from_faces(faces_rev);
    oriented_edges = _ps_edges_from_faces(ps_orient_all_faces_outward(verts, faces_rev));
    cls = poly_classify(p_rev, 1, 1e-6, 1, false);
    ids = ps_classify_edge_ids(cls, len(raw_edges));
    sites = ps_edge_sites(p_rev, classify = cls);

    assert(raw_edges != oriented_edges, "reversed input should change raw edge ordering for this regression test");
    assert_int_eq(len(sites), len(raw_edges), "reversed poly edge site count");

    for (site = sites) {
        ei = site[0];
        assert(site[9] == raw_edges[ei], str("edge site should preserve raw edge order ei=", ei, " got=", site[9], " expected=", raw_edges[ei]));
        assert_int_eq(site[11], ids[ei], "edge family id should match raw-order classify ids");
    }
}

module test_ps_vertex_sites__cube_records_match_vertex_structure() {
    p = hexahedron();
    verts = poly_verts(p);
    cls = poly_classify(p, 1, 1e-6, 1, false);
    ids = ps_classify_vert_ids(cls, len(verts));
    counts = ps_classify_counts(cls);
    sites = ps_vertex_sites(p, classify = cls);

    assert_int_eq(len(sites), len(verts), "vertex site count should match vertex count");

    for (site = sites) {
        vi = site[0];

        assert_int_eq(site[8], 3, "cube vertex valence should be 3");
        assert_int_eq(len(site[9]), site[8], "neighbor index count should match valence");
        assert_int_eq(len(site[10]), site[8], "neighbor point count should match valence");
        assert(abs(norm(site[1]) - site[6]) < 1e-9, str("vertex center/radius mismatch vi=", vi));
        assert_vec3_near(site[7], [0, 0, -site[6]], 1e-9, "vertex poly-center local");
        for (p_local = site[10])
            assert(abs(norm(p_local) - site[5]) < 1e-9, str("vertex neighbor edge length mismatch vi=", vi, " p=", p_local));
        assert_int_eq(site[11], ids[vi], "site vertex family id");
        assert_int_eq(site[12], counts[0], "site face family count");
        assert_int_eq(site[13], counts[1], "site edge family count");
        assert_int_eq(site[14], counts[2], "site vertex family count");
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
    tris = _ps_seg_face_tris3([0,1,2,3,4], pts3, 1e-9);
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

module test_ps_face_segments__default_matches_nonzero() {
    pts3 = [[0,9,0], [-5,-5,0], [8,3,0], [-8,3,0], [5,-5,0]];
    segs_default = ps_face_segments(pts3, eps=1e-9);
    segs_nonzero = ps_face_segments(pts3, "nonzero", 1e-9);
    assert(segs_default == segs_nonzero, "ps_face_segments default should match nonzero fill");
}

module test_ps_face_arrangement__pentagram_counts() {
    p = poly_antiprism(5, 2);
    pts3 = [for (i = poly_faces(p)[1]) poly_verts(p)[i]];
    arr = ps_face_arrangement(pts3, 1e-9);
    crossings = arr[1];
    nodes = arr[2];
    spans = arr[3];
    cells = arr[4];

    assert_int_eq(len(crossings), 5, "pentagram crossings");
    assert_int_eq(len(nodes), 10, "pentagram arrangement nodes");
    assert_int_eq(len(spans), 15, "pentagram arrangement spans");
    assert_int_eq(len(cells), 7, "pentagram arrangement cells");

    assert_int_eq(len([for (n = nodes) if (n[1] == "source_vertex") 1]), 5, "pentagram source-vertex count");
    assert_int_eq(len([for (n = nodes) if (n[1] == "crossing") 1]), 5, "pentagram crossing-node count");
    assert_int_eq(len([for (s = spans) if (s[6] == "source") 1]), len(spans), "pentagram span kinds");
}

module test_ps_face_boundary_model__pentagram_counts() {
    p = poly_antiprism(5, 2);
    pts3 = [for (i = poly_faces(p)[1]) poly_verts(p)[i]];
    pts2 = [for (pt = pts3) [pt[0], pt[1]]];
    nz = ps_face_boundary_model(pts3, "nonzero", 1e-9);
    eo = ps_face_boundary_model(pts3, "evenodd", 1e-9);

    assert_int_eq(len(nz[1]), 1, "pentagram nonzero filled cell count");
    assert_int_eq(len(nz[2]), 1, "pentagram nonzero boundary loop count");
    assert_int_eq(len(nz[3]), 10, "pentagram nonzero boundary span count");

    assert_int_eq(len(eo[1]), 5, "pentagram evenodd filled cell count");
    assert_int_eq(len(eo[2]), 2, "pentagram evenodd boundary loop count");
    assert_int_eq(len(eo[3]), 15, "pentagram evenodd boundary span count");

    for (bm = [nz, eo])
        for (span = bm[3]) {
            edge_idx = span[2];
            t0 = span[3];
            t1 = span[4];
            a = pts2[edge_idx];
            b = pts2[(edge_idx + 1) % len(pts2)];
            p0 = [a[0] + (b[0] - a[0]) * t0, a[1] + (b[1] - a[1]) * t0];
            p1 = [a[0] + (b[0] - a[0]) * t1, a[1] + (b[1] - a[1]) * t1];
            assert(norm(p0 - span[0][0]) < 1e-6, str("boundary span start should match source params edge=", edge_idx, " t0=", t0));
            assert(norm(p1 - span[0][1]) < 1e-6, str("boundary span end should match source params edge=", edge_idx, " t1=", t1));
        }
}

module test_ps_face_boundary_span_sites__pentagram_attach_adjacent_face_context() {
    p = poly_antiprism(5, 2);
    place_on_faces(p) {
        if ($ps_face_idx == 1) {
            sites = _ps_face_boundary_span_sites(
                $ps_face_pts3d_local,
                $ps_face_idx,
                $ps_poly_faces_idx,
                $ps_poly_verts_local,
                $ps_face_neighbors_idx,
                $ps_face_dihedrals,
                "nonzero",
                1e-9
            );

            assert_int_eq(len(sites), 10, "pentagram nonzero boundary span site count");

            for (site = sites) {
                ei = site[8];
                assert_int_eq(site[14], $ps_face_neighbors_idx[ei], "boundary span adjacent face id");
                assert(abs(site[15] - $ps_face_dihedrals[ei]) < 1e-6, str("boundary span dihedral mismatch edge=", ei));
                assert(!is_undef(site[16]), str("boundary span adjacent-face normal should be defined edge=", ei));
                assert(site[17] != 0, str("boundary span filled side should be nonzero edge=", ei));
                assert(!is_undef(site[18]), str("boundary span adjacent-face direction should be defined edge=", ei));
                assert(abs(site[18][0]) < 1e-6, str("boundary span adjacent-face direction should stay in local yz plane edge=", ei));
                assert(site[18][2] > 0, str("boundary span adjacent-face direction should point to current-face +Z edge=", ei));
            }
        }
    }
}

module test_ps_face_boundary_span_direction__projects_source_edge_into_face_plane() {
    ex = [1, 0, 0];
    ey = [0, 1, 0];
    ez = [0, 0, 1];
    source_ex_raw = v_norm([1, 0, 1]);
    source_ex_proj = v_norm(_ps_seg_project_to_plane(source_ex_raw, ez));
    adj_face_normal_local = v_norm([0, 1, 1]);
    dir_span_local = _ps_seg_boundary_span_adj_face_dir_span_local(source_ex_proj, ex, ey, ez, adj_face_normal_local);

    assert(abs(source_ex_proj[2]) < 1e-9, "projected source edge should lie in the current face plane");
    assert(abs(dir_span_local[0]) < 1e-9, "adjacent-face direction should stay in local yz plane after source-edge projection");
    assert(dir_span_local[2] > 0, "adjacent-face direction should keep the +Z branch");
}

module test_ps_face_boundary_span_sites__anti_tet_hex_is_span_directional() {
    p = poly_truncate(tetrahedron(), t = -0.5);
    place_on_faces(p) {
        if ($ps_face_idx == 0) {
            sites = _ps_face_boundary_span_sites(
                $ps_face_pts3d_local,
                $ps_face_idx,
                $ps_poly_faces_idx,
                $ps_poly_verts_local,
                $ps_face_neighbors_idx,
                $ps_face_dihedrals,
                "nonzero",
                1e-9
            );
            source_edges = len($ps_face_pts2d);
            repeated_edges = [
                for (ei = [0:1:source_edges-1])
                    if (len([for (s = sites) if (s[8] == ei) 1]) > 1)
                        ei
            ];
            mixed_dir_edges = [
                for (ei = repeated_edges)
                    let(
                        has_inc = len([for (s = sites) if (s[8] == ei && s[10] > s[9]) 1]) > 0,
                        has_dec = len([for (s = sites) if (s[8] == ei && s[10] < s[9]) 1]) > 0
                    )
                    if (has_inc && has_dec)
                        ei
            ];

            assert_int_eq(len(sites), 12, "anti-tet hex nonzero boundary span site count");
            assert(len(repeated_edges) > 0, "anti-tet hex should reuse source edges across multiple boundary spans");
            assert(len(mixed_dir_edges) > 0, "anti-tet hex should need per-span source-edge directionality");

            for (site = sites) {
                ei = site[8];
                assert_int_eq(site[14], $ps_face_neighbors_idx[ei], "anti-tet boundary span adjacent face id");
                assert(abs(site[15] - $ps_face_dihedrals[ei]) < 1e-6, str("anti-tet boundary span dihedral mismatch edge=", ei));
                assert(site[17] != 0, str("anti-tet boundary span filled side should be nonzero edge=", ei));
                assert(!is_undef(site[18]), str("anti-tet boundary span adjacent-face direction should be defined edge=", ei));
                assert(abs(site[18][0]) < 1e-6, str("anti-tet boundary span adjacent-face direction should stay in local yz plane edge=", ei));
                assert(site[18][2] > 0, str("anti-tet boundary span adjacent-face direction should point to current-face +Z edge=", ei));
            }

            se1_span_dirs = [
                for (site = sites)
                    if (site[8] == 1)
                        site[18]
            ];
            se1_has_pos_y = len([for (dir = se1_span_dirs) if (dir[1] > 0) 1]) > 0;
            se1_has_neg_y = len([for (dir = se1_span_dirs) if (dir[1] < 0) 1]) > 0;
            assert(se1_has_pos_y && se1_has_neg_y, "anti-tet spans from one long source edge should distinguish central vs end branches");
        }
    }
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
    assert(len(segs_evenodd) > 0, "synthetic star cutter should generate some cut geometry");
    assert(len(segs_nonzero) > len(segs_evenodd), str("nonzero star cutter should yield more cut segments than evenodd evenodd=", len(segs_evenodd), " nonzero=", len(segs_nonzero)));
}

module run_TestPlacement() {
    test_place_on_faces__family_ids_and_counts_from_classify();
    test_place_on_edges__family_ids_and_counts_from_classify();
    test_place_on_vertices__family_ids_and_counts_from_classify();
    test_place_on_faces__auto_classify_matches_precomputed();
    test_ps_face_sites__cube_records_match_face_structure();
    test_ps_edge_sites__cube_records_match_edge_structure();
    test_ps_edge_sites__cube_uses_adjacent_face_normal_bisector();
    test_ps_edge_sites__preserves_raw_edge_order_for_classify_ids();
    test_ps_vertex_sites__cube_records_match_vertex_structure();
    test_place_on_all__cube_single_family();
    test_place_on_edges__no_auto_classify_by_default();
    test_place_on_face_segments__star_face_split();
    test_place_on_faces__local_z_origin_consistent_for_face_and_poly_verts();
    test_seg_cycle_probe_point__concave_inside();
    test_seg_face_tris3__concave_area_preserved();
    test_seg_face_tris3__star_area_matches_segments();
    test_ps_face_segments__default_matches_nonzero();
    test_ps_face_arrangement__pentagram_counts();
    test_ps_face_boundary_model__pentagram_counts();
    test_ps_face_boundary_span_sites__pentagram_attach_adjacent_face_context();
    test_ps_face_boundary_span_direction__projects_source_edge_into_face_plane();
    test_ps_face_boundary_span_sites__anti_tet_hex_is_span_directional();
    test_ps_face_visible_segments__cube_face_unchanged();
    test_ps_face_visible_segments__star_antiprism_side_reduced();
    test_ps_face_visible_segments__cells_preserve_parent_winding();
    test_ps_face_geom_cut_segments__respects_fill_mode();
}

run_TestPlacement();
