use <../../polysymmetrica/core/placement.scad>
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

module run_TestPlacement() {
    test_place_on_faces__family_ids_and_counts_from_classify();
    test_place_on_edges__family_ids_and_counts_from_classify();
    test_place_on_vertices__family_ids_and_counts_from_classify();
    test_place_on_faces__auto_classify_matches_precomputed();
    test_place_on_all__cube_single_family();
    test_place_on_edges__no_auto_classify_by_default();
    test_place_on_face_segments__star_face_split();
    test_place_on_face_segments__default_nonzero_keeps_filled_star();
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
    test_ps_face_geom_cut_segments__respects_fill_mode();
    test_seg_merge_face_cut_group__preserves_disjoint_spans();
    test_seg_merge_face_cut_group__merges_touching_spans();
}

run_TestPlacement();
