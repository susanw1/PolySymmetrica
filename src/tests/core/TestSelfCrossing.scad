use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/placement.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/core/segments.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/models/tetrahedron.scad>

EPS = 1e-8;
MODE = "nonzero";
STAR_FACE_IDX = 1;
TRI_FACE_IDX = 12;
ANTI_FACE_IDX = 0;

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_list_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a - b) <= eps, str(msg, " expected=", b, " got=", a));
}

function _test_punch_poly() =
    poly_antiprism(7, 3, angle = 15);

function _test_punch_poly_angle0() =
    poly_antiprism(7, 3, angle = 0);

function _test_antitet_poly() =
    poly_truncate(tetrahedron(), t = -0.5);

function _test_face_site(poly, face_idx) =
    ps_face_sites(poly)[face_idx];

function _test_source_counts(records, source_idx_pos, n) =
    [
        for (ei = [0:1:n-1])
            len([for (r = records) if (r[source_idx_pos] == ei) 1])
    ];

function _test_replay_kind_count(sites, kind) =
    len([for (s = sites) if (ps_replay_site_foreign_kind(s) == kind) 1]);

function _test_coincident_intrusion_verts_local() =
    [
        [-2, -2, 0], [2, -2, 0], [2, 2, 0], [-2, 2, 0],
        [0, -2, -1], [0, 0, 1], [0, 2, -1],
        [0, -2, -1], [0, 0, 1], [0, 2, -1]
    ];

function _test_coincident_intrusion_faces_idx() =
    [
        [0, 1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ];

module test_ps_face_arrangement__7_3_15_star_has_stable_structure() {
    site = _test_face_site(_test_punch_poly(), STAR_FACE_IDX);
    arr = ps_face_arrangement(site[11]);

    assert_int_eq(len(arr[1]), 14, "star face crossing count");
    assert_int_eq(len(arr[2]), 21, "star face arrangement node count");
    assert_int_eq(len(arr[3]), 35, "star face arrangement span count");
    assert_int_eq(len(arr[4]), 16, "star face arrangement cell count");
    assert_list_eq(
        _test_source_counts(arr[3], 3, len(site[10])),
        [5, 5, 5, 5, 5, 5, 5],
        "star face split spans should distribute evenly across source edges"
    );
    assert_list_eq(
        [for (c = arr[1]) [c[0], c[2]]],
        [[0, 2], [0, 3], [0, 4], [0, 5], [1, 3], [1, 4], [1, 5], [1, 6], [2, 4], [2, 5], [2, 6], [3, 5], [3, 6], [4, 6]],
        "star face crossing source-edge pairs"
    );
}

module test_ps_face_boundary_model__7_3_15_star_has_true_nonzero_boundary() {
    site = _test_face_site(_test_punch_poly(), STAR_FACE_IDX);
    bm = ps_face_boundary_model(site[11], MODE);
    segments = ps_face_segments(site[11], MODE);

    assert_int_eq(len(bm[1]), 1, "star face nonzero filled cell count");
    assert_int_eq(len(bm[2]), 1, "star face nonzero boundary loop count");
    assert_int_eq(len(bm[3]), 14, "star face nonzero boundary span count");
    assert_list_eq(
        _test_source_counts(bm[3], 2, len(site[10])),
        [2, 2, 2, 2, 2, 2, 2],
        "star face boundary spans should distribute evenly across source edges"
    );
    assert_int_eq(len(segments), 1, "star face nonzero should produce one visible filled segment");
    assert_list_eq(
        segments[0][2],
        [0, 4, 5, 2, 3, 0, 1, 5, 6, 3, 4, 1, 2, 6],
        "star face filled segment source-edge lineage"
    );
}

module test_ps_face_filled_boundary_source_edges__7_3_15_star_groups_surviving_spans() {
    site = _test_face_site(_test_punch_poly(), STAR_FACE_IDX);
    source_edges = ps_face_filled_boundary_source_edges(site[11], MODE);

    assert_int_eq(len(source_edges), 7, "star face should expose one filled-boundary record per source edge");
    assert_list_eq([for (e = source_edges) e[0]], [0, 1, 2, 3, 4, 5, 6], "star face source-edge ids");
    assert_list_eq([for (e = source_edges) len(e[2])], [2, 2, 2, 2, 2, 2, 2], "star face surviving span count per source edge");
    assert_list_eq(
        [for (e = source_edges) [for (span = e[2]) span[8]]],
        [[-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1]],
        "star face filled side per surviving span"
    );
    assert_near(source_edges[0][2][0][3], 0, EPS, "star face source edge 0 first span t0");
    assert_near(source_edges[0][2][0][4], 0.356896, 1e-6, "star face source edge 0 first span t1");
    assert_near(source_edges[0][2][1][3], 0.643104, 1e-6, "star face source edge 0 second span t0");
    assert_near(source_edges[0][2][1][4], 1, EPS, "star face source edge 0 second span t1");
}

module test_ps_face_geom_cut_entries__7_3_15_triangle_records_foreign_cutters() {
    site = _test_face_site(_test_punch_poly(), TRI_FACE_IDX);
    cuts = ps_face_geom_cut_entries(site[10], site[0], site[13], site[12], mode = MODE, filter_parent = true);

    assert_int_eq(len(cuts), 6, "triangle punch-through cut count");
    assert_list_eq(
        [for (c = cuts) c[1]],
        [3, 8, 9, 10, 14, 15],
        "triangle punch-through cutter face ids"
    );
    assert(
        min([for (c = cuts) c[2]]) > 90 && max([for (c = cuts) c[2]]) < 140,
        str("triangle punch-through cut dihedrals should stay in expected range: ", [for (c = cuts) c[2]])
    );
}

module test_ps_face_foreign_intrusion_records__7_3_15_triangle_wraps_exact_face_cuts() {
    site = _test_face_site(_test_punch_poly(), TRI_FACE_IDX);
    records = ps_face_foreign_intrusion_records(site[10], site[0], site[13], site[12], mode = MODE, filter_parent = true);

    assert_int_eq(len(records), 6, "triangle punch-through intrusion count");
    assert_list_eq(
        [for (r = records) ps_intrusion_kind(r)],
        ["face_plane_cut", "face_plane_cut", "face_plane_cut", "face_plane_cut", "face_plane_cut", "face_plane_cut"],
        "triangle intrusion kinds"
    );
    assert_list_eq(
        [for (r = records) ps_intrusion_target_face_idx(r)],
        [TRI_FACE_IDX, TRI_FACE_IDX, TRI_FACE_IDX, TRI_FACE_IDX, TRI_FACE_IDX, TRI_FACE_IDX],
        "triangle intrusion target face ids"
    );
    assert_list_eq(
        [for (r = records) ps_intrusion_foreign_kind(r)],
        ["face", "face", "face", "face", "face", "face"],
        "triangle intrusion foreign kinds"
    );
    assert_list_eq(
        [for (r = records) ps_intrusion_foreign_idx(r)],
        [3, 8, 9, 10, 14, 15],
        "triangle intrusion foreign face ids"
    );
    assert_list_eq(
        [for (r = records) ps_intrusion_confidence(r)],
        ["exact", "exact", "exact", "exact", "exact", "exact"],
        "triangle intrusion confidence"
    );
    assert_list_eq([for (r = records) len(ps_intrusion_segment2d_local(r))], [2, 2, 2, 2, 2, 2], "triangle intrusion segment arities");
    assert(
        min([for (r = records) ps_intrusion_dihedral(r)]) > 90 && max([for (r = records) ps_intrusion_dihedral(r)]) < 140,
        str("triangle intrusion dihedrals should stay in expected range: ", [for (r = records) ps_intrusion_dihedral(r)])
    );
}

module test_ps_face_foreign_intrusion_records__preserves_coincident_foreign_face_provenance() {
    face_pts2d = [[-2, -2], [2, -2], [2, 2], [-2, 2]];
    records = ps_face_foreign_intrusion_records(
        face_pts2d,
        0,
        _test_coincident_intrusion_faces_idx(),
        _test_coincident_intrusion_verts_local(),
        mode = MODE,
        filter_parent = true
    );

    assert_int_eq(len(records), 2, "coincident cuts from distinct foreign faces should both survive");
    assert_list_eq([for (r = records) ps_intrusion_foreign_idx(r)], [1, 2], "coincident intrusion foreign face ids");
    assert_list_eq(
        [for (r = records) ps_intrusion_segment2d_local(r)],
        [ps_intrusion_segment2d_local(records[0]), ps_intrusion_segment2d_local(records[0])],
        "coincident intrusion segments remain geometrically identical"
    );
}

module test_ps_face_foreign_face_replay_sites__7_3_15_triangle_builds_target_local_frames() {
    site = _test_face_site(_test_punch_poly(), TRI_FACE_IDX);
    replay = ps_face_foreign_face_replay_sites(site[10], site[0], site[13], site[12], site[9], mode = MODE, filter_parent = true);

    assert_int_eq(len(replay), 6, "triangle foreign face replay site count");
    assert_list_eq(
        [for (s = replay) ps_replay_site_foreign_idx(s)],
        [3, 8, 9, 10, 14, 15],
        "triangle replay foreign face ids"
    );
    assert_list_eq(
        [for (s = replay) ps_replay_site_foreign_kind(s)],
        ["face", "face", "face", "face", "face", "face"],
        "triangle replay foreign kinds"
    );
    assert_list_eq(
        [for (s = replay) len(ps_replay_site_intrusion_segment2d_local(s))],
        [2, 2, 2, 2, 2, 2],
        "triangle replay keeps target-local intrusion segment"
    );
    assert_list_eq(
        [for (s = replay) len(ps_replay_site_face_pts2d(s))],
        [3, 3, 3, 3, 3, 3],
        "triangle replay face point arities"
    );

    for (s = replay) {
        assert_near(norm(ps_replay_site_ex_local(s)), 1, EPS, "replay ex is unit");
        assert_near(norm(ps_replay_site_ey_local(s)), 1, EPS, "replay ey is unit");
        assert_near(norm(ps_replay_site_ez_local(s)), 1, EPS, "replay ez is unit");
        assert_near(v_dot(ps_replay_site_ex_local(s), ps_replay_site_ey_local(s)), 0, EPS, "replay ex/ey orthogonal");
        assert_near(v_dot(ps_replay_site_ex_local(s), ps_replay_site_ez_local(s)), 0, EPS, "replay ex/ez orthogonal");
        assert_near(v_dot(ps_replay_site_ey_local(s), ps_replay_site_ez_local(s)), 0, EPS, "replay ey/ez orthogonal");
        assert(max([for (p = ps_replay_site_face_pts3d_local(s)) abs(p[2])]) <= 1e-6, "planar replay face lies near local z=0");
    }
}

module test_ps_face_foreign_proxy_replay_sites__7_3_15_triangle_includes_edge_and_vertex_candidates() {
    site = _test_face_site(_test_punch_poly(), TRI_FACE_IDX);
    replay = ps_face_foreign_proxy_replay_sites(site[10], site[0], site[13], site[12], site[9], mode = MODE, filter_parent = true);

    assert_int_eq(_test_replay_kind_count(replay, "face"), 6, "triangle proxy replay exact face count");
    assert(_test_replay_kind_count(replay, "edge") > 0, "triangle proxy replay should include edge candidates");
    assert(_test_replay_kind_count(replay, "vertex") > 0, "triangle proxy replay should include vertex candidates");

    for (s = replay) {
        assert_near(norm(ps_replay_site_ex_local(s)), 1, EPS, "proxy replay ex is unit");
        assert_near(norm(ps_replay_site_ey_local(s)), 1, EPS, "proxy replay ey is unit");
        assert_near(norm(ps_replay_site_ez_local(s)), 1, EPS, "proxy replay ez is unit");
        if (ps_replay_site_foreign_kind(s) == "face")
            assert(ps_replay_site_intrusion_confidence(s) == "exact", "face replay confidence");
        if (ps_replay_site_foreign_kind(s) != "face")
            assert(ps_replay_site_intrusion_confidence(s) == "candidate", "edge/vertex replay confidence");
    }
}

module test_ps_face_visible_segments__7_3_15_triangle_splits_into_visible_cells() {
    site = _test_face_site(_test_punch_poly(), TRI_FACE_IDX);
    visible = ps_face_visible_segments(site[10], site[0], site[13], site[12], mode = MODE, filter_parent = true);

    assert_int_eq(len(visible), 3, "triangle punch-through visible segment count");
    assert_list_eq([for (seg = visible) len(seg[0])], [8, 3, 3], "triangle visible segment arities");
    assert_list_eq(
        [for (seg = visible) seg[3]],
        [
            ["parent", "parent", "parent", "cut", "cut", "cut", "cut", "cut"],
            ["parent", "cut", "cut"],
            ["parent", "parent", "cut"]
        ],
        "triangle visible segment edge-kind lineage"
    );
}

module test_ps_face_visible_segments__7_3_0_triangle_catches_meeting_cut_edges() {
    site = _test_face_site(_test_punch_poly_angle0(), TRI_FACE_IDX);
    cuts = ps_face_geom_cut_entries(site[10], site[0], site[13], site[12], mode = MODE, filter_parent = true);
    visible = ps_face_visible_segments(site[10], site[0], site[13], site[12], mode = MODE, filter_parent = true);

    assert_int_eq(len(cuts), 6, "angle=0 triangle punch-through cut count");
    assert_list_eq(
        [for (c = cuts) c[1]],
        [3, 8, 9, 10, 14, 15],
        "angle=0 triangle punch-through cutter face ids"
    );
    assert_int_eq(len(visible), 5, "angle=0 triangle should split into five visible cells");
    assert_list_eq([for (seg = visible) len(seg[0])], [4, 3, 3, 3, 3], "angle=0 triangle visible segment arities");
    assert_list_eq(
        [for (seg = visible) seg[3]],
        [
            ["parent", "parent", "cut", "cut"],
            ["cut", "parent", "cut"],
            ["parent", "cut", "cut"],
            ["parent", "parent", "cut"],
            ["cut", "parent", "parent"]
        ],
        "angle=0 triangle visible segment edge-kind lineage"
    );
}

module test_ps_face_filled_boundary_source_edges__7_3_0_triangle_is_simple_boundary() {
    site = _test_face_site(_test_punch_poly_angle0(), TRI_FACE_IDX);
    source_edges = ps_face_filled_boundary_source_edges(site[11], MODE);

    assert_int_eq(len(source_edges), 3, "simple triangle should expose three filled-boundary source edges");
    assert_list_eq([for (e = source_edges) e[0]], [0, 1, 2], "simple triangle source-edge ids");
    assert_list_eq([for (e = source_edges) len(e[2])], [1, 1, 1], "simple triangle surviving span count per source edge");
    assert_list_eq(
        [for (e = source_edges) [for (span = e[2]) span[8]]],
        [[-1], [-1], [-1]],
        "simple triangle filled side per source edge"
    );
    assert_list_eq(
        [for (e = source_edges) [for (span = e[2]) [span[3], span[4]]]],
        [[[0, 1]], [[0, 1]], [[0, 1]]],
        "simple triangle boundary span ranges are oriented to boundary traversal"
    );
}

module test_place_on_face_filled_boundary_source_edges__7_3_15_star_exposes_context() {
    place_on_faces(_test_punch_poly()) {
        if ($ps_face_idx == STAR_FACE_IDX) {
            place_on_face_filled_boundary_source_edges(MODE) {
                assert_int_eq($ps_boundary_source_edge_count, 7, "placed source-edge record count");
                assert_int_eq($ps_boundary_source_edge_span_count, 2, "placed source-edge span count");
                assert(
                    $ps_boundary_source_edge_idx >= 0 && $ps_boundary_source_edge_idx < 7,
                    str("placed source-edge idx in range: ", $ps_boundary_source_edge_idx)
                );
                assert_int_eq(
                    len($ps_boundary_source_edge_boundary_span_idxs),
                    $ps_boundary_source_edge_span_count,
                    "placed source-edge boundary-span id arity"
                );
                assert_int_eq(
                    len($ps_boundary_source_edge_sides),
                    $ps_boundary_source_edge_span_count,
                    "placed source-edge filled-side arity"
                );
                assert_int_eq(
                    $ps_boundary_source_edge_frame_side,
                    -1,
                    "placed source-edge frame normalizes filled side to local right"
                );
                assert_int_eq(
                    len($ps_boundary_source_edge_span_t_ranges_local),
                    $ps_boundary_source_edge_span_count,
                    "placed source-edge frame-local t-range arity"
                );
                assert_int_eq(
                    len($ps_boundary_source_edge_span_sides_local),
                    $ps_boundary_source_edge_span_count,
                    "placed source-edge frame-local filled-side arity"
                );
            }
        }
    }
}

module test_place_on_face_filled_boundary_source_edges__antitet_uses_span_direction() {
    place_on_faces(_test_antitet_poly()) {
        if ($ps_face_idx == ANTI_FACE_IDX) {
            place_on_face_filled_boundary_source_edges(MODE) {
                if ($ps_boundary_source_edge_idx == 1) {
                    assert_list_eq(
                        $ps_boundary_source_edge_span_sides_local,
                        [-1, 1, 1],
                        "antitet long source edge has middle/end spans on opposite frame sides"
                    );
                }
                if ($ps_boundary_source_edge_idx == 0) {
                    assert_list_eq(
                        $ps_boundary_source_edge_span_sides_local,
                        [-1],
                        "antitet short end source edge normalizes from source-param direction"
                    );
                }
            }
        }
    }
}

module test_place_on_face_foreign_intrusions__7_3_15_triangle_exposes_context() {
    place_on_faces(_test_punch_poly()) {
        if ($ps_face_idx == TRI_FACE_IDX) {
            place_on_face_foreign_intrusions(mode = MODE) {
                assert_int_eq($ps_intrusion_count, 6, "triangle intrusion iterator count");
                assert($ps_intrusion_idx >= 0 && $ps_intrusion_idx < $ps_intrusion_count, "triangle intrusion iterator idx bounds");
                assert_int_eq($ps_intrusion_target_face_idx, TRI_FACE_IDX, "triangle intrusion iterator target face id");
                assert($ps_intrusion_kind == "face_plane_cut", "triangle intrusion iterator kind");
                assert($ps_intrusion_foreign_kind == "face", "triangle intrusion iterator foreign kind");
                assert($ps_intrusion_confidence == "exact", "triangle intrusion iterator confidence");
                assert_int_eq(len($ps_intrusion_segment2d_local), 2, "triangle intrusion iterator segment arity");
            }
        }
    }
}

module test_place_on_face_foreign_face_replay_sites__7_3_15_triangle_exposes_context() {
    place_on_faces(_test_punch_poly()) {
        if ($ps_face_idx == TRI_FACE_IDX) {
            place_on_face_foreign_face_replay_sites(mode = MODE, coords = "parent") {
                assert_int_eq($ps_replay_count, 6, "triangle replay iterator count");
                assert($ps_replay_idx >= 0 && $ps_replay_idx < $ps_replay_count, "triangle replay iterator idx bounds");
                assert($ps_replay_kind == "foreign_face", "triangle replay iterator kind");
                assert($ps_replay_foreign_kind == "face", "triangle replay iterator foreign kind");
                assert($ps_replay_intrusion_confidence == "exact", "triangle replay iterator confidence");
                assert_int_eq(len($ps_replay_intrusion_segment2d_local), 2, "triangle replay intrusion segment arity");
                assert_int_eq(len($ps_replay_face_pts2d), 3, "triangle replay face point arity");
                assert_near(norm($ps_replay_ex_local), 1, EPS, "triangle replay ex unit");
                assert_near(norm($ps_replay_ey_local), 1, EPS, "triangle replay ey unit");
                assert_near(norm($ps_replay_ez_local), 1, EPS, "triangle replay ez unit");
            }
        }
    }
}

module _test_assert_triangle_proxy_face_child(expected_child_idx) {
    assert($ps_proxy_count > 6, str("triangle proxy iterator should include face plus candidates, count=", $ps_proxy_count));
    assert($ps_proxy_idx >= 0 && $ps_proxy_idx < $ps_proxy_count, "triangle proxy iterator idx bounds");
    assert($ps_proxy_kind == "foreign_face", "triangle proxy kind");
    assert($ps_proxy_source_kind == "face", "triangle proxy source kind");
    assert_int_eq($ps_proxy_child_idx, expected_child_idx, "triangle face proxy child slot");
    assert_int_eq($ps_proxy_target_face_idx, TRI_FACE_IDX, "triangle proxy target face id");
    assert_int_eq(len($ps_proxy_intrusion_segment2d_local), 2, "triangle proxy intrusion segment arity");
    assert($ps_proxy_intrusion_confidence == "exact", "triangle proxy intrusion confidence");
    assert_int_eq(len($ps_proxy_face_pts2d), 3, "triangle proxy face point arity");
    assert_near(norm($ps_proxy_ex_local), 1, EPS, "triangle proxy ex unit");
    assert_near(norm($ps_proxy_ey_local), 1, EPS, "triangle proxy ey unit");
    assert_near(norm($ps_proxy_ez_local), 1, EPS, "triangle proxy ez unit");
}

module _test_assert_triangle_proxy_edge_child(expected_child_idx) {
    assert($ps_proxy_idx >= 0 && $ps_proxy_idx < $ps_proxy_count, "triangle edge proxy idx bounds");
    assert($ps_proxy_kind == "foreign_edge", "triangle edge proxy kind");
    assert($ps_proxy_source_kind == "edge", "triangle edge proxy source kind");
    assert_int_eq($ps_proxy_child_idx, expected_child_idx, "triangle edge proxy child slot");
    assert_int_eq($ps_proxy_target_face_idx, TRI_FACE_IDX, "triangle edge proxy target face id");
    assert($ps_proxy_intrusion_confidence == "candidate", "triangle edge proxy confidence");
    assert_int_eq($ps_edge_idx, $ps_proxy_source_idx, "triangle edge proxy child edge id");
    assert_list_eq($ps_edge_pts_local, $ps_proxy_edge_pts_local, "triangle edge proxy child edge points");
    assert_list_eq($ps_edge_verts_idx, $ps_proxy_edge_verts_idx, "triangle edge proxy child edge vertices");
    assert_near(norm($ps_proxy_ex_local), 1, EPS, "triangle edge proxy ex unit");
    assert_near(norm($ps_proxy_ey_local), 1, EPS, "triangle edge proxy ey unit");
    assert_near(norm($ps_proxy_ez_local), 1, EPS, "triangle edge proxy ez unit");
}

module _test_assert_triangle_proxy_vertex_child(expected_child_idx) {
    assert($ps_proxy_idx >= 0 && $ps_proxy_idx < $ps_proxy_count, "triangle vertex proxy idx bounds");
    assert($ps_proxy_kind == "foreign_vertex", "triangle vertex proxy kind");
    assert($ps_proxy_source_kind == "vertex", "triangle vertex proxy source kind");
    assert_int_eq($ps_proxy_child_idx, expected_child_idx, "triangle vertex proxy child slot");
    assert_int_eq($ps_proxy_target_face_idx, TRI_FACE_IDX, "triangle vertex proxy target face id");
    assert($ps_proxy_intrusion_confidence == "candidate", "triangle vertex proxy confidence");
    assert_int_eq($ps_vertex_idx, $ps_proxy_source_idx, "triangle vertex proxy child vertex id");
    assert_int_eq($ps_vertex_valence, $ps_proxy_vertex_valence, "triangle vertex proxy child valence");
    assert_list_eq($ps_vertex_neighbors_idx, $ps_proxy_vertex_neighbors_idx, "triangle vertex proxy child neighbors");
    assert_near(norm($ps_proxy_ex_local), 1, EPS, "triangle vertex proxy ex unit");
    assert_near(norm($ps_proxy_ey_local), 1, EPS, "triangle vertex proxy ey unit");
    assert_near(norm($ps_proxy_ez_local), 1, EPS, "triangle vertex proxy ez unit");
}

module _test_assert_triangle_proxy_face_child_element_context(expected_child_idx) {
    _test_assert_triangle_proxy_face_child(expected_child_idx);

    assert_int_eq($ps_face_idx, $ps_proxy_source_idx, "triangle proxy child face id should be foreign face id");
    assert_int_eq(len($ps_face_pts2d), len($ps_proxy_face_pts2d), "triangle proxy child face point arity");
    assert_list_eq($ps_face_pts2d, $ps_proxy_face_pts2d, "triangle proxy child face points should match proxy face points");
    assert_list_eq($ps_face_pts3d_local, $ps_proxy_face_pts3d_local, "triangle proxy child local 3d face points should match proxy face points");
    assert_list_eq($ps_poly_faces_idx[$ps_face_idx], $ps_proxy_face_verts_idx, "triangle proxy child face vertex ids should match proxy source face");
    assert_int_eq(len($ps_face_neighbors_idx), len($ps_face_pts2d), "triangle proxy child neighbor arity");
    assert_int_eq(len($ps_face_dihedrals), len($ps_face_pts2d), "triangle proxy child dihedral arity");
}

module test_place_on_face_foreign_proxy_sites__7_3_15_triangle_dispatches_face_child() {
    place_on_faces(_test_punch_poly()) {
        if ($ps_face_idx == TRI_FACE_IDX) {
            place_on_face_foreign_proxy_sites(mode = MODE, coords = "parent") {
                _test_assert_triangle_proxy_face_child(0);
            }

            place_on_face_foreign_proxy_sites(mode = MODE, coords = "parent", face_child = 1, edge_child = 0, vertex_child = 2) {
                assert($ps_proxy_source_kind == "edge", "remapped child 0 should receive edge proxies");
                _test_assert_triangle_proxy_face_child(1);
                assert($ps_proxy_source_kind == "vertex", "remapped child 2 should receive vertex proxies");
            }
        }
    }
}

module test_place_on_face_foreign_proxy_sites__7_3_15_triangle_element_child_uses_source_face_context() {
    place_on_faces(_test_punch_poly()) {
        if ($ps_face_idx == TRI_FACE_IDX) {
            place_on_face_foreign_proxy_sites(mode = MODE) {
                _test_assert_triangle_proxy_face_child_element_context(0);
            }
        }
    }
}

module test_place_on_face_foreign_proxy_sites__7_3_15_triangle_dispatches_edge_and_vertex_children() {
    place_on_faces(_test_punch_poly()) {
        if ($ps_face_idx == TRI_FACE_IDX) {
            place_on_face_foreign_proxy_sites(mode = MODE) {
                _test_assert_triangle_proxy_face_child_element_context(0);
                _test_assert_triangle_proxy_edge_child(1);
                _test_assert_triangle_proxy_vertex_child(2);
            }
        }
    }
}

module test_face_local_iterators__parent_coords_preserve_metadata() {
    place_on_faces(_test_punch_poly()) {
        if ($ps_face_idx == STAR_FACE_IDX) {
            place_on_face_filled_boundary_source_edges(mode = MODE, coords = "parent") {
                assert_int_eq($ps_boundary_source_edge_count, 7, "parent-coords source-edge record count");
                assert_int_eq(len($ps_boundary_source_edge_segment2d_local), 2, "parent-coords source edge segment arity");
                assert_int_eq(
                    len($ps_boundary_source_edge_span_segments2d_local),
                    $ps_boundary_source_edge_span_count,
                    "parent-coords source-edge span segment arity"
                );
            }

            place_on_face_boundary_spans(mode = MODE, coords = "parent") {
                assert($ps_boundary_span_count > 0, "parent-coords boundary span count");
                assert_int_eq(len($ps_boundary_span_segment2d_local), 2, "parent-coords boundary span segment arity");
            }
        }
    }
}

module run_TestSelfCrossing() {
    test_ps_face_arrangement__7_3_15_star_has_stable_structure();
    test_ps_face_boundary_model__7_3_15_star_has_true_nonzero_boundary();
    test_ps_face_filled_boundary_source_edges__7_3_15_star_groups_surviving_spans();
    test_ps_face_geom_cut_entries__7_3_15_triangle_records_foreign_cutters();
    test_ps_face_foreign_intrusion_records__7_3_15_triangle_wraps_exact_face_cuts();
    test_ps_face_foreign_intrusion_records__preserves_coincident_foreign_face_provenance();
    test_ps_face_foreign_face_replay_sites__7_3_15_triangle_builds_target_local_frames();
    test_ps_face_foreign_proxy_replay_sites__7_3_15_triangle_includes_edge_and_vertex_candidates();
    test_ps_face_visible_segments__7_3_15_triangle_splits_into_visible_cells();
    test_ps_face_visible_segments__7_3_0_triangle_catches_meeting_cut_edges();
    test_ps_face_filled_boundary_source_edges__7_3_0_triangle_is_simple_boundary();
    test_place_on_face_filled_boundary_source_edges__7_3_15_star_exposes_context();
    test_place_on_face_filled_boundary_source_edges__antitet_uses_span_direction();
    test_place_on_face_foreign_intrusions__7_3_15_triangle_exposes_context();
    test_place_on_face_foreign_face_replay_sites__7_3_15_triangle_exposes_context();
    test_place_on_face_foreign_proxy_sites__7_3_15_triangle_dispatches_face_child();
    test_place_on_face_foreign_proxy_sites__7_3_15_triangle_element_child_uses_source_face_context();
    test_place_on_face_foreign_proxy_sites__7_3_15_triangle_dispatches_edge_and_vertex_children();
    test_face_local_iterators__parent_coords_preserve_metadata();
}

run_TestSelfCrossing();
