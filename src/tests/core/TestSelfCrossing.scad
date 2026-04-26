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

module run_TestSelfCrossing() {
    test_ps_face_arrangement__7_3_15_star_has_stable_structure();
    test_ps_face_boundary_model__7_3_15_star_has_true_nonzero_boundary();
    test_ps_face_filled_boundary_source_edges__7_3_15_star_groups_surviving_spans();
    test_ps_face_geom_cut_entries__7_3_15_triangle_records_foreign_cutters();
    test_ps_face_visible_segments__7_3_15_triangle_splits_into_visible_cells();
    test_ps_face_visible_segments__7_3_0_triangle_catches_meeting_cut_edges();
    test_ps_face_filled_boundary_source_edges__7_3_0_triangle_is_simple_boundary();
    test_place_on_face_filled_boundary_source_edges__7_3_15_star_exposes_context();
    test_place_on_face_filled_boundary_source_edges__antitet_uses_span_direction();
}

run_TestSelfCrossing();
