use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/placement.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/core/segments.scad>

EPS = 1e-8;
MODE = "nonzero";
STAR_FACE_IDX = 1;
TRI_FACE_IDX = 12;

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_list_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

function _test_punch_poly() =
    poly_antiprism(7, 3, angle = 15);

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

module run_TestSelfCrossing() {
    test_ps_face_arrangement__7_3_15_star_has_stable_structure();
    test_ps_face_boundary_model__7_3_15_star_has_true_nonzero_boundary();
    test_ps_face_geom_cut_entries__7_3_15_triangle_records_foreign_cutters();
    test_ps_face_visible_segments__7_3_15_triangle_splits_into_visible_cells();
}

run_TestSelfCrossing();
