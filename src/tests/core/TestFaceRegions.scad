use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/face_regions.scad>
use <../../polysymmetrica/core/placement.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/models/platonics_all.scad>

module assert_near(a, b, eps=1e-9, msg="") {
    assert(abs(a - b) <= eps, str(msg, " expected=", b, " got=", a));
}

module assert_le(a, b, msg="") {
    assert(a <= b, str(msg, " expected ", a, " <= ", b));
}

module assert_ge(a, b, msg="") {
    assert(a >= b, str(msg, " expected ", a, " >= ", b));
}

module test_ps_face_cut_join_dihed__is_complement() {
    assert_near(ps_face_cut_join_dihed(102.857142857), 257.142857143, 1e-6, "join dihedral should complement cutter dihedral");
    assert_near(ps_face_cut_join_dihed(120), 240, 1e-9, "join dihedral complement");
}

module test_ps_face_cut_profile2d_from_cutter_normal__is_planar_across_span() {
    prof = ps_face_cut_profile2d_from_cutter_normal(
        -1,
        2,
        [1, 0],
        [-0.8, 0, 0.6],
        0.2
        ,
        126.86989764584402
    );
    assert(len(prof) == 2, "normal-derived cut profile should define one straight cut plane across the active z span");
    assert_near(prof[0][0], 0.1, 1e-9, "profile should stay on the kept side over the full z span");
    assert_near(prof[0][1], -1, 1e-9, "lower sample z");
    assert_near(prof[1][0], 1.6, 1e-9, "profile should use cut dihedral for slope magnitude while staying on the kept side");
    assert_near(prof[1][1], 2, 1e-9, "upper sample z");
}

module test_ps_face_visible_cell_mask_loop__cut_edges_pull_inward() {
    cell = [
        [[0, 0], [4, 0], [4, 2], [0, 2]],
        undef,
        undef,
        ["parent", "cut", "parent", "cut"],
        [undef, 0, undef, 1]
    ];
    loop = ps_face_visible_cell_mask_loop(cell, 0.2);
    assert(len(loop) == 4, "mask loop arity");
    assert(loop[1][0] < 4, "right cut edge should pull inward");
    assert(loop[3][0] > 0, "left cut edge should pull inward");
}

module test_ps_face_region_inset_at_z__zero_at_face_plane() {
    assert_near(ps_face_region_inset_at_z(120, 0), 0, 1e-9, "face-region inset at z=0");
}

module test_ps_clip_to_face_region_ctx__smoke() {
    place_on_faces(hexahedron()) {
        if ($ps_face_idx == 0)
            ps_clip_to_face_region_ctx(-1, 1)
                cube([4, 4, 4], center = true);
    }
}

module test_ps_clip_to_visible_face_segments_ctx__smoke() {
    place_on_faces(poly_prism(5, 2)) {
        if ($ps_vertex_count == 4)
            ps_clip_to_visible_face_segments_ctx(-1, 1, cut_clearance = 0.1, mode = "nonzero")
                cube([6, 6, 6], center = true);
    }
}

module test_ps_clip_to_visible_face_segments_ctx__cut_bands_smoke() {
    place_on_faces(poly_prism(5, 2)) {
        if ($ps_vertex_count == 4)
            ps_clip_to_visible_face_segments_ctx(-1, 1, cut_clearance = 0.1, mode = "nonzero", apply_cut_bands = true)
                cube([6, 6, 6], center = true);
    }
}

module run_TestFaceRegions() {
    test_ps_face_cut_join_dihed__is_complement();
    test_ps_face_cut_profile2d_from_cutter_normal__is_planar_across_span();
    test_ps_face_visible_cell_mask_loop__cut_edges_pull_inward();
    test_ps_face_region_inset_at_z__zero_at_face_plane();
    test_ps_clip_to_face_region_ctx__smoke();
    test_ps_clip_to_visible_face_segments_ctx__smoke();
    test_ps_clip_to_visible_face_segments_ctx__cut_bands_smoke();
}

run_TestFaceRegions();
