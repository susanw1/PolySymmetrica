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

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

function area2(poly) =
    let(n = len(poly))
    (n < 3) ? 0 :
    sum([
        for (i = [0:1:n-1])
            let(a = poly[i], b = poly[(i + 1) % n])
            a[0] * b[1] - b[0] * a[1]
    ]);

function _poly_same_pts(poly_a, poly_b, eps=1e-6) =
    (len(poly_a) == len(poly_b)) &&
    len([
        for (i = [0:1:len(poly_a)-1])
            if (ps_point_eq(poly_a[i], poly_b[i], eps)) 1
    ]) == len(poly_a);

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

module test_ps_fr_loop_from_halfplanes__clips_square() {
    probe = [0, 0];
    lines = [
        [[1, 0], -1, true],
        [[0, 1], -1, true],
        [[-1, 0], -1, true],
        [[0, -1], -1, true]
    ];
    poly = _ps_fr_loop_from_halfplanes(lines, probe);
    assert(len(poly) == 4, "clipped square arity");
    assert_near(abs(area2(poly)), 8, 1e-6, "clipped square area2");
}

module test_ps_face_visible_cell_loop_at_z_clipped__shrinks_bad_antiprism_cell() {
    p = poly_antiprism(n=7, p=3, angle=15);
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cell = vis[0];
            loop_adj = _ps_fr_visible_cell_loop_at_z(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, 1.6, -0.4, 1.6, 1.1, 1e-8);
            loop_clip = ps_face_visible_cell_loop_at_z_clipped(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, 1.6, -0.4, 1.6, 1.1, 1e-8);

            assert(len(loop_clip) >= 3, "clipped loop should remain polygonal");
            assert(abs(area2(loop_clip)) < abs(area2(loop_adj)), "clipped loop should be tighter than adjacent-line loop for bad antiprism cell");
        }
    }
}

module test_ps_face_visible_cell_region_planes__include_slab_and_edges() {
    p = poly_antiprism(n=7, p=3, angle=15);
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cell = vis[0];
            planes = ps_face_visible_cell_region_planes(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, -0.4, 1.6, -0.4, 1.6, 1.1, 1e-8);
            assert(len(planes) == len(cell[0]) + 2, "region planes should include one plane per edge plus z slab");
        }
    }
}

module test_ps_face_visible_cell_region_planes__can_append_run_end_planes() {
    p = poly_antiprism(n=7, p=3, angle=15);
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cell = vis[0];
            base_planes = ps_face_visible_cell_region_planes(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, -0.4, 1.6, -0.4, 1.6, 1.1, 1e-8);
            run_end_entries = ps_face_visible_cell_cut_run_end_entries(cell);
            planes_with_runs = ps_face_visible_cell_region_planes(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, -0.4, 1.6, -0.4, 1.6, 1.1, 1e-8, true);

            assert_int_eq(len(planes_with_runs), len(base_planes) + len(run_end_entries), "run-end planes should append to the base region-plane set");
        }
    }
}

module test_ps_face_visible_cell_cut_run_end_entries__single_span_orientation() {
    cell = [
        [[0, 0], [2, 0], [2, 1], [0, 1]],
        undef,
        [0, 1, 2, 3],
        ["parent", "cut", "parent", "parent"],
        [undef, 7, undef, undef],
        [undef, 0, undef, undef]
    ];
    ends = ps_face_visible_cell_cut_run_end_entries(cell);
    assert(len(ends) == 2, "single cut span should have start and end planes");

    start = ends[0];
    finish = ends[1];

    assert_near(start[0][0][0], 0, 1e-9, "start normal x");
    assert_near(start[0][0][1], 1, 1e-9, "start normal y");
    assert_near(start[0][1], 0, 1e-9, "start plane should pass through run start");
    assert(start[0][2], "start plane keep side should be >= along edge direction");
    assert(start[1] == 1 && start[2] == true && start[3] == 0, "start entry metadata");

    assert_near(finish[0][0][0], 0, 1e-9, "end normal x");
    assert_near(finish[0][0][1], 1, 1e-9, "end normal y");
    assert_near(finish[0][1], 1, 1e-9, "end plane should pass through run end");
    assert(!finish[0][2], "end plane keep side should be <= along edge direction");
    assert(finish[1] == 1 && finish[2] == false && finish[3] == 0, "end entry metadata");
}

module test_ps_face_visible_cell_cut_run_end_entries__split_spans_stay_distinct() {
    p = poly_antiprism(n=7, p=3, angle=15);
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cell = vis[0];
            ends = ps_face_visible_cell_cut_run_end_entries(cell);
            run_ids = cell[5];
            cut_runs = [
                for (rid = run_ids)
                    if (!is_undef(rid))
                        rid
            ];
            distinct_runs = [
                for (ri = [0:1:len(cut_runs)-1])
                    if (len([for (rj = [0:1:ri-1]) if (cut_runs[rj] == cut_runs[ri]) 1]) == 0)
                        cut_runs[ri]
            ];

            assert_int_eq(len(ends), 2 * len(distinct_runs), "each split cut span should contribute two run-end planes");
            for (rid = distinct_runs) {
                rid_ends = [for (e = ends) if (e[3] == rid) e];
                assert_int_eq(len(rid_ends), 2, str("run should keep distinct start/end planes rid=", rid));
            }
        }
    }
}

module test_ps_face_visible_cell_loop_at_z_from_region_planes__matches_clipped_loop() {
    p = poly_antiprism(n=7, p=3, angle=15);
    zs = [-0.4, 0.4, 1.6];
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cell = vis[0];

            for (z = zs) {
                loop_clip = ps_face_visible_cell_loop_at_z_clipped(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, z, -0.4, 1.6, 1.1, 1e-8);
                loop_planes = ps_face_visible_cell_loop_at_z_from_region_planes(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, z, -0.4, 1.6, -0.4, 1.6, 1.1, 1e-8);

                assert(len(loop_planes) == len(loop_clip), str("plane-sliced loop arity should match clipped loop at z=", z));
                assert(_poly_same_pts(loop_planes, loop_clip, 1e-5), str("plane-sliced loop should match clipped loop at z=", z));
            }
        }
    }
}

module test_ps_fr_atom_can_use_clipped_loops__rejects_inner_edge_atoms() {
    p = poly_antiprism(n=7, p=3, angle=15);
    levels = [-0.4, 1.6, 1.8];
    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            vis = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            atoms = _ps_fr_visible_cell_atoms(vis[0], 1e-8);

            eligible = [
                for (atom = atoms)
                    let(
                        clip_loops = [
                            for (z = levels)
                                [ps_face_visible_cell_loop_at_z_clipped(atom, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, z, -0.4, 1.6, 1.1, 1e-8), z]
                        ],
                        clip_sizes = [for (lz = clip_loops) len(lz[0])],
                        clip_stable = min(clip_sizes) >= 3 &&
                            len([for (n = clip_sizes) if (n != clip_sizes[0]) 1]) == 0
                    )
                    if (_ps_fr_atom_has_inner_edges(atom) && clip_stable)
                        _ps_fr_atom_can_use_clipped_loops(atom, clip_loops, 1e-8)
            ];

            assert(len(eligible) > 0, "fixture should include inner-edge atoms with superficially stable clipped-loop arity");
            assert_int_eq(len([for (ok = eligible) if (ok) 1]), 0, "inner-edge atoms must not take clipped-loop lofting");
        }
    }
}

module test_ps_face_cut_profile2d_from_cutter_normal__matches_f2_f5_join() {
    p = poly_antiprism(n=7, p=3, angle=15);

    place_on_faces(p) {
        if ($ps_face_idx == 2) {
            cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            place_on_face_visible_segments("nonzero", 1e-8, true) {
                if ($ps_vis_seg_idx == 0) {
                    cell = [$ps_vis_seg_pts2d, undef, $ps_vis_seg_edge_ids, $ps_vis_seg_edge_kinds, $ps_vis_seg_cut_entry_ids];
                    probe = _ps_seg_cycle_probe_point(cell[0], 1e-8);
                    inward_n = _ps_fr_cell_edge_inward_n(cell[0], 1, probe, 1e-8);
                    cid = $ps_vis_seg_cut_entry_ids[1];
                    cutter_face = cut_entries[cid][1];
                    cut_dihed = cut_entries[cid][2];
                    cutter_n3 = ps_face_frame_normal($ps_poly_verts_local, $ps_poly_faces_idx[cutter_face]);
                    prof = ps_face_cut_profile2d_from_cutter_normal(-0.4, 1.2, inward_n, cutter_n3, 1.1, cut_dihed, 1e-8);

                    assert_near(cut_dihed, 94.23407613227815, 1e-4, "f2/c0/e1 cut dihedral");
                    assert_near(ps_face_cut_join_dihed(cut_dihed), 265.7659238677219, 1e-4, "f2/c0/e1 join dihedral");
                    assert_near(prof[0][0], 0.55, 1e-5, "f2/c0/e1 lower u");
                    assert_near(prof[0][1], -0.4, 1e-9, "f2/c0/e1 lower z");
                    assert_near(prof[1][0], 2.0359176277869236, 1e-5, "f2/c0/e1 upper u");
                    assert_near(prof[1][1], 1.2, 1e-9, "f2/c0/e1 upper z");
                }
            }
        }
    }

    place_on_faces(p) {
        if ($ps_face_idx == 5) {
            cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, 1e-8, "nonzero", true);
            place_on_face_visible_segments("nonzero", 1e-8, true) {
                if ($ps_vis_seg_idx == 1) {
                    cell = [$ps_vis_seg_pts2d, undef, $ps_vis_seg_edge_ids, $ps_vis_seg_edge_kinds, $ps_vis_seg_cut_entry_ids];
                    probe = _ps_seg_cycle_probe_point(cell[0], 1e-8);
                    inward_n = _ps_fr_cell_edge_inward_n(cell[0], 1, probe, 1e-8);
                    cid = $ps_vis_seg_cut_entry_ids[1];
                    cutter_face = cut_entries[cid][1];
                    cut_dihed = cut_entries[cid][2];
                    cutter_n3 = ps_face_frame_normal($ps_poly_verts_local, $ps_poly_faces_idx[cutter_face]);
                    prof = ps_face_cut_profile2d_from_cutter_normal(-0.4, 1.2, inward_n, cutter_n3, 1.1, cut_dihed, 1e-8);

                    assert_near(cut_dihed, 94.23407613227815, 1e-4, "f5/c1/e1 cut dihedral");
                    assert_near(ps_face_cut_join_dihed(cut_dihed), 265.7659238677219, 1e-4, "f5/c1/e1 join dihedral");
                    assert_near(prof[0][0], 0.55, 1e-5, "f5/c1/e1 lower u");
                    assert_near(prof[0][1], -0.4, 1e-9, "f5/c1/e1 lower z");
                    assert_near(prof[1][0], 2.0359176277869236, 1e-5, "f5/c1/e1 upper u");
                    assert_near(prof[1][1], 1.2, 1e-9, "f5/c1/e1 upper z");
                }
            }
        }
    }
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
    test_ps_fr_loop_from_halfplanes__clips_square();
    test_ps_face_visible_cell_loop_at_z_clipped__shrinks_bad_antiprism_cell();
    test_ps_face_visible_cell_region_planes__include_slab_and_edges();
    test_ps_face_visible_cell_region_planes__can_append_run_end_planes();
    test_ps_face_visible_cell_cut_run_end_entries__single_span_orientation();
    test_ps_face_visible_cell_cut_run_end_entries__split_spans_stay_distinct();
    test_ps_face_visible_cell_loop_at_z_from_region_planes__matches_clipped_loop();
    test_ps_fr_atom_can_use_clipped_loops__rejects_inner_edge_atoms();
    test_ps_face_cut_profile2d_from_cutter_normal__matches_f2_f5_join();
    test_ps_face_region_inset_at_z__zero_at_face_plane();
    test_ps_clip_to_face_region_ctx__smoke();
    test_ps_clip_to_visible_face_segments_ctx__smoke();
    test_ps_clip_to_visible_face_segments_ctx__cut_bands_smoke();
}

run_TestFaceRegions();
