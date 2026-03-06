use <../../polysymmetrica/core/placement.scad>
use <../../polysymmetrica/core/classify.scad>
use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/archimedians_all.scad>

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
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
            place_on_face_segments() {
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

module run_TestPlacement() {
    test_place_on_faces__family_ids_and_counts_from_classify();
    test_place_on_edges__family_ids_and_counts_from_classify();
    test_place_on_vertices__family_ids_and_counts_from_classify();
    test_place_on_faces__auto_classify_matches_precomputed();
    test_place_on_all__cube_single_family();
    test_place_on_edges__no_auto_classify_by_default();
    test_place_on_face_segments__star_face_split();
}
