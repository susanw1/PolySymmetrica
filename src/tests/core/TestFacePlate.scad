use <../../polysymmetrica/examples/printing/face_plate.scad>

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module test_face_plate_requires_segmented_loops__simple_face_zero_inset_is_false() {
    pts = [[-10, -10], [10, -10], [10, 10], [-10, 10]];
    assert(!_face_plate_requires_segmented_loops(pts, 0, 1e-8), "simple face with zero inset should not force segmentation");
}

module test_face_plate_requires_segmented_loops__star_face_zero_inset_is_true() {
    pts = [for (i = [0:4]) [10 * cos(-144 * i), 10 * sin(-144 * i)]];
    assert(_face_plate_requires_segmented_loops(pts, 0, 1e-8), "self-intersecting star face with zero inset should force segmentation");
}

module test_face_plate_requires_segmented_loops__simple_face_positive_inset_is_true() {
    pts = [[-10, -10], [10, -10], [10, 10], [-10, 10]];
    assert(_face_plate_requires_segmented_loops(pts, 0.01, 1e-8), "positive inset should use segmented loop path");
}

module test_face_plate_requires_segmented_loops__star_face_segments_even_with_zero_inset() {
    pts = [for (i = [0:4]) [10 * cos(-144 * i), 10 * sin(-144 * i)]];
    diheds = [60, 80, 100, 120, 140];
    insets = [for (_ = diheds) 0];
    body_loops = _face_plate_requires_segmented_loops(pts, 0, 1e-8)
        ? _segmented_body_loops(pts, diheds, insets, 1e-8)
        : [[pts, diheds]];
    assert_int_eq(len(body_loops), 6, "star face should still split into simple body loops at zero inset");
}

module run_TestFacePlate() {
    test_face_plate_requires_segmented_loops__simple_face_zero_inset_is_false();
    test_face_plate_requires_segmented_loops__star_face_zero_inset_is_true();
    test_face_plate_requires_segmented_loops__simple_face_positive_inset_is_true();
    test_face_plate_requires_segmented_loops__star_face_segments_even_with_zero_inset();
}

run_TestFacePlate();
