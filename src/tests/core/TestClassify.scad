use <../testing_util.scad>
use <../../polysymmetrica/core/classify.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/archimedians_all.scad>
use <../../polysymmetrica/models/catalans_all.scad>

function _classify_counts(poly, detail=0) =
    let(cls = poly_classify(poly, detail))
    [len(cls[0]), len(cls[1]), len(cls[2])];

module test_classify__platonics_single_family() {
    plats = [
        tetrahedron(),
        hexahedron(),
        octahedron(),
        dodecahedron(),
        icosahedron()
    ];

    for (p = plats) {
        counts = _classify_counts(p, 0);
        assert_int_eq(counts[0], 1, "platonics: face families = 1");
        assert_int_eq(counts[1], 1, "platonics: edge families = 1");
        assert_int_eq(counts[2], 1, "platonics: vertex families = 1");
    }
}

module test_classify__truncated_octa_families() {
    p = poly_truncate(octahedron());
    counts = _classify_counts(p, 0);
    assert_int_eq(counts[0], 2, "trunc octa: face families = 2");
    assert_int_eq(counts[1], 1, "trunc octa: edge families = 1");
    assert_int_eq(counts[2], 1, "trunc octa: vertex families = 1");
}

module test_classify__rhombi_families() {
    r = rhombicuboctahedron();
    g = great_rhombicuboctahedron();
    counts_r = _classify_counts(r, 0);
    counts_g = _classify_counts(g, 0);

    assert_int_eq(counts_r[0], 2, "rhombicubocta: face families = 2");
    assert_int_eq(counts_r[1], 1, "rhombicubocta: edge families = 1");
    assert_int_eq(counts_r[2], 1, "rhombicubocta: vertex families = 1");

    assert_int_eq(counts_g[0], 3, "great rhombicubocta: face families = 3");
    assert_int_eq(counts_g[2], 1, "great rhombicubocta: vertex families = 1");
}

module test_classify__rhombi_duals_face_family() {
    d1 = deltoidal_icositetrahedron();
    d2 = disdyakis_dodecahedron();
    counts1 = _classify_counts(d1, 0);
    counts2 = _classify_counts(d2, 0);

    assert_int_eq(counts1[0], 1, "deltoidal icositetra: face families = 1");
    assert_int_eq(counts2[0], 1, "disdyakis dodeca: face families = 1");
}


module run_TestClassify() {
    test_classify__platonics_single_family();
    test_classify__truncated_octa_families();
    test_classify__rhombi_families();
    test_classify__rhombi_duals_face_family();
}

