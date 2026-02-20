use <../../core/truncation.scad>
use <../../core/params.scad>
use <../../core/duals.scad>
use <../../core/render.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

/**
Quick snub demo (snub cube + snub dodecahedron).
*/

spacing = 100;

// Basic Snubs
p = poly_snub(hexahedron(), df = 0.05, angle = 15, c = undef, de = 0.1, handedness = 1);
//p = poly_snub(dodecahedron());
demo(p, name="snub");
//poly_describe(p, detail = 3);


//for (h = [1, -1]) {
//    // Regular default snubs
//    translate([-spacing, h*spacing, 0]) demo(poly_snub(tetrahedron(), handedness = h), name="snub_tet");
//    translate([0, h*spacing, 0]) demo(poly_snub(hexahedron(), handedness = h), name="snub_cube");
//    translate([spacing, h*spacing, 0]) demo(poly_snub(dodecahedron(), handedness = h), name="snub_dodecahedron");
//
//    // Cuboctahedron (rectified octahedron): baseline and overrides.
//    b1 = poly_rectify(octahedron());
//    translate([spacing*2, h*spacing, 0]) demo(poly_snub(b1, handedness = h), name="snub_cubocta (auto)");
//
//    b2 = poly_dual(b1);
//    translate([spacing*3, h*spacing, 0]) demo(poly_snub(b2, handedness = h), name="snub_rh_dod");
//
//    b3 = (poly_rectify(icosahedron()));
//    translate([spacing*4, h*spacing, 0]) demo(poly_snub(b3, handedness = h), name="snub_icosidod");
//
//    b4 = poly_dual(b3);
//    translate([spacing*5, h*spacing, 0]) demo(poly_snub(b4, handedness = h), name="snub_rh_triaconta");
//}

//// params_overrides examples
//// Format:
//// [
////   ["face", "family", <face_family_id>, ["angle", deg], ["df", value]],
////   ["vert", "family", <vert_family_id>, ["c", value]]   // or ["de", value]
//// ]
////
//// Family ids come from poly_classify(poly, detail=1). For cuboctahedron there
//// are two face families (triangles/squares) and typically one vertex family.
//
//p0 = poly_rectify(octahedron()); // cuboctahedron
//seed_rows = ps_snub_default_params_overrides(p0);
//
////// Baseline (fully automatic defaults).
//translate([0, 0, 0]) demo(
//    poly_snub(p0),
//    name="snub cubocta (baseline auto)"
//);
////
////// Exploratory override: tweak family 0 only (manual values).
//translate([spacing, 0, 0]) demo(
//    poly_snub(
//        p0,
//        angle=15,
//        c=0.04,
//        df=0.035,
//        params_overrides=[
//            ["face", "family", 0, ["angle", 20], ["df", 0.03]]
//        ]
//    ),
//    name="snub cubocta (explore face#0)"
//);
////
////// Exploratory override: tweak family 1 only (manual values).
//translate([spacing*2, 0, 0]) demo(
//    poly_snub(
//        p0,
//        angle=15,
//        c=0.04,
//        df=0.035,
//        params_overrides=[
//            ["face", "family", 1, ["angle", 10], ["df", 0.045]]
//        ]
//    ),
//    name="snub cubocta (explore face#1)"
//);
////
////// Solver-seeded family override: starts from computed defaults.
////// This avoids hard-coding "magic" global constants.
//translate([spacing*3, 0, 0]) demo(
//    poly_snub(
//        p0,
//        params_overrides=concat(
//            seed_rows,
//            [
//                ["face", "family", 0, ["df", ps_params_get(seed_rows, "face", "df", 0, 0) * 0.9], ["angle", ps_params_get(seed_rows, "face", "angle", 0, 0) + 2]],
//                ["face", "family", 1, ["df", ps_params_get(seed_rows, "face", "df", 0, 1) * 1.1], ["angle", ps_params_get(seed_rows, "face", "angle", 0, 1) - 2]]
//            ]
//        )
//    ),
//    name="snub cubocta (solver-seeded overrides)"
//);

//// Experiment: snub a rhombicosidodecahedron
//b4 = poly_cantellate(icosahedron());
//translate([0, 0, 0]) poly_render(b4,40);
//translate([spacing, 0, 0]) demo(poly_snub(b4, angle=20, c=0.01, df = 0.01), ir = 40, name="snub_rhombicosidodecahedron");
