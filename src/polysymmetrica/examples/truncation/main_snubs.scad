use <../../core/truncation.scad>
use <../../core/duals.scad>
use <../../core/render.scad>
use <../../models/platonics_all.scad>
use <../truncation/util_demo.scad>

/**
Quick snub demo (snub cube + snub dodecahedron).
*/

spacing = 100;

//// Basic Snubs
//p = poly_snub(hexahedron());
////p = poly_snub(dodecahedron());
//demo(p, name="snub");
//poly_describe(p, detail = 3);


//// Regular default snubs
//translate([-spacing, 0, 0]) demo(poly_snub(tetrahedron()), name="snub_tet");
//translate([0, 0, 0]) demo(poly_snub(hexahedron()), name="snub");
//translate([spacing, 0, 0]) demo(poly_snub(dodecahedron()), name="snub_dodecahedron");

//// Cuboctahedron (rectified octahedron): compare family preference in default solve.
//b1 = poly_rectify(octahedron());
//translate([spacing*2, 0, 0]) demo(poly_snub(b1), name="snub_cubocta (auto)");
//translate([spacing*3, 0, 0]) demo(poly_snub(b1, family_k=3), name="snub_cubocta (family 3)");
//translate([spacing*4, 0, 0]) demo(poly_snub(b1, family_k=4), name="snub_cubocta (family 4)");
//
//b2 = poly_dual(b1);
//translate([spacing*5, 0, 0]) demo(poly_snub(b2), name="snub_dual(cubocta)");




//b3 = poly_dual(poly_rectify(icosahedron()));
//translate([spacing*5, 0, 0]) demo(poly_snub(b3, angle=20, c=0.01, df = 0.01), name="snub_rh_triac");


//// params_overrides examples
//// Format:
//// [
////   ["face", <face_family_id>, ["angle", deg], ["df", value]],
////   ["vert", <vert_family_id>, ["c", value]]   // or ["de", value]
//// ]
////
//// Family ids come from poly_classify(poly, detail=1). For cuboctahedron there
//// are two face families (triangles/squares) and typically one vertex family.
//
//p0 = poly_rectify(octahedron()); // cuboctahedron
//seed = _ps_snub_default_params(p0, family_k=3);
//seed_df = seed[0];
//seed_angle = seed[1];
//seed_c = seed[2];
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
//        angle=seed_angle,
//        c=seed_c,
//        df=seed_df,
//        params_overrides=[
//            ["face", "family", 0, ["df", seed_df * 0.9], ["angle", seed_angle + 2]],
//            ["face", "family", 1, ["df", seed_df * 1.1], ["angle", seed_angle - 2]],
//            ["vert", "family", 0, ["c", seed_c]]
//        ]
//    ),
//    name="snub cubocta (solver-seeded overrides)"
//);

// Experiment: snub a rhombicosidodecahedron
b4 = poly_cantellate(icosahedron());
translate([0, 0, 0]) poly_render(b4,40);
translate([spacing, 0, 0]) demo(poly_snub(b4, angle=20, c=0.01, df = 0.01), ir = 40, name="canttrunc_cuboctahedron");
