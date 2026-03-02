use <../../core/prisms.scad>
use <../../core/truncation.scad>
use <../../core/duals.scad>

use <util_demo.scad>

/*
Prism transform playground.

Focus: explore duals + transform combinations on prism families.
*/

sx = 125;
sy = 125;

base_shapes = [
    ["prism(6)", function() poly_prism(6)],
    ["antiprism(6)", function() poly_antiprism(6)],
    ["prism(8, tall)", function() poly_prism(8, height_scale=2.3)],
    ["antiprism(7, twist)", function() poly_antiprism(7, angle=8)]
];

transform_ops = [
    ["identity", function(p) p],
    ["dual", function(p) poly_dual(p)],
    ["truncate", function(p) poly_truncate(p)],
    ["chamfer", function(p) poly_chamfer(p)],
    ["cantellate", function(p) poly_cantellate(p)]
];

// Face families (via external classify probe):
//   fid 0 = cap n-gons
//   fid 1 = side triangles
// Calibrated externally via main_prism_param_calc.scad.
AP6_CANT_ROWS = [
    ["face", "family", 0, ["df", 0.32]],
    ["face", "family", 1, ["df", 0.16]]
];
AP7_CANT_ROWS = [
    ["face", "family", 0, ["df", 0.30]],
    ["face", "family", 1, ["df", 0.10]]
];

PR6_CT_T = 0.06;
PR6_CT_C = 0.03;
AP6_CT_T = 0.06;
AP6_CT_C = 0.03;

for (op = with_index(transform_ops)) {
    for (shape = with_index(base_shapes)) {
        p0 = shape[1][1]();
        p1 = (op[1][0] != "cantellate")
            ? op[1][1](p0)
            : (shape[1][0] == "antiprism(6)")
                ? poly_cantellate(p0, params_overrides=AP6_CANT_ROWS)
            : (shape[1][0] == "antiprism(7, twist)")
                ? poly_cantellate(p0, params_overrides=AP7_CANT_ROWS)
            : poly_cantellate(p0);
        label = str(shape[1][0], " (", op[1][0], ")");
        translate([shape[0] * sx, op[0] * sy, 0]) demo(p1, name=label);
    }
}

// Single focused cantitruncate examples for reference:
// fixed params sampled via main_prism_param_calc.scad
p_pr = poly_prism(6);
p_ap = poly_antiprism(6);

translate([0, len(transform_ops) * sy + 20, 0])
    demo(poly_cantitruncate(p_pr, t=PR6_CT_T, c=PR6_CT_C, params_overrides=undef),
         name="prism(6) cantitruncate fixed");

translate([sx, len(transform_ops) * sy + 20, 0])
    demo(poly_cantitruncate(p_ap, t=AP6_CT_T, c=AP6_CT_C, params_overrides=undef),
         name="antiprism(6) cantitruncate fixed");

translate([2 * sx, len(transform_ops) * sy + 20, 0])
    demo(poly_snub(p_pr, c=0.07, df=0.05, angle=10),
         name="prism(6) snub explicit");

translate([3 * sx, len(transform_ops) * sy + 20, 0])
    demo(poly_snub(p_ap, c=0.04, df=0.008, angle=15),
         name="antiprism(6) snub explicit");
