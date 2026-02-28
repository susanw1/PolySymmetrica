use <../../core/prisms.scad>
use <../../core/truncation.scad>
use <../../core/duals.scad>
use <../../core/solvers.scad>
use <../../core/params.scad>

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

for (op = with_index(transform_ops)) {
    for (shape = with_index(base_shapes)) {
        p0 = shape[1][1]();
        p1 = op[1][1](p0);
        label = str(shape[1][0], " (", op[1][0], ")");
        translate([shape[0] * sx, op[0] * sy, 0]) demo(p1, name=label);
    }
}

// Single focused cantitruncate examples for reference:
// - regular-ish dominant-family solve on prism(6)
// - same on antiprism(6)
p_pr = poly_prism(6);
p_ap = poly_antiprism(6);
rows_pr = solve_cantitruncate_dominant_edges_params(p_pr, 6);
rows_ap = solve_cantitruncate_dominant_edges_params(p_ap, 6);

translate([0, len(transform_ops) * sy + 20, 0])
    demo(poly_cantitruncate(p_pr, t=0, c=0, params_overrides=rows_pr),
         name="prism(6) cantitruncate dominant-size=6");

translate([sx, len(transform_ops) * sy + 20, 0])
    demo(poly_cantitruncate(p_ap, t=0, c=0, params_overrides=rows_ap),
         name="antiprism(6) cantitruncate dominant-size=6");

translate([2 * sx, len(transform_ops) * sy + 20, 0])
    demo(poly_snub(p_pr, c=0.07, df=0.10, angle=15),
         name="prism(6) snub explicit");

translate([3 * sx, len(transform_ops) * sy + 20, 0])
    demo(poly_snub(p_ap, c=0.04, df=0.008, angle=15),
         name="antiprism(6) snub explicit");
