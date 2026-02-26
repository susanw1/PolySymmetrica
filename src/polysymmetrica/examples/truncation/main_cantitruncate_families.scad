use <../../core/truncation.scad>
use <../../core/solvers.scad>
use <../../core/render.scad>
use <../../models/platonics_all.scad>
use <util_demo.scad>

/**
Cantitruncate with per-face-family c values (dominant family).
Example on cuboctahedron (rectified octahedron).
*/

base = poly_rectify(octahedron());

// Dominant squares (size 4), include per-edge-family c to improve planarity.
params_sq = solve_cantitruncate_dominant_edges_params(base, 4);
p_sq = poly_cantitruncate(base, t=0, c=0, params_overrides=params_sq);

// Dominant triangles (size 3)
params_tri = solve_cantitruncate_dominant_edges_params(base, 3);
p_tri = poly_cantitruncate(base, t=0, c=0, params_overrides=params_tri);

translate([-100, 0, 0]) demo(base);
translate([0, 0, 0]) demo(p_sq);
translate([100, 0, 0]) demo(p_tri);
