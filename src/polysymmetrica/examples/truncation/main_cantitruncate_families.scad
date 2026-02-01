use <../../core/truncation.scad>
use <../../core/render.scad>
use <../../models/regular_all.scad>
use <util_demo.scad>

/**
Cantitruncate with per-face-family c values (dominant family).
Example on cuboctahedron (rectified octahedron).
*/

base = poly_rectify(octahedron());

// Dominant squares (size 4), include per-edge-family c to improve planarity
sol_sq = solve_cantitruncate_dominant_edges(base, 4);
p_sq = poly_cantitruncate_families(base, sol_sq[0], sol_sq[1], c_edge_by_pair=sol_sq[2]);

// Dominant triangles (size 3)
sol_tri = solve_cantitruncate_dominant_edges(base, 3);
p_tri = poly_cantitruncate_families(base, sol_tri[0], sol_tri[1], c_edge_by_pair=sol_tri[2]);

translate([-100, 0, 0]) demo(base);
translate([0, 0, 0]) demo(p_sq);
translate([100, 0, 0]) demo(p_tri);
