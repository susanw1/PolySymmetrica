use <../core/truncation.scad>
use <../models/regular_all.scad>

// Great rhombi* via cantitruncate uniform solver (fast, edge-face metrics only)
function great_rhombicuboctahedron(steps=6, rounds=3, w_edge=1, w_square=1) =
    let(
        base = hexahedron(),
        sol = solve_cantitruncate_uniform_fast_refine(base, 0, 1, 0, 1, steps, rounds, w_edge, w_square)
    )
    let(_ = echo("great_rhombicuboctahedron sol", sol))
    poly_cantitruncate(base, sol[0], sol[1]);

function great_rhombicosidodecahedron(steps=6, rounds=3, w_edge=1, w_square=1) =
    let(
        base = dodecahedron(),
        sol = solve_cantitruncate_uniform_fast_refine(base, 0, 1, 0, 1, steps, rounds, w_edge, w_square)
    )
    let(_ = echo("great_rhombicosidodecahedron sol", sol))
    poly_cantitruncate(base, sol[0], sol[1]);
