use <../core/truncation.scad>
use <../models/regular_all.scad>

// Great rhombi* via cantitruncate uniform solver (fast, edge-face metrics only)
function great_rhombicuboctahedron() =
    let(
        base = hexahedron(),
        sol = solve_cantitruncate_trig(base)
    )
    poly_cantitruncate(base, sol[0], sol[1]);

function great_rhombicosidodecahedron() =
    let(
        base = dodecahedron(),
        sol = solve_cantitruncate_trig(base)
    )
    poly_cantitruncate(base, sol[0], sol[1]);
