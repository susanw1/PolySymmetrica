use <../core/truncation.scad>
use <../models/platonics_all.scad>

/**
Archimedean solids derived from Platonic bases.
Snubs are left as stubs for now.
*/

// ---- Truncations (t = 1/3) ----

function truncated_tetrahedron() = poly_truncate(tetrahedron());
function truncated_cube() = poly_truncate(hexahedron());
function truncated_octahedron() = poly_truncate(octahedron());
function truncated_dodecahedron() = poly_truncate(dodecahedron());
function truncated_icosahedron() = poly_truncate(icosahedron());

// ---- Rectifications ----

function cuboctahedron() = poly_rectify(hexahedron());
function icosidodecahedron() = poly_rectify(dodecahedron());

// ---- Cantellations (small rhombi*) ----

function rhombicuboctahedron() = poly_cantellate(hexahedron(), 1/3);
function rhombicosidodecahedron() = poly_cantellate(dodecahedron(), 1/3);

// ---- Cantitruncations (great rhombi*) ----

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

// ---- Snubs (stubs) ----

function snub_cube() = undef; // TODO
function snub_dodecahedron() = undef; // TODO
