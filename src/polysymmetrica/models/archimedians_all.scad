use <../core/truncation.scad>
use <../core/cantellation.scad>
use <../models/platonics_all.scad>

/**
Archimedean solids derived from Platonic bases.
Snubs are left as stubs for now.
*/

// ---- Truncations ----

function truncated_tetrahedron() = poly_truncate(tetrahedron());
function truncated_cube() = poly_truncate(hexahedron());
function truncated_octahedron() = poly_truncate(octahedron());
function truncated_dodecahedron() = poly_truncate(dodecahedron());
function truncated_icosahedron() = poly_truncate(icosahedron());

// ---- Rectifications ----

function cuboctahedron() = poly_rectify(hexahedron());
function icosidodecahedron() = poly_rectify(dodecahedron());

// ---- Cantellations (small rhombi*) ----

function rhombicuboctahedron() =
    let(
        base = hexahedron(),
        df = cantellate_square_df(base, 0, 1, 40, 0)
    )
    poly_cantellate(base, df);

function rhombicosidodecahedron() =
    let(
        base = dodecahedron(),
        df = cantellate_square_df(base, 0, 1, 40, 0)
    )
    poly_cantellate(base, df);

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

/**
Return all Archimedean solids as [[name, poly], ...].
Note: this constructs every poly when called; snubs are undef.
*/
function archimedians_all() = [
    ["truncated_tetrahedron", truncated_tetrahedron()],
    ["truncated_cube", truncated_cube()],
    ["truncated_octahedron", truncated_octahedron()],
    ["truncated_dodecahedron", truncated_dodecahedron()],
    ["truncated_icosahedron", truncated_icosahedron()],
    ["cuboctahedron", cuboctahedron()],
    ["icosidodecahedron", icosidodecahedron()],
    ["rhombicuboctahedron", rhombicuboctahedron()],
    ["rhombicosidodecahedron", rhombicosidodecahedron()],
    ["great_rhombicuboctahedron", great_rhombicuboctahedron()],
    ["great_rhombicosidodecahedron", great_rhombicosidodecahedron()],
    ["snub_cube", snub_cube()],
    ["snub_dodecahedron", snub_dodecahedron()]
];
