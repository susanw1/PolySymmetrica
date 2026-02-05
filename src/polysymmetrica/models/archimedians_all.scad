use <../core/truncation.scad>
use <../core/solvers.scad>
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

function snub_cube() = poly_snub(hexahedron());
function snub_dodecahedron() = poly_snub(dodecahedron());

/**
Return all Archimedean solids as [[name, fn], ...],
where fn is a zero-arg function returning the poly.
Snubs return undef.
*/
function archimedians_all() = [
    ["truncated_tetrahedron", function() truncated_tetrahedron()],
    ["truncated_cube", function() truncated_cube()],
    ["truncated_octahedron", function() truncated_octahedron()],
    ["truncated_dodecahedron", function() truncated_dodecahedron()],
    ["truncated_icosahedron", function() truncated_icosahedron()],
    ["cuboctahedron", function() cuboctahedron()],
    ["icosidodecahedron", function() icosidodecahedron()],
    ["rhombicuboctahedron", function() rhombicuboctahedron()],
    ["rhombicosidodecahedron", function() rhombicosidodecahedron()],
    ["great_rhombicuboctahedron", function() great_rhombicuboctahedron()],
    ["great_rhombicosidodecahedron", function() great_rhombicosidodecahedron()],
    ["snub_cube", function() snub_cube()],
    ["snub_dodecahedron", function() snub_dodecahedron()]
];
