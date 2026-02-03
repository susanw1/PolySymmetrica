use <../core/duals.scad>
use <../models/archimedians_all.scad>

/**
Catalan solids derived as duals of Archimedean solids.
Snub duals are left as stubs for now.
*/

// ---- Duals of truncations ----

function triakis_tetrahedron() = poly_dual(truncated_tetrahedron());
function triakis_octahedron() = poly_dual(truncated_cube());
function tetrakis_hexahedron() = poly_dual(truncated_octahedron());
function triakis_icosahedron() = poly_dual(truncated_dodecahedron());
function pentakis_dodecahedron() = poly_dual(truncated_icosahedron());

// ---- Duals of rectifications ----

function rhombic_dodecahedron() = poly_dual(cuboctahedron());
function rhombic_triacontahedron() = poly_dual(icosidodecahedron());

// ---- Duals of cantellations (small rhombi*) ----

function deltoidal_icositetrahedron() = poly_dual(rhombicuboctahedron());
function deltoidal_hexecontahedron() = poly_dual(rhombicosidodecahedron());

// ---- Duals of cantitruncations (great rhombi*) ----

function disdyakis_dodecahedron() = poly_dual(great_rhombicuboctahedron());
function disdyakis_triacontahedron() = poly_dual(great_rhombicosidodecahedron());

// ---- Snub duals (stubs) ----

function pentagonal_icositetrahedron() = undef; // TODO (dual of snub cube)
function pentagonal_hexecontahedron() = undef; // TODO (dual of snub dodecahedron)

// Archimedean -> Catalan name mapping (for pairing in demos)
function archimedean_catalan_name_pairs() = [
    ["truncated_tetrahedron", "triakis_tetrahedron"],
    ["truncated_cube", "triakis_octahedron"],
    ["truncated_octahedron", "tetrakis_hexahedron"],
    ["truncated_dodecahedron", "triakis_icosahedron"],
    ["truncated_icosahedron", "pentakis_dodecahedron"],
    ["cuboctahedron", "rhombic_dodecahedron"],
    ["icosidodecahedron", "rhombic_triacontahedron"],
    ["rhombicuboctahedron", "deltoidal_icositetrahedron"],
    ["rhombicosidodecahedron", "deltoidal_hexecontahedron"],
    ["great_rhombicuboctahedron", "disdyakis_dodecahedron"],
    ["great_rhombicosidodecahedron", "disdyakis_triacontahedron"],
    ["snub_cube", "pentagonal_icositetrahedron"],
    ["snub_dodecahedron", "pentagonal_hexecontahedron"]
];

function archimedean_to_catalan_name(n) =
    let(
        pairs = archimedean_catalan_name_pairs(),
        idxs = search([n], [for (p = pairs) p[0]], 1)
    )
    (len(idxs) == 0 || idxs[0] < 0) ? undef : pairs[idxs[0]][1];

/**
Return all Catalan solids as [[name, poly], ...].
Note: this constructs every poly when called; snub duals are undef.
*/
function catalans_all() = [
    ["triakis_tetrahedron", triakis_tetrahedron()],
    ["triakis_octahedron", triakis_octahedron()],
    ["tetrakis_hexahedron", tetrakis_hexahedron()],
    ["triakis_icosahedron", triakis_icosahedron()],
    ["pentakis_dodecahedron", pentakis_dodecahedron()],
    ["rhombic_dodecahedron", rhombic_dodecahedron()],
    ["rhombic_triacontahedron", rhombic_triacontahedron()],
    ["deltoidal_icositetrahedron", deltoidal_icositetrahedron()],
    ["deltoidal_hexecontahedron", deltoidal_hexecontahedron()],
    ["disdyakis_dodecahedron", disdyakis_dodecahedron()],
    ["disdyakis_triacontahedron", disdyakis_triacontahedron()],
    ["pentagonal_icositetrahedron", pentagonal_icositetrahedron()],
    ["pentagonal_hexecontahedron", pentagonal_hexecontahedron()]
];
