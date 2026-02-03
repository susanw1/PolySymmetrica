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
