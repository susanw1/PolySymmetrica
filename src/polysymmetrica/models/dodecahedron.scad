// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

/**
Defines the dodecahedron - derived from icosahedron() using poly_dual().
*/

use <../core/duals.scad>
use <icosahedron.scad>

function dodecahedron() = poly_dual(icosahedron());
