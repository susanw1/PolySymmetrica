// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

/**
Defines the hexahedron (simple cube!) - derived from icosahedron() using poly_dual().
*/

use <../core/duals.scad>
use <octahedron.scad>

function hexahedron() = poly_dual(octahedron());
