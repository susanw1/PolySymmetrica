// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

/**
This is a handy utility file to pull in all 5 Platonic solid definitions, including named functions for dodecahedron and hexahedron (cube).
*/

include <tetrahedron.scad>
include <octahedron.scad>
include <icosahedron.scad>

include <hexahedron.scad>
include <dodecahedron.scad>

/**
Return all Platonic solids as [[name, poly], ...].
Note: this constructs every poly when called.
*/
function platonics_all() = [
    ["tetrahedron", tetrahedron()],
    ["hexahedron", hexahedron()],
    ["octahedron", octahedron()],
    ["dodecahedron", dodecahedron()],
    ["icosahedron", icosahedron()]
];
