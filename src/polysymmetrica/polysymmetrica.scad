// ---------------------------------------------------------------------------
// PolySymmetrica — Umbrella Include File
// Version: 0.1.0
//
// This file re-exports the full PolySymmetrica engine so users can
// simply:
//      use <polysymmetrica.scad>
// instead of importing individual modules.
//
// All public API symbols are provided via this file.
// ---------------------------------------------------------------------------

// Core math + utilities
use <core/funcs.scad>
use <core/placement.scad>
use <core/duals.scad>

// Polyhedral models (primitive symmetries + dual models)
use <models/tetrahedron.scad>
use <models/octahedron.scad>
use <models/icosahedron.scad>

// Users may add further models (Archimedeans, Catalans, compounds)
// in the models/ directory, and import them here.
//
