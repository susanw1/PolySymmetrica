// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <../../core/placement.scad>
use <../../models/regular_all.scad>

// Inter-radius - sets the size of the shapes
IR = 30;

DISPLAY_OFFSET_EXAMPLE = 200; // set to 0 to overlay edge examples onto other examples

COLORS = [ "black", "red", "green", "yellow", "blue", "magenta", "cyan", "white", "orange" ];
module color_map(L, i) { color(COLORS[L[i]]) children(); }

translate([0,0,DISPLAY_OFFSET_EXAMPLE]) {
    // Triangular-faced regular polyhedra
    translate([-100, 0, 0])
    place_on_faces(tetrahedron(), IR) {
        color_map([1, 4, 2, 3], $ps_facet_idx)
            translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }

    translate([0, 0, 0])
    place_on_faces(octahedron(), IR) {
        color_map([1, 4, 2, 3, 2, 3, 1, 4], $ps_facet_idx)
            translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }


    translate([100, 0, 0])
    place_on_faces(icosahedron(), IR) {
        color_map([1, 4, 2, 3, 8, 3, 8, 2, 3, 2,    1, 8, 4, 4, 1, 8, 4, 3, 2, 1 ], $ps_facet_idx)
            translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }

    // Duals
    translate([0, 100, 0])
    place_on_faces(hexahedron(), IR) {
        color_map([1, 8, 4, 3, 7, 2 ], $ps_facet_idx)
            translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }

    translate([100, 100, 0])
    place_on_faces(dodecahedron(), IR) {
        color_map([ 1, 4, 3, 4, 4, 2, 3, 2, 2, 3, 1, 1 ], $ps_facet_idx)
            translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }
}
