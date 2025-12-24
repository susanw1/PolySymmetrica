// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <../../core/placement.scad>
use <../../models/regular_all.scad>


/**
Demonstrates the place_on_faces(), place_on_vertices(), place_on_edges() modules for the five regular polyhedra.
*/

DISPLAY_OFFSET_PLACE_ON_FACES_EXAMPLE = -100; // set to 0 to overlay face examples onto vertex example
DISPLAY_OFFSET_PLACE_ON_EDGES_EXAMPLE = 100; // set to 0 to overlay edge examples onto vertex example

// Inter-radius - sets the size of the shapes
IR = 30;

//////////////////////////////////////////////
// PLACE ON FACETS

// Triangular-faced regular polyhedra

translate([0,0,DISPLAY_OFFSET_PLACE_ON_FACES_EXAMPLE]) {
    translate([-100, 0, 0])
    color("green")
    place_on_faces(tetrahedron(), IR) {
        translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }

    translate([0, 0, 0])
    color("blue")
    place_on_faces(octahedron(), IR) {
        translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }

    translate([100, 0, 0])
    color("red")
    place_on_faces(icosahedron(), IR) {
        translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }

    // Duals
    translate([0, 100, 0])
    color("cyan")
    place_on_faces(hexahedron(), IR) {
        translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }

    translate([100, 100, 0])
    color("pink")
    place_on_faces(dodecahedron(), IR) {
        translate([0,0,-1]) cylinder(h = 1, r = $ps_facet_radius, $fn = $ps_vertex_count);
    }
}


//////////////////////////////////////////////
// PLACE ON VERTICES (Gold)

// Triangular-faced regular polyhedra 

translate([-100, 0, 0]) 
color("gold")
place_on_vertices(tetrahedron(), IR) {
    cylinder(h = 5, r = 5, $fn = $ps_vertex_valence);
}

translate([0, 0, 0]) 
color("gold")
place_on_vertices(octahedron(), IR) {
    cylinder(h = 5, r = 5, $fn = $ps_vertex_valence);
}

translate([100, 0, 0]) 
color("gold")
place_on_vertices(icosahedron(), IR) {
    cylinder(h = 5, r = 5, $fn = $ps_vertex_valence);
}

// Duals
translate([0, 100, 0]) 
color("gold")
place_on_vertices(hexahedron(), IR) {
    cylinder(h = 5, r = 5, $fn = $ps_vertex_valence);
}

translate([100, 100, 0]) 
color("gold")
place_on_vertices(dodecahedron(), IR) {
    cylinder(h = 5, r = 5, $fn = $ps_vertex_valence);
}


//////////////////////////////////////////////
// PLACE ON EDGES (Silver)

// Triangular-faced regular polyhedra 

translate([0,0,DISPLAY_OFFSET_PLACE_ON_EDGES_EXAMPLE]) {
    translate([-100, 0, 0]) 
    color("silver")
    place_on_edges(tetrahedron(), IR) {
        rotate([0,90,0]) translate([0,0,-$ps_edge_len/2]) cylinder(h = $ps_edge_len, r = 2, $fn = 20);
    }

    translate([0, 0, 0]) 
    color("silver")
    place_on_edges(octahedron(), IR) {
        rotate([0,90,0]) translate([0,0,-$ps_edge_len/2]) cylinder(h = $ps_edge_len, r = 2, $fn = 20);
    }

    translate([100, 0, 0]) 
    color("silver")
    place_on_edges(icosahedron(), IR) {
        rotate([0,90,0]) translate([0,0,-$ps_edge_len/2]) cylinder(h = $ps_edge_len, r = 2, $fn = 20);
    }

    // Duals
    translate([0, 100, 0]) 
    color("silver")
    place_on_edges(hexahedron(), IR) {
        rotate([0,90,0]) translate([0,0,-$ps_edge_len/2]) cylinder(h = $ps_edge_len, r = 2, $fn = 20);
    }

    translate([100, 100, 0]) 
    color("silver")
    place_on_edges(dodecahedron(), IR) {
        rotate([0,90,0]) translate([0,0,-$ps_edge_len/2]) cylinder(h = $ps_edge_len, r = 2, $fn = 20);
    }
}