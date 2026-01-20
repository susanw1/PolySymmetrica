// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <placement.scad>

module poly_render(poly, inter_radius = 1) {
    let(
        k = (inter_radius * poly_e_over_ir(poly)),
        pts = poly_verts(poly) * k,
        faces_rhr = orient_all_faces_outward(pts, poly_faces(poly)),
        // OpenSCAD expects LHR (clockwise from outside), so flip RHR faces.
        faces_lhr = [for (f = faces_rhr) _ps_reverse(f)]
    )
    polyhedron(
        points = pts,
        faces  = faces_lhr
    );
}


module poly_describe(poly, name = undef, detail = 1) {
    echo ("==========", is_undef(name)? "Polyhedron" : name, "============");
    verts = poly_verts(poly);
    faces = poly_faces(poly);

    echo ("#vertices: ", len(verts));
    echo ("#faces: ", len(faces));
    echo ("#edges: to be implemented");
    echo ("e_over_ir: ", poly_e_over_ir(poly));
    
    if (detail > 0) {
        place_on_faces(poly, 1) {
            let(
                basic1 = str("facet#", $ps_facet_idx, ": (", $ps_vertex_count, " verts) vert_idx: ", faces[$ps_facet_idx]),
                more2 = str(" poly_rad: ", $ps_face_midradius, " verts: ", [for (f = faces[$ps_facet_idx]) verts[f]])
            )
            echo (str("  ", basic1, detail > 1? more2 : ""));
        }
        place_on_vertices(poly, 1) {
            let(
                basic1 = str("vertex#", $ps_vertex_idx, ": (", $ps_vertex_valence, " edges) neighbours_idx: ", $ps_vertex_neighbors_idx),
                more2 = str("vert_idx: ", verts[$ps_vertex_idx])
            )
            echo (str("  ", basic1, detail > 1? more2 : ""));
//        echo (str("  vertex#", $ps_vertex_idx, ": (", $ps_vertex_valence, " edges), poly_rad: ", $ps_vert_radius,
//            " vert_idx: ", verts[$ps_vertex_idx], ", neighbours_idx: ", $ps_vertex_neighbors_idx));
        }
        place_on_edges(poly, 1) {
            let(
                basic1 = str("edge#", $ps_edge_idx, ": len: ", $ps_edge_len, " verts_idx: ", $ps_edge_verts_idx, " faces_idx: ", $ps_edge_adj_faces_idx),
                more2 = str(" ")
            )
            echo (str("  ", basic1, detail > 1? more2 : ""));
//            echo (str("  edge#", $ps_edge_idx, ": len: ", $ps_edge_len, ", poly_rad: ", $ps_edge_midradius,
//                ", verts_idx: ", $ps_edge_verts_idx, ", faces_idx: ", $ps_edge_adj_faces_idx));
        }
    }
    echo ("======================");
}


//// TESTING

use <../models/octahedron.scad>
use <../core/duals.scad>
use <../core/truncation.scad>

poly_describe(octahedron(), "Octahedron");
//poly_describe(poly_dual(poly_truncate(octahedron())), "Dual of Trunc Octa");

poly_render(octahedron(), 30);
