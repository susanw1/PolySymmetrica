// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2026 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>
use <placement.scad>

module poly_render(poly, inter_radius = 1) {
    let(
        k = (inter_radius * poly_e_over_ir(poly)),
        pts = poly_verts(poly) * k,
        // Faces are expected to be LHR (clockwise from outside).
        faces_lhr = poly_faces(poly)
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
    echo ("#edges: ", len(_ps_edges_from_faces(faces)));
    echo ("e_over_ir: ", poly_e_over_ir(poly));
    
    if (detail > 0) {
        place_on_faces(poly, 1) {
            let(
                face_verts = faces[$ps_face_idx],
                basic1 = str("face#", $ps_face_idx, ": (", $ps_vertex_count, " verts) vert_idx: ", face_verts),
                more2 = str(" poly_rad: ", $ps_face_midradius, " verts: ", [for (f = face_verts) verts[f]],
                    " lengths: ", [for (fi = [0:1:len(face_verts)-1]) norm(verts[face_verts[fi]]-verts[face_verts[(fi+1)%len(face_verts)]]) ]),
                max_plane_err = _ps_face_planarity_err(verts, face_verts),
                more3 = str(" max_plane_err: ", max_plane_err)
            )
            echo (str("  ", basic1, detail > 1 ? more2 : "", detail > 2 ? more3 : ""));
        }
        place_on_vertices(poly, 1) {
            let(
                basic1 = str("vertex#", $ps_vertex_idx, ": (", $ps_vertex_valence, " edges) neighbours_idx: ", $ps_vertex_neighbors_idx),
                more2 = str("vert_idx: ", verts[$ps_vertex_idx])
            )
            echo (str("  ", basic1, detail > 1? more2 : ""));
        }
        place_on_edges(poly, 1) {
            let(
                basic1 = str("edge#", $ps_edge_idx, ": len: ", $ps_edge_len, " verts_idx: ", $ps_edge_verts_idx, " faces_idx: ", $ps_edge_adj_faces_idx),
                more2 = str(" ")
            )
            echo (str("  ", basic1, detail > 1? more2 : ""));
        }
    }
    echo ("======================");
}
