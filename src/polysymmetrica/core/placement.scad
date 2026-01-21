// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>

// ---- Generic face-placement driver ----
module place_on_faces(poly, inter_radius = 1, edge_len = undef) {
    exp_edge_len = is_undef(edge_len)? inter_radius * poly_e_over_ir(poly) : edge_len;
    scale = exp_edge_len;

    verts = poly_verts(poly);
    faces = poly_faces(poly);
    faces0 = ps_orient_all_faces_outward(verts, faces);
    edges = _ps_edges_from_faces(faces0);
    edge_faces = ps_edge_faces_table(faces0, edges);
    face_n = [ for (f = faces0) ps_face_normal(verts, f) ];

    for (fi = [0 : 1 : len(faces)-1]) {
        f      = faces[fi];
        center = poly_face_center(poly, fi, scale);
        ex     = poly_face_ex(poly, fi, scale);
        ey     = poly_face_ey(poly, fi, scale);
        ez     = poly_face_ez(poly, fi, scale);

        face_midradius = norm(center);  // distance from origin to face centre
        rad_vec = [ for (vid = f) norm(verts[vid] * scale - center) ];
        facet_radius = sum(rad_vec) / len(rad_vec);

        // Vector from face centre to polyhedral centre (which is at world [0,0,0]), expressed in LOCAL coords.
        // World-space vector is -center.
        poly_center_local = [
            v_dot(-center, ex),
            v_dot(-center, ey),
            v_dot(-center, ez)
        ];

        // Face boundary vertices in LOCAL coords (about the face centre)
        face_verts_local = [
            for (vid = f)
                let(p = verts[vid] * scale - center)
                    [ v_dot(p, ex), v_dot(p, ey), v_dot(p, ez) ]
        ];
        // 2D projection for polygon(), in face plane coords
        face_pts2d = [
            for (p = face_verts_local)
                [ p[0], p[1] ]
        ];

        // Per-facet metadata (local-space friendly) - mean average values where faces are irregular
        $ps_facet_idx         = fi;                 // index of this facet, 0..N-1
        $ps_edge_len          = exp_edge_len;       // (mean) length of edge
        $ps_vertex_count      = len(face_pts2d);    // vertex count for this facet (length of the $ps_face_pts2d list)
        $ps_face_midradius    = face_midradius;     // (mean) distance of the face polygon centre from polyhedral centre
        $ps_facet_radius      = facet_radius;       // (mean) distance from facet centre to vertices
        $ps_poly_center_local = poly_center_local;  // polyhedral centre in local coords (for regular faces, [0, 0, -$face_midradius])
        $ps_face_pts2d        = face_pts2d;         // [[x,y]...] for polygon()
        
        // Adjacent faces per edge (aligned with face vertex order)
        $ps_facet_neighbors_idx = [
            for (k = [0:1:len(f)-1])
                let(
                    v0 = f[k],
                    v1 = f[(k+1)%len(f)],
                    ei = ps_find_edge_index(edges, v0, v1),
                    adj = edge_faces[ei]
                )
                (len(adj) < 2) ? undef : ((adj[0] == fi) ? adj[1] : adj[0])
        ];
        // Dihedral angles per edge (degrees, aligned with face vertex order)
        $ps_facet_dihedrals = [
            for (k = [0:1:len(f)-1])
                let(
                    v0 = f[k],
                    v1 = f[(k+1)%len(f)],
                    ei = ps_find_edge_index(edges, v0, v1),
                    adj = edge_faces[ei],
                    n0 = face_n[fi],
                    n1 = (len(adj) < 2) ? n0 : face_n[(adj[0] == fi) ? adj[1] : adj[0]],
                    dotn = v_dot(n0, n1),
                    c = (dotn > 1) ? 1 : ((dotn < -1) ? -1 : dotn)
                )
                180 - acos(c)
        ];

        multmatrix(ps_frame_matrix(center, ex, ey, ez))
            children();
    }
}

// ---- Place children on all vertices of a polyhedron ----
module place_on_vertices(poly, inter_radius = 1, edge_len = undef) {
    target_edge_len = is_undef(edge_len)? inter_radius * poly_e_over_ir(poly) : edge_len;
    scale = target_edge_len;

    verts = poly_verts(poly);
    faces = poly_faces(poly);
    edges = _ps_edges_from_faces(faces);

    for (vi = [0 : 1 : len(verts)-1]) {

        v0 = verts[vi] * scale;              // world-space vertex position
        ez = v_norm(v0);                     // outward: from centre to vertex

        // Pick a neighbour to define an in-plane direction
        ni = poly_vertex_neighbor(poly, vi);
        v1 = verts[ni] * scale;
        neighbor_dir = v1 - v0;

        // Project neighbour_dir onto plane perpendicular to ez
        proj = neighbor_dir - ez * v_dot(neighbor_dir, ez);
        proj_len = norm(proj);
        ex = proj_len == 0 ? [1,0,0] : proj / proj_len;  // fallback shouldn't trigger on regulars
        ey = v_cross(ez, ex);

        center      = v0;             // vertex position
        vert_radius = norm(center);   // distance from poly centre

        // New: full neighbor list from edges
        neighbors_idx = [
            for (e = edges)
                if (e[0] == vi) e[1]
                else if (e[1] == vi) e[0]
        ];
        valence = len(neighbors_idx);

        // Neighbor vectors in LOCAL coords
        neighbor_pts_local = [
            for (nj = neighbors_idx)
                let(pw = verts[nj] * scale - v0)
                    [ v_dot(pw, ex), v_dot(pw, ey), v_dot(pw, ez) ]
        ];

        // Metadata for children (local-space friendly)
        $ps_vertex_idx              = vi;
        $ps_vertex_valence          = valence;
        $ps_vertex_neighbors_idx    = neighbors_idx;
        $ps_vertex_neighbor_pts_local = neighbor_pts_local;

        $ps_edge_len                = target_edge_len;      // (target edge length parameter)
        $ps_vert_radius             = vert_radius;
        $ps_poly_center_local       = [0, 0, -vert_radius];  // by construction

        multmatrix(ps_frame_matrix(center, ex, ey, ez))
            children();
    }
}


// ---- Place children on all edges of a polyhedron ----
module place_on_edges(poly, inter_radius = 1, edge_len = undef) {
    scale = is_undef(edge_len)? inter_radius * poly_e_over_ir(poly) : edge_len;
    verts = poly_verts(poly);
    faces = poly_faces(poly);
    edges = _ps_edges_from_faces(faces);

    for (ei = [0 :  1 : len(edges)-1]) {

        e  = edges[ei];
        v0 = verts[e[0]] * scale;
        v1 = verts[e[1]] * scale;

        center = (v0 + v1) / 2;          // edge midpoint (world)
        ex     = v_norm(v1 - v0);        // along edge
        ez     = v_norm(center);         // radial, outward from origin
        ey     = v_cross(ez, ex);        // completes basis

        edge_midradius   = norm(center);
        edge_len_actual  = norm(v1 - v0);

        // vector from edge centre to poly centre in local coords
        poly_center_local = [
            v_dot(-center, ex),
            v_dot(-center, ey),
            v_dot(-center, ez)
        ];

        // Edge endpoints in LOCAL coords
        edge_pts_local = [[-edge_len_actual/2, 0, 0], [edge_len_actual/2, 0, 0]];

        // Adjacent faces (indices)
        adj_faces_idx = [
            for (fi = [0 :  1 : len(faces)-1])
                if (ps_face_has_edge(faces[fi], e[0], e[1])) fi
        ];

        // Metadata for children (edge-local)
        $ps_edge_idx            = ei;
        $ps_edge_len            = edge_len_actual;      // actual length of this edge (vs supplied edge_len = scaling factor arg)
        $ps_edge_midradius      = edge_midradius;
        $ps_poly_center_local   = poly_center_local;

        $ps_edge_pts_local      = edge_pts_local;
        $ps_edge_verts_idx      = e;
        $ps_edge_adj_faces_idx  = adj_faces_idx;

        multmatrix(ps_frame_matrix(center, ex, ey, ez))
            children();
    }
}


module face_debug() {
    // Face index
    color("white") translate([0,0,2])
        text(str($ps_facet_idx), size=5, halign="center", valign="center");

    // Local axes
    color("red")   cube([8,1,1], center=false);
    color("green") rotate([0,0,90]) cube([8,1,1], center=false);
    color("blue")  rotate([0,-90,0]) cube([8,1,1], center=false);

    // Radial line to centre
    color("yellow") cylinder(h = -$ps_poly_center_local[2], r = 0.5, center=false);
}


