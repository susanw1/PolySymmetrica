use <funcs.scad>

// ---- Generic face-placement driver ----
module place_on_faces(poly, edge_len) {
    scale = edge_len / poly_unit_edge(poly);
    verts = poly_verts(poly);
    
    for (fi = [0 : len(poly_faces(poly))-1]) {
        f      = poly_faces(poly)[fi];
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

        // Per-facet metadata (local-space friendly)
        $ps_facet_idx         = fi;
        $ps_edge_len          = edge_len;
        $ps_face_midradius    = face_midradius;
        $ps_facet_radius      = facet_radius;
        $ps_poly_center_local = poly_center_local;  
        
        multmatrix(frame_matrix(center, ex, ey, ez))
            children();
    }
}

// ---- Generic inter-radius driver ----
module place_on_faces_ir(poly, inter_radius) {
    edge_len = inter_radius * poly_e_over_ir(poly);
    place_on_faces(poly, edge_len) children();
}


// ---- Place children on all vertices of a polyhedron ----
module place_on_vertices(poly, edge_len) {

    scale = edge_len / poly_unit_edge(poly);
    verts = poly_verts(poly);

    for (vi = [0 : len(verts)-1]) {

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

        // Metadata for children (local-space friendly)
        $ps_vertex_idx        = vi;
        $ps_edge_len          = edge_len;
        $ps_vert_radius       = vert_radius;
        $ps_poly_center_local = [0, 0, -vert_radius];  // by construction

        multmatrix(frame_matrix(center, ex, ey, ez))
            children();
    }
}

// Generic inter-radius driver for vertices
module place_on_vertices_ir(poly, inter_radius) {
    edge_len = inter_radius * poly_e_over_ir(poly);
    place_on_vertices(poly, edge_len) children();
}


// ---- Place children on all edges of a polyhedron ----
module place_on_edges(poly, edge_len) {

    scale = edge_len / poly_unit_edge(poly);
    verts = poly_verts(poly);
    edges = _ps_edges_from_faces(poly_faces(poly));

    for (ei = [0 : len(edges)-1]) {

        e  = edges[ei];
        v0 = verts[e[0]] * scale;
        v1 = verts[e[1]] * scale;

        center = (v0 + v1) / 2;          // edge midpoint (world)
        ex     = v_norm(v1 - v0);        // along edge
        ez     = v_norm(center);         // radial, outward from origin
        ey     = v_cross(ez, ex);        // completes basis

        edge_midradius    = norm(center);

        // vector from edge centre to poly centre in local coords
        poly_center_local = [
            v_dot(-center, ex),
            v_dot(-center, ey),
            v_dot(-center, ez)
        ];

        // Metadata for children (edge-local)
        $ps_edge_idx          = ei;
        $ps_edge_len          = edge_len;
        $ps_edge_midradius    = edge_midradius;
        $ps_poly_center_local = poly_center_local;

        multmatrix(frame_matrix(center, ex, ey, ez))
            children();
    }
}

module place_on_edges_ir(poly, inter_radius) {
    edge_len = inter_radius * poly_e_over_ir(poly);
    place_on_edges(poly, edge_len) children();
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



