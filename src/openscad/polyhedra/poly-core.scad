use <funcs.scad>
use <polygon-sym.scad>
use <edge-mount.scad>


// ---- Poly descriptor accessors ----
function poly_verts(poly)      = poly[0];
function poly_faces(poly)      = poly[1];
function poly_unit_edge(poly)  = poly[2];
function poly_e_over_ir(poly)  = poly[3];

// ---- Generic face-frame helpers ----
//function poly_face_center(poly, fi, scale) =
//    let(f  = poly_faces(poly)[fi],
//        vs = poly_verts(poly),
//        v0 = vs[f[0]] * scale,
//        v1 = vs[f[1]] * scale,
//        v2 = vs[f[2]] * scale)
//    (v0 + v1 + v2) / 3;

function poly_face_center(poly, fi, scale) =
    let(
        f   = poly_faces(poly)[fi],
        vs  = poly_verts(poly),
        xs  = [ for (vid = f) vs[vid][0] * scale ],
        ys  = [ for (vid = f) vs[vid][1] * scale ],
        zs  = [ for (vid = f) vs[vid][2] * scale ]
    )
    [
        sum(xs) / len(f),
        sum(ys) / len(f),
        sum(zs) / len(f)
    ];


function poly_face_ez(poly, fi, scale) =
    let(f  = poly_faces(poly)[fi],
        vs = poly_verts(poly),
        v0 = vs[f[0]] * scale,
        v1 = vs[f[1]] * scale,
        v2 = vs[f[2]] * scale)
    v_norm(v_cross(v1 - v0, v2 - v0));   // outward normal


function poly_face_ex(poly, fi, scale) =
    let(f      = poly_faces(poly)[fi],
        vs     = poly_verts(poly),
        center = poly_face_center(poly, fi, scale),
        v0     = vs[f[0]] * scale)
    v_norm(v0 - center);   // local +X points to vertex 0

function poly_face_ey(poly, fi, scale) =
    v_cross(
        poly_face_ez(poly, fi, scale),
        poly_face_ex(poly, fi, scale)
    );

function frame_matrix(center, ex, ey, ez) = [
    [ex[0], ey[0], ez[0], center[0]],
    [ex[1], ey[1], ez[1], center[1]],
    [ex[2], ey[2], ez[2], center[2]],
    [0,      0,     0,     1]
];



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
        
        // Vector from face centre to poly centre (which is at world [0,0,0]), expressed in LOCAL coords.
        // World-space vector is -center.
        poly_center_local = [
            v_dot(-center, ex),
            v_dot(-center, ey),
            v_dot(-center, ez)
        ];

        // Per-facet metadata (local-space friendly)
        $ph_facet_idx        = fi;
        $ph_edge_len         = edge_len;
        $ph_face_midradius    = face_midradius;
        $ph_facet_radius      = facet_radius;
        $ph_poly_center_local = poly_center_local;
        
        multmatrix(frame_matrix(center, ex, ey, ez))
            children();
    }
}

// ---- Generic inter-radius driver ----
module place_on_faces_ir(poly, inter_radius) {
    edge_len = inter_radius * poly_e_over_ir(poly);
    place_on_faces(poly, edge_len) children();
}


module face_debug() {
    // Face index
    color("white") translate([0,0,2])
        text(str($ph_facet_idx), size=5, halign="center", valign="center");

    // Local axes
    color("red")   cube([8,1,1], center=false);
    color("green") rotate([0,0,90]) cube([8,1,1], center=false);
    color("blue")  rotate([0,-90,0]) cube([8,1,1], center=false);

    // Radial line to centre
    color("yellow") cylinder(h = -$ph_poly_center_local[2], r = 0.5, center=false);
}



// Return index of some neighbour vertex of vi (from the face list)
function poly_vertex_neighbor(poly, vi) =
    let(
        faces = poly_faces(poly),
        // collect "next" vertex after vi in any face that contains it
        candidates = [
            for (f = faces)
                for (k = [0 : len(f)-1])
                    if (f[k] == vi) f[(k+1) % len(f)]
        ]
    ) candidates[0];  // first one is enough

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
        $ph_vertex_idx        = vi;
        $ph_edge_len          = edge_len;
        $ph_vert_radius       = vert_radius;
        $ph_poly_center_local = [0, 0, -vert_radius];  // by construction

        multmatrix(frame_matrix(center, ex, ey, ez))
            children();
    }
}

// Generic inter-radius driver for vertices
module place_on_vertices_ir(poly, inter_radius) {
    edge_len = inter_radius * poly_e_over_ir(poly);
    place_on_vertices(poly, edge_len) children();
}



// Build unique undirected edge list from faces
function poly_edges(poly) =
    let(
        faces = poly_faces(poly),
        raw_edges = [
            for (fi = [0 : len(faces)-1]) 
                let(f = faces[fi])
                    for (k = [0 : len(f)-1]) 
                        let(
                            a = f[k],
                            b = f[(k+1) % len(f)],
                            e = (a < b) ? [a, b] : [b, a]
                        ) e
        ],
        uniq_edges = [
            for (i = [0 : len(raw_edges)-1])
                let(ei = raw_edges[i])
                    // include if this is the first occurrence
                    if (sum([
                            for (j = [0 : 1 : i-1])
                                edge_equal(raw_edges[j], ei) ? 1 : 0
                        ]) == 0) ei
        ]
    )
    uniq_edges;

// ---- Place children on all edges of a polyhedron ----
module place_on_edges(poly, edge_len) {

    scale = edge_len / poly_unit_edge(poly);
    verts = poly_verts(poly);
    edges = poly_edges(poly);

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
        $ph_edge_idx          = ei;
        $ph_edge_len          = edge_len;
        $ph_edge_midradius    = edge_midradius;
        $ph_poly_center_local = poly_center_local;

        multmatrix(frame_matrix(center, ex, ey, ez))
            children();
    }
}

module place_on_edges_ir(poly, inter_radius) {
    edge_len = inter_radius * poly_e_over_ir(poly);
    place_on_edges(poly, edge_len) children();
}
