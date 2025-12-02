use <funcs.scad>
use <placement.scad>

// Ensure face orientation so normal points outward (centroid·normal > 0)
function orient_face_outward(verts, f) =
    let(
        c = face_centroid(verts, f),
        n = face_normal(verts, f)
    )
    (v_dot(c, n) >= 0)
        ? f
        : [ for (i = [len(f)-1 : -1 : 0]) f[i] ];  // reversed

function orient_all_faces_outward(verts, faces) =
    [ for (f = faces) orient_face_outward(verts, f) ];

// Build unique undirected edges from faces
function edges_from_faces(faces) =
    let(
        raw_edges = [
            for (fi = [0 : len(faces)-1])
                let(f = faces[fi])
                    for (k = [0 : len(f)-1])
                        let(
                            a = f[k],
                            b = f[(k+1) % len(f)],
                            e = (a < b) ? [a,b] : [b,a]
                        ) e
        ],
        uniq_edges = [
            for (i = [0 : len(raw_edges)-1])
                let(ei = raw_edges[i])
                    if (sum([
                            for (j = [0 : 1 : i-1])
                                edge_equal(raw_edges[j], ei) ? 1 : 0
                        ]) == 0) ei
        ]
    )
    uniq_edges;

// Does face f contain undirected edge {a,b}?
function face_has_edge(f, a, b) =
    sum([
        for (k = [0 : len(f)-1])
            let(
                x = f[k],
                y = f[(k+1) % len(f)]
            )
            ((x==a && y==b) || (x==b && y==a)) ? 1 : 0
    ]) > 0;
        
        
// For each edge, list the indices of the faces incident to it (should be 2)
function edge_faces_table(faces, edges) =
    [
        for (ei = [0 : len(edges)-1])
            let(e = edges[ei])
            [
                for (fi = [0 : len(faces)-1])
                    if (face_has_edge(faces[fi], e[0], e[1])) fi
            ]
    ];  


function poly_face_centers_unscaled(poly) =
    let(
        faces = poly_faces(poly),
        verts = poly_verts(poly)
    )
    [ for (fi = [0 : len(faces)-1])
        face_centroid(verts, faces[fi])
    ];


function vertex_incident_faces(poly, vi) =
    let(faces = poly_faces(poly))
    [
        for (fi = [0 : len(faces)-1])
            let(f = faces[fi])
            if (sum([
                for (k = [0 : len(f)-1])
                    f[k] == vi ? 1 : 0
            ]) > 0)
                fi
    ];


function edges_incident_to_vertex(edges, v) =
    [
        for (ei = [0 : 1 : len(edges)-1])
            let(e = edges[ei])
            if (e[0] == v || e[1] == v) ei
    ];

function find_edge_index(edges, a, b) =
    let(
        e = (a < b) ? [a,b] : [b,a],
        idxs = [for (i = [0 : len(edges)-1]) if (edge_equal(edges[i], e)) i]
    )
    idxs[0];   // assume the edge exists

function next_face_around_vertex1(v, f_cur, f_prev, faces, edges, edge_faces) =
    let(
        incEdges = edges_incident_to_vertex(edges, v),

        candidates = [
            for (ei = incEdges)
                let(ef = edge_faces[ei])
                // edge has exactly 2 faces and one is f_cur
                if (len(ef) == 2 && (ef[0] == f_cur || ef[1] == f_cur))
                    (ef[0] == f_cur ? ef[1] : ef[0])
        ],

        filtered = [
            for (cf = candidates)
                if (cf != f_prev) cf
        ]
    )
    filtered[0];  // there should always be exactly 1 in a convex poly

function next_face_around_vertex(v, f_cur, f_prev, faces, edges, edge_faces) =
    let(
        f = faces[f_cur],
        n = len(f),

        // position of v in this face
        pos = [for (k = [0 : n-1]) if (f[k] == v) k],
        k0  = pos[0],                 // there should be exactly one

        // the two neighbours of v in this face
        k_prev = (k0 - 1 + n) % n,
        k_next = (k0 + 1) % n,
        v_prev = f[k_prev],
        v_next = f[k_next],

        // edges (v -> v_next) and (v_prev -> v)
        ei1 = find_edge_index(edges, v,      v_next),
        ei2 = find_edge_index(edges, v_prev, v     ),

        ef1 = edge_faces[ei1],
        ef2 = edge_faces[ei2],

        cand1 = (ef1[0] == f_cur ? ef1[1] : ef1[0]),
        cand2 = (ef2[0] == f_cur ? ef2[1] : ef2[0]),

        candidates = [cand1, cand2],
        filtered   = [for (cf = candidates) if (cf != f_prev) cf]
    )
    filtered[0];   // in a convex poly this is unique

function faces_around_vertex_rec(v, f_cur, f_prev, f_start,
                                 faces, edges, edge_faces, acc = []) =
    let(next = next_face_around_vertex(v, f_cur, f_prev, faces, edges, edge_faces))
        (next == f_start)
            ? concat(acc, [f_cur])
            : faces_around_vertex_rec(
                  v,
                  next,
                  f_cur,
                  f_start,
                  faces, edges, edge_faces,
                  concat(acc, [f_cur])
          );


function faces_around_vertex(poly, v, edges, edge_faces) =
    let(
        faces = poly_faces(poly),
        inc   = vertex_incident_faces(poly, v),
        start = inc[0]
    )
    faces_around_vertex_rec(v, start, -1, start, faces, edges, edge_faces);

               
function dual_faces(poly, centers) =
    let(
        faces      = poly_faces(poly),
        edges      = edges_from_faces(faces),
        edge_faces = edge_faces_table(faces, edges),
        verts      = poly_verts(poly)
    )
    [
        for (vi = [0 : len(verts)-1])
            faces_around_vertex(poly, vi, edges, edge_faces)
    ];


function dual_unit_edge_and_e_over_ir(verts, faces) =
    let(
        edges    = edges_from_faces(faces),
        e0       = edges[0],
        vA       = verts[e0[0]],
        vB       = verts[e0[1]],
        unit_e   = norm(vB - vA),
        mid      = (vA + vB) / 2,
        ir       = norm(mid),
        e_over_ir = unit_e / ir
    )
    [unit_e, e_over_ir];


// ---- Generic dual of a convex polyhedron ----
// poly = [ verts, faces, unit_edge, e_over_ir ]
function poly_dual(poly) =
    let(
        // 1) new vertices: face centres of original
        centers      = poly_face_centers_unscaled(poly),

        // 2) raw dual faces: one per original vertex
        faces_raw    = dual_faces(poly, centers),

        // 3) orient faces so normals point outward
        faces_orient = orient_all_faces_outward(centers, faces_raw),

        // 4) compute new unit_edge and e_over_ir
        ue_eir       = dual_unit_edge_and_e_over_ir(centers, faces_orient),
        unit_e_dual  = ue_eir[0],
        e_over_ir_dual = ue_eir[1]
    )
    [ centers, faces_orient, unit_e_dual, e_over_ir_dual ];


//function dodecahedron() = poly_dual(icosahedron());
//function tetrahedron_dual() = poly_dual(tetrahedron());
//function hexahedron() = poly_dual(octahedron());
//
//
//translate([100, 0, 0]) place_on_faces_ir(hexahedron(), 30) 
//    circle(r = $ps_facet_radius, $fn = 4);
//
//place_on_faces_ir(dodecahedron(), 30) 
////    face_debug();
//    circle(r = $ps_facet_radius, $fn = 5);
//
////dodeca_faces_sym(40) {
////    circle(r = $ps_facet_radius, $fn = 5);
////}
