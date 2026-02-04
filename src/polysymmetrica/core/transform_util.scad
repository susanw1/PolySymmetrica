// ---------------------------------------------------------------------------
// PolySymmetrica - transform utilities
// Shared helpers for transform-style operators.

use <funcs.scad>
use <duals.scad>

// Index k where edge f[k]->f[k+1] matches (v0->v1), or -1 if not found.
function _ps_face_edge_index(f, v0, v1) =
    let(
        n = len(f),
        hits = [for (k = [0:1:n-1]) if (f[k] == v0 && f[(k+1)%n] == v1) k]
    )
    (len(hits) == 0) ? -1 : hits[0];

function _ps_face_edge_site(base, k, near_next=false) =
    base + 2*k + (near_next ? 1 : 0);

// Return [s_near_v0, s_near_v1] for edge (v0,v1) on face fidx.
function _ps_face_edge_sites_for_face_edge(faces, fidx, v0, v1, base) =
    let(
        f = faces[fidx],
        k_dir = _ps_face_edge_index(f, v0, v1),
        k = (k_dir >= 0) ? k_dir : _ps_face_edge_index(f, v1, v0),
        flip = (k_dir < 0)
    )
    assert(k >= 0, "cantitruncate: edge not found in face")
    [
        flip ? _ps_face_edge_site(base, k, true) : _ps_face_edge_site(base, k, false),
        flip ? _ps_face_edge_site(base, k, false) : _ps_face_edge_site(base, k, true)
    ];

// Common poly prep: oriented faces, edges, edge->faces, normals, poly0.
function _ps_poly_base(poly) =
    let(
        verts0 = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts0, poly_faces(poly)),
        edges = _ps_edges_from_faces(faces0),
        edge_faces = ps_edge_faces_table(faces0, edges),
        face_n = [ for (f = faces0) ps_face_normal(verts0, f) ],
        poly0 = make_poly(verts0, faces0, poly_e_over_ir(poly))
    )
    [verts0, faces0, edges, edge_faces, face_n, poly0];

// Edge points: two per edge, located at fraction t along each edge.
function _ps_edge_points(verts0, edges, t) =
    [
        for (ei = [0:1:len(edges)-1])
            let(
                a = edges[ei][0],
                b = edges[ei][1],
                A = verts0[a],
                B = verts0[b]
            )
            [ A + t * (B - A), B + t * (A - B) ]
    ];

// Face cycles from face-edge sites (2n-gons).
function _ps_face_cycles_from_face_edge_sites(faces0, face_edge_offsets) =
    [
        for (fi = [0:1:len(faces0)-1])
            let(n = len(faces0[fi]), base = face_edge_offsets[fi])
            [
                for (k = [0:1:n-1])
                    let(
                        s0 = _ps_face_edge_site(base, k, false),
                        s1 = _ps_face_edge_site(base, k, true)
                    )
                    each [[1, s0], [1, s1]]
            ]
    ];

// Edge cycles from face-edge sites (quads).
function _ps_edge_cycles_from_face_edge_sites(faces0, edges, edge_faces, face_edge_offsets) =
    [
        for (ei = [0:1:len(edges)-1])
            let(
                e = edges[ei],
                fpair = edge_faces[ei],
                f0 = fpair[0],
                f1 = fpair[1],
                v0 = e[0],
                v1 = e[1],
                b0 = face_edge_offsets[f0],
                b1 = face_edge_offsets[f1],
                s_pair0 = _ps_face_edge_sites_for_face_edge(faces0, f0, v0, v1, b0),
                s_pair1 = _ps_face_edge_sites_for_face_edge(faces0, f1, v1, v0, b1)
            )
            [
                [1, s_pair0[0]],
                [1, s_pair0[1]],
                [1, s_pair1[0]],
                [1, s_pair1[1]]
            ]
    ];

// Vertex cycles from face-edge sites (2*valence-gons).
function _ps_vert_cycles_from_face_edge_sites(verts0, faces0, edges, edge_faces, face_edge_offsets, poly0) =
    [
        for (vi = [0:1:len(verts0)-1])
            let(fc = faces_around_vertex(poly0, vi, edges, edge_faces))
            [
                for (fi = fc)
                    let(
                        f = faces0[fi],
                        n = len(f),
                        pos = _ps_index_of(f, vi),
                        k_prev = (pos - 1 + n) % n,
                        base = face_edge_offsets[fi],
                        s_prev = _ps_face_edge_site(base, k_prev, true),
                        s_next = _ps_face_edge_site(base, pos, false)
                    )
                    each [[1, s_prev], [1, s_next]]
            ]
    ];
