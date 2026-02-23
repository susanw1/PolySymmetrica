// ---------------------------------------------------------------------------
// PolySymmetrica - Solver helpers
// Parameter solvers for mixed-face transforms.

use <funcs.scad>
use <duals.scad>

function _ps_unique_ints(a) =
    let(
        s = _ps_sort([for (x = a) x])
    )
    [for (i = [0:1:len(s)-1]) if (i == 0 || s[i] != s[i-1]) s[i]];

// Build sorted unique [a,b] pairs from a list of sizes.
function _ps_pairs_from_sizes(sizes) =
    let(
        s = _ps_sort([for (x = sizes) x]),
        uniq = [for (i = [0:1:len(s)-1]) if (i == 0 || s[i] != s[i-1]) s[i]]
    )
    [for (i = [0:1:len(uniq)-1]) for (j = [i:1:len(uniq)-1]) [uniq[i], uniq[j]]];

function _ps_map_face_c(face_len, c_by_size, default_c=0) =
    let(
        idxs = [for (i = [0:1:len(c_by_size)-1]) if (c_by_size[i][0] == face_len) i]
    )
    (len(idxs) == 0) ? default_c : c_by_size[idxs[0]][1];

function _ps_map_edge_c(face_len, adj_len, c_by_pair, default_c=0) =
    let(
        a = min(face_len, adj_len),
        b = max(face_len, adj_len),
        idxs = [for (i = [0:1:len(c_by_pair)-1]) if (c_by_pair[i][0] == a && c_by_pair[i][1] == b) i]
    )
    (len(idxs) == 0) ? default_c : c_by_pair[idxs[0]][2];

// Solve per-face-family c values by matching dominant family to trig solution.
function solve_cantitruncate_dominant(poly, dominant_size, edge_idx=undef) =
    let(
        verts = poly_verts(poly),
        faces0 = ps_orient_all_faces_outward(verts, poly_faces(poly)),
        poly0 = make_poly(verts, faces0, poly_e_over_ir(poly)),
        edges = _ps_edges_from_faces(faces0),
        face_sizes = [for (f = faces0) len(f)],
        size_set = _ps_unique_ints(face_sizes),
        f_idx = [for (i = [0:1:len(face_sizes)-1]) if (face_sizes[i] == dominant_size) i][0],
        sol = solve_cantitruncate_trig(poly0, f_idx, edge_idx),
        t = sol[0],
        c_dom = sol[1],
        c_by_size = [
            for (sz = size_set)
                let(f_i = [for (i = [0:1:len(face_sizes)-1]) if (face_sizes[i] == sz) i][0])
                (sz == dominant_size)
                    ? [sz, c_dom]
                    : let(sol_i = solve_cantitruncate_trig(poly0, f_i, edge_idx)) [sz, sol_i[1]]
        ]
    )
    [t, c_by_size];

// Solve per-face-family c values and per-edge-family overrides.
function solve_cantitruncate_dominant_edges(poly, dominant_size, edge_idx=undef) =
    let(
        sol = solve_cantitruncate_dominant(poly, dominant_size, edge_idx),
        t = sol[0],
        c_by_size = sol[1],
        c_edge_by_pair = [
            for (p = _ps_pairs_from_sizes([for (c = c_by_size) c[0]]))
                let(
                    c0 = _ps_map_face_c(p[0], c_by_size, 0),
                    c1 = _ps_map_face_c(p[1], c_by_size, 0),
                    c_edge = (c0 + c1) / 2
                )
                [p[0], p[1], c_edge]
        ]
    )
    [t, c_by_size, c_edge_by_pair];

// Trig-based solver for regular bases (one edge type).
// Uses face interior angle and dihedral to compute t and c directly.
function solve_cantitruncate_trig(poly, face_idx=0, edge_idx=undef) =
    let(
        verts = poly_verts(poly),
        faces = ps_orient_all_faces_outward(verts, poly_faces(poly)),
        f = faces[face_idx],
        n = len(f),
        v0 = f[0],
        v_prev = f[(n-1)%n],
        v_next = f[1],
        a0 = v_norm(verts[v_prev] - verts[v0]),
        a1 = v_norm(verts[v_next] - verts[v0]),
        phi = acos(ps_clamp(v_dot(a0, a1), -1, 1)),
        t = 1 / (2 * (1 + sin(phi/2))),
        edges = _ps_edges_from_faces(faces),
        ei = is_undef(edge_idx) ? ps_find_edge_index(edges, v0, v_next) : edge_idx,
        fpair = ps_edge_faces_table(faces, edges)[ei],
        f0 = fpair[0],
        f1 = fpair[1],
        n0 = ps_face_normal(verts, faces[f0]),
        n1 = ps_face_normal(verts, faces[f1]),
        // alpha is angle between outward normals
        alpha = acos(ps_clamp(v_dot(n0, n1), -1, 1)),
        a = norm(verts[v_next] - verts[v0]),
        ir = min([for (e = edges) norm((verts[e[0]] + verts[e[1]]) / 2)]),
        // across-face distance uses sin(alpha/2); for cube alpha=90 so sin=cos
        d_f = (1 - 2*t) * a / (2 * sin(alpha/2)),
        c = abs(d_f) / ir
    )
    [t, c];
