// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Validation helpers (structural + geometric)

use <funcs.scad>

// ---- basic helpers ----

function _ps_face_has_distinct_indices(f) =
    let(n = len(f))
    n == len([ for (i=[0:n-1]) if (len([for (j=[0:n-1]) if (f[j]==f[i]) 1]) == 1) 1 ]);

function _ps_faces_min_arity(faces, min_k=3) =
    min([ for (f=faces) len(f) ]) >= min_k;

function _ps_face_edges(f) =
    let(n = len(f))
    [ for (i=[0:1:n-1]) [f[i], f[(i+1)%n]] ];

function _ps_face_has_edge(f, a, b) =
    let(n = len(f))
    max([
        for (i=[0:1:n-1])
            let(u=f[i], v=f[(i+1)%n])
            ((u==a && v==b) || (u==b && v==a)) ? 1 : 0
    ]) == 1;

function _ps_edge_len_ok(verts, e, eps) =
    norm(verts[e[1]] - verts[e[0]]) > eps;

function _ps_faces_edges_nonzero(verts, faces, eps) =
    min([
        for (f=faces)
            min([ for (e=_ps_face_edges(f)) _ps_edge_len_ok(verts, e, eps) ? 1 : 0 ])
    ]) == 1;

// ---- planarity & intersections ----

function _ps_face_planar(verts, f, eps) =
    let(
        n = face_normal(verts, f),
        nlen = norm(n),
        d = v_dot(n, verts[f[0]])
    )
    (nlen > eps) && (min([
        for (vid = f)
            abs(v_dot(n, verts[vid]) - d) <= eps ? 1 : 0
    ]) == 1);

function _ps_faces_planar(verts, faces, eps) =
    min([ for (f = faces) _ps_face_planar(verts, f, eps) ? 1 : 0 ]) == 1;

function _ps_orient2d(a, b, c) =
    (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

function _ps_seg_intersect(a, b, c, d, eps=1e-9) =
    let(
        o1 = _ps_orient2d(a, b, c),
        o2 = _ps_orient2d(a, b, d),
        o3 = _ps_orient2d(c, d, a),
        o4 = _ps_orient2d(c, d, b)
    )
    (o1*o2 < -eps) && (o3*o4 < -eps);

function _ps_face_self_intersections(verts, f, eps=1e-9) =
    let(
        n = face_normal(verts, f),
        c = face_centroid(verts, f),
        a0 = verts[f[0]] - c,
        proj = a0 - n * v_dot(a0, n),
        proj_len = norm(proj),
        ex = proj_len == 0 ? [1,0,0] : proj / proj_len,
        ey = v_cross(n, ex),
        pts2d = [
            for (vid = f)
                let(vdir = verts[vid] - c)
                [ v_dot(vdir, ex), v_dot(vdir, ey) ]
        ],
        m = len(pts2d)
    )
    len([
        for (i = [0:1:m-1])
            for (j = [i+1:1:m-1])
                if (abs(i-j) > 1 && !(i==0 && j==m-1))
                    let(
                        a = pts2d[i],
                        b = pts2d[(i+1)%m],
                        c2 = pts2d[j],
                        d2 = pts2d[(j+1)%m]
                    )
                    if (_ps_seg_intersect(a, b, c2, d2, eps)) 1
    ]);

function _ps_faces_no_self_intersections(verts, faces, eps) =
    max([ for (f = faces) _ps_face_self_intersections(verts, f, eps) ]) == 0;

// ---- manifoldness & convexity ----

function _ps_edges_manifold(verts, faces) =
    let(edges = _ps_edges_from_faces(faces))
    min([
        for (e = edges)
            let(count = sum([ for (f = faces) _ps_face_has_edge(f, e[0], e[1]) ? 1 : 0 ]))
            count == 2 ? 1 : 0
    ]) == 1;

function _ps_faces_outward(verts, faces, eps=1e-9) =
    min([
        for (f = faces)
            let(c = face_centroid(verts, f),
                n = face_normal(verts, f))
            (v_dot(c, n) >= -eps) ? 1 : 0
    ]) == 1;

function _ps_poly_convex(verts, faces, eps=1e-9) =
    let(faces_out = orient_all_faces_outward(verts, faces))
    min([
        for (f = faces_out)
            let(n = face_normal(verts, f),
                d = v_dot(n, verts[f[0]]))
            min([ for (v = verts) v_dot(n, v) <= d + eps ? 1 : 0 ]) == 1 ? 1 : 0
    ]) == 1;

// ---- public validators ----

// mode: "struct" | "closed" | "convex" | "star_ok"
function poly_valid(poly, mode="closed", eps=1e-9) =
    let(
        verts = poly_verts(poly),
        faces = poly_faces(poly),
        base_ok =
            len(verts) >= 4 &&
            len(faces) >= 4 &&
            all_faces_valid(verts, faces) &&
            _ps_faces_min_arity(faces, 3) &&
            min([for (f=faces) _ps_face_has_distinct_indices(f) ? 1 : 0]) == 1 &&
            _ps_faces_edges_nonzero(verts, faces, eps) &&
            _ps_faces_planar(verts, faces, eps),
        manifold_ok = _ps_edges_manifold(verts, faces),
        outward_ok = _ps_faces_outward(verts, faces, eps),
        convex_ok = _ps_poly_convex(verts, faces, eps),
        no_self_intersect = _ps_faces_no_self_intersections(verts, faces, eps)
    )
    (mode == "struct") ? base_ok :
    base_ok &&
    manifold_ok &&
    (mode == "convex" ? (outward_ok && convex_ok && no_self_intersect) :
     mode == "closed" ? no_self_intersect : true);

module assert_poly_valid_mode(poly, mode="closed", eps=1e-9) {
    assert(poly_valid(poly, mode, eps), str("poly invalid (mode=", mode, ")"));
}
