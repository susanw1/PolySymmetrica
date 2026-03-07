// ---------------------------------------------------------------------------
// PolySymmetrica - Face segmentation helpers
// Extracts simple face segments from possibly self-intersecting face loops.

use <funcs.scad>

function _ps_seg_cross2(a, b) = a[0] * b[1] - a[1] * b[0];
function _ps_seg_interp3(a, b, t) = a + (b - a) * t;

function _ps_seg_proper_intersection(a, b, c, d, eps=1e-9) =
    let(
        r = b - a,
        s = d - c,
        den = _ps_seg_cross2([r[0], r[1]], [s[0], s[1]]),
        qmp = c - a,
        ta = (abs(den) <= eps) ? undef : (_ps_seg_cross2([qmp[0], qmp[1]], [s[0], s[1]]) / den),
        tb = (abs(den) <= eps) ? undef : (_ps_seg_cross2([qmp[0], qmp[1]], [r[0], r[1]]) / den),
        ok = !is_undef(ta) && !is_undef(tb) &&
             (ta > eps) && (ta < 1 - eps) &&
             (tb > eps) && (tb < 1 - eps)
    )
    ok ? [ta, tb, _ps_seg_interp3(a, b, ta)] : undef;

function _ps_seg_nonadj(n, i, j) =
    abs(i - j) > 1 && !(i == 0 && j == n - 1);

function _ps_seg_hits(pts3d, eps=1e-9) =
    let(n = len(pts3d))
    [
        for (i = [0:1:n-1])
            for (j = [i+1:1:n-1])
                if (_ps_seg_nonadj(n, i, j))
                    let(
                        a = pts3d[i],
                        b = pts3d[(i+1)%n],
                        c = pts3d[j],
                        d = pts3d[(j+1)%n],
                        hit = _ps_seg_proper_intersection(a, b, c, d, eps)
                    )
                    if (!is_undef(hit))
                        [i, hit[0], j, hit[1], hit[2]]
    ];

function _ps_seg_sort_uniq(vals, eps=1e-9) =
    let(sv = _ps_sort(vals))
    [
        for (i = [0:1:len(sv)-1])
            if (i == 0 || abs(sv[i] - sv[i-1]) > eps) sv[i]
    ];

function _ps_seg_edge_ts(pts3d, hits, eps=1e-9) =
    let(n = len(pts3d))
    [
        for (i = [0:1:n-1])
            _ps_seg_sort_uniq(
                concat(
                    [0, 1],
                    [for (h = hits) if (h[0] == i) h[1]],
                    [for (h = hits) if (h[2] == i) h[3]]
                ),
                eps
            )
    ];

function _ps_seg_unique_nodes(raw_nodes, eps=1e-9, i=0, acc=[]) =
    (i >= len(raw_nodes)) ? acc :
    let(
        p = raw_nodes[i],
        hit = [for (j = [0:1:len(acc)-1]) if (norm([acc[j][0]-p[0], acc[j][1]-p[1]]) <= eps) j],
        acc2 = (len(hit) == 0) ? concat(acc, [p]) : acc
    )
    _ps_seg_unique_nodes(raw_nodes, eps, i + 1, acc2);

function _ps_seg_node_index(nodes, p, eps=1e-9) =
    let(hit = [for (i = [0:1:len(nodes)-1]) if (norm([nodes[i][0]-p[0], nodes[i][1]-p[1]]) <= eps) i])
    (len(hit) == 0) ? undef : hit[0];

function _ps_seg_segments_from_ts(pts3d, edge_ts, nodes, eps=1e-9) =
    let(n = len(pts3d))
    [
        for (ei = [0:1:n-1])
            let(
                a = pts3d[ei],
                b = pts3d[(ei+1)%n],
                ts = edge_ts[ei]
            )
            for (k = [0:1:len(ts)-2])
                let(
                    t0 = ts[k],
                    t1 = ts[k+1],
                    p0 = _ps_seg_interp3(a, b, t0),
                    p1 = _ps_seg_interp3(a, b, t1),
                    u = _ps_seg_node_index(nodes, p0, eps),
                    v = _ps_seg_node_index(nodes, p1, eps)
                )
                if (!is_undef(u) && !is_undef(v) && u != v && (t1 - t0) > eps)
                    [u, v, ei, t0, t1]
    ];

function _ps_seg_min_pair_idx(pairs, idx=0, best=0) =
    (idx >= len(pairs)) ? best :
    _ps_seg_min_pair_idx(pairs, idx + 1, (pairs[idx][0] < pairs[best][0]) ? idx : best);

function _ps_seg_remove_at(list, idx) =
    [for (i = [0:1:len(list)-1]) if (i != idx) list[i]];

function _ps_seg_sort_pairs(pairs, acc=[]) =
    (len(pairs) == 0) ? acc :
    let(mi = _ps_seg_min_pair_idx(pairs))
    _ps_seg_sort_pairs(_ps_seg_remove_at(pairs, mi), concat(acc, [pairs[mi]]));

function _ps_seg_next_hedge(nodes, hedges, cur_hi) =
    let(
        h = hedges[cur_hi],
        v = h[1],
        outs = [for (hi = [0:1:len(hedges)-1]) if (hedges[hi][0] == v) hi],
        pairs = [
            for (hi = outs)
                let(
                    to = hedges[hi][1],
                    dx = nodes[to][0] - nodes[v][0],
                    dy = nodes[to][1] - nodes[v][1]
                )
                [atan2(dy, dx), hi]
        ],
        ord = [for (p = _ps_seg_sort_pairs(pairs)) p[1]],
        rev_hits = [
            for (k = [0:1:len(hedges)-1])
                if (
                    k != cur_hi &&
                    hedges[k][2] == h[2] &&
                    hedges[k][0] == h[1] &&
                    hedges[k][1] == h[0]
                ) k
        ],
        rev = (len(rev_hits) == 0) ? undef : rev_hits[0],
        pos = search(rev, ord),
        m = len(ord)
    )
    (is_undef(rev) || m == 0 || len(pos) == 0) ? undef : ord[(pos[0] - 1 + m) % m];

function _ps_seg_trace_cycle(nodes, hedges, start_hi, cur_hi, acc_h, acc_n, depth, max_depth) =
    (depth > max_depth) ? [[], []] :
    let(
        h = hedges[cur_hi],
        v = h[1],
        acc_h2 = concat(acc_h, [cur_hi]),
        acc_n2 = concat(acc_n, [v]),
        next_hi = _ps_seg_next_hedge(nodes, hedges, cur_hi)
    )
    is_undef(next_hi) ? [[], []] :
    (next_hi == start_hi) ? [acc_h2, acc_n2] :
    _ps_list_contains(acc_h2, next_hi) ? [[], []] :
    _ps_seg_trace_cycle(nodes, hedges, start_hi, next_hi, acc_h2, acc_n2, depth + 1, max_depth);

function _ps_seg_poly_area2(pts2d) =
    let(n = len(pts2d))
    (n < 3) ? 0 :
    sum([for (i = [0:1:n-1]) let(j = (i+1)%n) pts2d[i][0]*pts2d[j][1] - pts2d[j][0]*pts2d[i][1]]) / 2;

function _ps_seg_set_true_many(flags, ids) =
    [for (i = [0:1:len(flags)-1]) flags[i] || _ps_list_contains(ids, i)];

function _ps_seg_extract_cycles(nodes, hedges, input_area_sign, hi=0, visited=undef, cycles=[], eps=1e-9) =
    let(
        vis = is_undef(visited) ? [for (_i = [0:1:len(hedges)-1]) false] : visited
    )
    (hi >= len(hedges)) ? cycles :
    vis[hi] ? _ps_seg_extract_cycles(nodes, hedges, input_area_sign, hi + 1, vis, cycles, eps) :
    let(
        tr = _ps_seg_trace_cycle(nodes, hedges, hi, hi, [], [hedges[hi][0]], 0, len(hedges) + 2),
        hs = tr[0],
        ns_closed = tr[1],
        ns = (len(ns_closed) >= 2) ? [for (i = [0:1:len(ns_closed)-2]) ns_closed[i]] : [],
        pts2d = [for (ni = ns) [nodes[ni][0], nodes[ni][1]]],
        pts3d = [for (ni = ns) nodes[ni]],
        area = _ps_seg_poly_area2(pts2d),
        keep = (len(hs) >= 3) && (len(ns) >= 3) && (abs(area) > eps),
        edge_ids = [for (hidx = hs) hedges[hidx][3]],
        kinds = [for (_k = edge_ids) "parent"],
        cyc = [pts2d, pts3d, edge_ids, kinds, ns],
        vis2 = (len(hs) == 0) ? _ps_list_set(vis, hi, true) : _ps_seg_set_true_many(vis, hs),
        cycles2 = keep ? concat(cycles, [cyc]) : cycles
    )
    _ps_seg_extract_cycles(nodes, hedges, input_area_sign, hi + 1, vis2, cycles2, eps);

function _ps_seg_dedupe_cycles(cycles, i=0, keys=[], acc=[]) =
    (i >= len(cycles)) ? acc :
    let(
        c = cycles[i],
        key = _ps_sort(c[4]),
        seen = len([for (k = keys) if (k == key) 1]) > 0,
        keys2 = seen ? keys : concat(keys, [key]),
        acc2 = seen ? acc : concat(acc, [c])
    )
    _ps_seg_dedupe_cycles(cycles, i + 1, keys2, acc2);

function _ps_seg_point_in_poly_evenodd(pt, poly, eps=1e-9) =
    let(
        x = pt[0],
        y = pt[1],
        n = len(poly),
        hits = [
            for (i = [0:1:n-1])
                let(
                    j = (i + 1) % n,
                    xi = poly[i][0], yi = poly[i][1],
                    xj = poly[j][0], yj = poly[j][1],
                    cond = ((yi > y) != (yj > y)),
                    den = yj - yi,
                    xint = cond ? ((xj - xi) * (y - yi) / (abs(den) <= eps ? ((den >= 0) ? eps : -eps) : den) + xi) : 0
                )
                (cond && x < xint) ? 1 : 0
        ]
    )
    (sum(hits) % 2) == 1;

function _ps_seg_point_seg_dist(pt, a, b, eps=1e-12) =
    let(
        ab = [b[0] - a[0], b[1] - a[1]],
        ap = [pt[0] - a[0], pt[1] - a[1]],
        den = ab[0]*ab[0] + ab[1]*ab[1],
        t = (den <= eps) ? 0 : max(0, min(1, (ap[0]*ab[0] + ap[1]*ab[1]) / den)),
        q = [a[0] + t * ab[0], a[1] + t * ab[1]]
    )
    norm([pt[0] - q[0], pt[1] - q[1]]);

function _ps_seg_point_on_poly_boundary(pt, poly, eps=1e-9) =
    max([
        for (i = [0:1:len(poly)-1])
            let(a = poly[i], b = poly[(i+1)%len(poly)])
            (_ps_seg_point_seg_dist(pt, a, b, eps) <= eps) ? 1 : 0
    ]) == 1;

// Pick a robust interior sample for a cycle (works for concave simple polygons):
// try centroid first, then inward-offset probes from each edge midpoint.
function _ps_seg_cycle_probe_point(poly, eps=1e-9) =
    let(
        n = len(poly),
        cx = sum([for (p = poly) p[0]]) / n,
        cy = sum([for (p = poly) p[1]]) / n,
        centroid = [cx, cy],
        area = _ps_seg_poly_area2(poly),
        sgn = (area >= 0) ? 1 : -1,
        xs = [for (p = poly) p[0]],
        ys = [for (p = poly) p[1]],
        diag = norm([max(xs) - min(xs), max(ys) - min(ys)]),
        delta = max(10 * eps, 1e-4 * max(1, diag)),
        edge_probes = [
            for (i = [0:1:n-1])
                let(
                    a = poly[i],
                    b = poly[(i+1)%n],
                    dx = b[0] - a[0],
                    dy = b[1] - a[1],
                    L = norm([dx, dy]),
                    mx = (a[0] + b[0]) / 2,
                    my = (a[1] + b[1]) / 2,
                    nx = (L <= eps) ? 0 : sgn * (-dy) / L,
                    ny = (L <= eps) ? 0 : sgn * ( dx) / L
                )
                [mx + delta * nx, my + delta * ny]
        ],
        cands = concat([centroid], edge_probes),
        idxs = [
            for (i = [0:1:len(cands)-1])
                if (_ps_seg_point_in_poly_evenodd(cands[i], poly, eps) &&
                    !_ps_seg_point_on_poly_boundary(cands[i], poly, 10*eps))
                    i
        ]
    )
    (len(idxs) > 0) ? cands[idxs[0]] : centroid;

// Segment a face loop in face-local coordinates.
// Input:  face_pts3d_local = [[x,y,z], ...] in loop order.
// Return: [[seg_pts2d, seg_pts3d_local, seg_parent_edge_ids, seg_edge_kinds], ...]
function ps_face_segments(face_pts3d_local, mode="evenodd", eps=1e-8) =
    let(n = len(face_pts3d_local))
    (n < 3) ? [] :
    let(
        hits = _ps_seg_hits(face_pts3d_local, eps),
        ts = _ps_seg_edge_ts(face_pts3d_local, hits, eps),
        raw_nodes = [
            for (ei = [0:1:n-1])
                let(
                    a = face_pts3d_local[ei],
                    b = face_pts3d_local[(ei+1)%n],
                    ts_e = ts[ei]
                )
                for (t = ts_e)
                    _ps_seg_interp3(a, b, t)
        ],
        nodes = _ps_seg_unique_nodes(raw_nodes, eps),
        segs = _ps_seg_segments_from_ts(face_pts3d_local, ts, nodes, eps),
        hedges = concat(
            [
                for (si = [0:1:len(segs)-1])
                    let(s = segs[si])
                    [s[0], s[1], si, s[2], s[3], s[4]]
            ],
            [
                for (si = [0:1:len(segs)-1])
                    let(s = segs[si])
                    [s[1], s[0], si, s[2], s[4], s[3]]
            ]
        ),
        outer2d = [for (p = face_pts3d_local) [p[0], p[1]]],
        input_area = _ps_seg_poly_area2(outer2d),
        input_sign = (input_area >= 0) ? 1 : -1,
        cycles = _ps_seg_extract_cycles(nodes, hedges, input_sign, 0, undef, [], eps),
        cycles_u = _ps_seg_dedupe_cycles(cycles),
        filtered = [
            for (c = cycles_u)
                let(
                    pts2d = c[0],
                    pts3d = c[1],
                    edge_ids = c[2],
                    kinds = c[3],
                    probe = _ps_seg_cycle_probe_point(pts2d, eps),
                    in_evenodd = _ps_seg_point_in_poly_evenodd(probe, outer2d, eps),
                    keep = (mode == "all") ? true : in_evenodd
                )
                if (keep) [pts2d, pts3d, edge_ids, kinds]
        ]
    )
    (len(filtered) == 0)
        ? [[[for (p = face_pts3d_local) [p[0], p[1]]], face_pts3d_local, [for (i = [0:1:n-1]) i], [for (_i = [0:1:n-1]) "parent"]]]
        : filtered;

// Nested face-segment iterator. Call inside place_on_faces(...).
module place_on_face_segments(mode="evenodd", eps=1e-8) {
    face_pts3d_local = is_undef($ps_face_pts3d_local)
        ? [for (p = $ps_face_pts2d) [p[0], p[1], 0]]
        : $ps_face_pts3d_local;
    segs = ps_face_segments(face_pts3d_local, mode, eps);
    for (si = [0:1:len(segs)-1]) {
        s = segs[si];
        $ps_seg_idx = si;
        $ps_seg_count = len(segs);
        $ps_seg_vertex_count = len(s[0]);
        $ps_seg_pts2d = s[0];
        $ps_seg_pts3d_local = s[1];
        $ps_seg_parent_face_edge_idx = s[2];
        $ps_seg_edge_kind = s[3];
        $ps_face_has_segments = len(segs) > 1;
        children();
    }
}

function _ps_seg_close2(a, b, eps=1e-8) =
    norm([a[0] - b[0], a[1] - b[1]]) <= eps;

function _ps_seg2_eq(s0, s1, eps=1e-8) =
    (_ps_seg_close2(s0[0], s1[0], eps) && _ps_seg_close2(s0[1], s1[1], eps)) ||
    (_ps_seg_close2(s0[0], s1[1], eps) && _ps_seg_close2(s0[1], s1[0], eps));

function _ps_seg_dedupe_segments(segs, eps=1e-8, i=0, acc=[]) =
    (i >= len(segs)) ? acc :
    let(
        s = segs[i],
        hit = len([for (a = acc) if (_ps_seg2_eq(a, s, eps)) 1]) > 0,
        acc2 = hit ? acc : concat(acc, [s])
    )
    _ps_seg_dedupe_segments(segs, eps, i + 1, acc2);

function _ps_seg_is_parent_edge(seg2d, face_pts2d, eps=1e-8) =
    let(
        p = seg2d[0],
        q = seg2d[1],
        n = len(face_pts2d)
    )
    max([
        for (i = [0:1:n-1])
            let(
                a = face_pts2d[i],
                b = face_pts2d[(i+1)%n],
                dx = b[0] - a[0],
                dy = b[1] - a[1],
                L = sqrt(dx*dx + dy*dy),
                ux = (L <= eps) ? 0 : dx / L,
                uy = (L <= eps) ? 0 : dy / L,
                // Signed area (2D cross) against parent edge line.
                cp = abs((p[0]-a[0]) * dy - (p[1]-a[1]) * dx),
                cq = abs((q[0]-a[0]) * dy - (q[1]-a[1]) * dx),
                colinear = (L > eps) && (cp <= eps * max(1, L)) && (cq <= eps * max(1, L)),
                tp = (p[0]-a[0]) * ux + (p[1]-a[1]) * uy,
                tq = (q[0]-a[0]) * ux + (q[1]-a[1]) * uy,
                s0 = min(tp, tq),
                s1 = max(tp, tq),
                ov0 = max(0, s0),
                ov1 = min(L, s1),
                overlaps = (ov1 - ov0) > eps
            )
            (colinear && overlaps) ? 1 : 0
    ]) == 1;

function _ps_seg_orient2(a, b, c) =
    _ps_seg_cross2([b[0] - a[0], b[1] - a[1]], [c[0] - a[0], c[1] - a[1]]);

function _ps_seg_point_in_tri2(p, a, b, c, eps=1e-9) =
    let(
        o0 = _ps_seg_orient2(a, b, p),
        o1 = _ps_seg_orient2(b, c, p),
        o2 = _ps_seg_orient2(c, a, p),
        has_neg = (o0 < -eps) || (o1 < -eps) || (o2 < -eps),
        has_pos = (o0 > eps) || (o1 > eps) || (o2 > eps)
    )
    !(has_neg && has_pos);

function _ps_seg_is_ear(idxs, i, pts2d, sign, eps=1e-9) =
    let(
        m = len(idxs),
        ip = (i - 1 + m) % m,
        in = (i + 1) % m,
        a = pts2d[idxs[ip]],
        b = pts2d[idxs[i]],
        c = pts2d[idxs[in]],
        convex = (sign * _ps_seg_orient2(a, b, c)) > eps,
        inside = len([
            for (j = [0:1:m-1])
                if (j != ip && j != i && j != in)
                    let(p = pts2d[idxs[j]])
                    if (_ps_seg_point_in_tri2(p, a, b, c, eps)) 1
        ]) > 0
    )
    convex && !inside;

function _ps_seg_ear_tris(idxs, pts2d, sign, eps=1e-9, guard=0) =
    let(m = len(idxs))
    (m < 3) ? [] :
    (m == 3) ? [[idxs[0], idxs[1], idxs[2]]] :
    (guard > max(64, m * m))
        ? [for (k = [1:1:m-2]) [idxs[0], idxs[k], idxs[k+1]]]
        : let(
            ear_pos = [for (i = [0:1:m-1]) if (_ps_seg_is_ear(idxs, i, pts2d, sign, eps)) i]
        )
        (len(ear_pos) == 0)
            ? [for (k = [1:1:m-2]) [idxs[0], idxs[k], idxs[k+1]]]
            : let(
                i0 = ear_pos[0],
                tri = [idxs[(i0 - 1 + m) % m], idxs[i0], idxs[(i0 + 1) % m]],
                idxs2 = _ps_seg_remove_at(idxs, i0)
            )
            concat([tri], _ps_seg_ear_tris(idxs2, pts2d, sign, eps, guard + 1));

function _ps_seg_triangulate_simple_poly_idx(pts2d, eps=1e-9) =
    let(
        n = len(pts2d),
        idxs = [for (i = [0:1:n-1]) i],
        area = _ps_seg_poly_area2(pts2d),
        sign = (area >= 0) ? 1 : -1
    )
    (n < 3) ? [] : _ps_seg_ear_tris(idxs, pts2d, sign, eps, 0);

// Triangulate a face in face-local coordinates.
// Unlike a simple fan, this path handles concave/self-intersecting loops by:
// 1) splitting into simple segments via ps_face_segments(..., "all"),
// 2) ear-clipping each simple segment.
function _ps_seg_face_tris3(face_idx_loop, poly_verts_local, eps=1e-8) =
    let(
        face_pts3d = [for (vi = face_idx_loop) poly_verts_local[vi]],
        segs = ps_face_segments(face_pts3d, "all", eps)
    )
    [
        for (s = segs)
            let(
                pts2d = s[0],
                pts3d = s[1],
                tris_idx = _ps_seg_triangulate_simple_poly_idx(pts2d, eps)
            )
            for (t = tris_idx)
                [pts3d[t[0]], pts3d[t[1]], pts3d[t[2]]]
    ];

function _ps_seg_unique_pts3(raw_pts, eps=1e-8, i=0, acc=[]) =
    (i >= len(raw_pts)) ? acc :
    let(
        p = raw_pts[i],
        hit = len([for (q = acc) if (norm(p - q) <= eps) 1]) > 0,
        acc2 = hit ? acc : concat(acc, [p])
    )
    _ps_seg_unique_pts3(raw_pts, eps, i + 1, acc2);

function _ps_seg_max_pair_idx(pairs, idx=0, best=0) =
    (idx >= len(pairs)) ? best :
    _ps_seg_max_pair_idx(pairs, idx + 1, (pairs[idx][0] > pairs[best][0]) ? idx : best);

function _ps_seg_farthest_pair(pts) =
    (len(pts) < 2) ? undef :
    (len(pts) == 2) ? [pts[0], pts[1]] :
    let(
        pairs = [
            for (i = [0:1:len(pts)-1])
                for (j = [i+1:1:len(pts)-1])
                    [norm(pts[i] - pts[j]), pts[i], pts[j]]
        ],
        p = pairs[_ps_seg_max_pair_idx(pairs)]
    )
    [p[1], p[2]];

function _ps_seg_plane_edge_hit(a, b, eps=1e-8) =
    let(
        za = a[2],
        zb = b[2],
        den = za - zb
    )
    (abs(za) <= eps && abs(zb) <= eps) ? undef :
    ((za > eps && zb > eps) || (za < -eps && zb < -eps)) ? undef :
    (abs(den) <= eps) ? undef :
    let(
        t = za / den,
        p = a + (b - a) * t
    )
    (t < -eps || t > 1 + eps) ? undef : p;

function _ps_seg_tri_plane_segment(a, b, c, eps=1e-8) =
    let(
        raw_pts = [
            for (p = [
                _ps_seg_plane_edge_hit(a, b, eps),
                _ps_seg_plane_edge_hit(b, c, eps),
                _ps_seg_plane_edge_hit(c, a, eps)
            ])
                if (!is_undef(p)) p
        ],
        pts = _ps_seg_unique_pts3(raw_pts, eps),
        pair = _ps_seg_farthest_pair(pts)
    )
    (is_undef(pair)) ? undef :
    [[pair[0][0], pair[0][1]], [pair[1][0], pair[1][1]]];

// Compute local 2D cut segments formed by intersections of other faces with this face plane.
// Inputs are expected from place_on_faces(...):
// - face_pts2d: current face loop in local 2D
// - face_idx: current face index
// - poly_faces_idx: full poly faces as index loops
// - poly_verts_local: full poly vertices in current face-local coordinates
function ps_face_geom_cut_segments(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, filter_parent=true) =
    (is_undef(face_pts2d) || is_undef(poly_faces_idx) || is_undef(poly_verts_local)) ? [] :
    let(
        raw = [
            for (fj = [0:1:len(poly_faces_idx)-1])
                if (fj != face_idx)
                    let(
                        f = poly_faces_idx[fj],
                        tris3 = _ps_seg_face_tris3(f, poly_verts_local, eps)
                    )
                    for (tri = tris3)
                        let(
                            a = tri[0],
                            b = tri[1],
                            c = tri[2],
                            seg = _ps_seg_tri_plane_segment(a, b, c, eps)
                        )
                        if (!is_undef(seg) && norm([seg[1][0]-seg[0][0], seg[1][1]-seg[0][1]]) > eps)
                            seg
        ],
        uniq = _ps_seg_dedupe_segments(raw, eps),
        out = [
            for (s = uniq)
                if (!(filter_parent && _ps_seg_is_parent_edge(s, face_pts2d, eps)))
                    s
        ]
    )
    out;

// Iterate geometry-derived cut segments for current face.
module place_on_face_geom_cut_segments(eps=1e-8, filter_parent=true) {
    assert(!is_undef($ps_face_pts2d), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_idx), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_poly_verts_local)");
    segs = ps_face_geom_cut_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, filter_parent);
    for (si = [0:1:len(segs)-1]) {
        $ps_face_cut_idx = si;
        $ps_face_cut_count = len(segs);
        $ps_face_cut_segment2d_local = segs[si];
        $ps_face_cut_segments2d_local = segs;
        children();
    }
}

module _ps_face_cut_strip(seg2d, face_thk, kerf=0.2, extend=0.5, z_pad=0.2, eps=1e-8) {
    a = seg2d[0];
    b = seg2d[1];
    dx = b[0] - a[0];
    dy = b[1] - a[1];
    len2 = sqrt(dx*dx + dy*dy);
    if (len2 > eps) {
        ux = dx / len2;
        uy = dy / len2;
        nx = -uy;
        ny = ux;
        hw = kerf / 2;
        a2 = [a[0] - ux * extend, a[1] - uy * extend];
        b2 = [b[0] + ux * extend, b[1] + uy * extend];
        pts = [
            [a2[0] + nx*hw, a2[1] + ny*hw],
            [b2[0] + nx*hw, b2[1] + ny*hw],
            [b2[0] - nx*hw, b2[1] - ny*hw],
            [a2[0] - nx*hw, a2[1] - ny*hw]
        ];
        linear_extrude(height = face_thk + 2*z_pad, center=true)
            polygon(points = pts);
    }
}

// Build a subtraction body from geometry-derived face cut segments.
// Use inside place_on_faces(...), typically in difference() with a face plate/face polygon.
module face_cut_stencil(face_thk, kerf=0.2, extend=0.5, z_pad=0.2, eps=1e-8, filter_parent=true) {
    assert(face_thk > 0, "face_cut_stencil: face_thk must be > 0");
    assert(!is_undef($ps_face_pts2d), "face_cut_stencil: requires place_on_faces context");
    assert(!is_undef($ps_face_idx), "face_cut_stencil: requires place_on_faces context");
    assert(!is_undef($ps_poly_faces_idx), "face_cut_stencil: requires place_on_faces context");
    assert(!is_undef($ps_poly_verts_local), "face_cut_stencil: requires place_on_faces context");

    segs = ps_face_geom_cut_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, filter_parent);
    union() {
        for (s = segs)
            _ps_face_cut_strip(s, face_thk, kerf, extend, z_pad, eps);
    }
}
