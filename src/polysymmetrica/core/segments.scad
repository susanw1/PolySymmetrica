// ---------------------------------------------------------------------------
// PolySymmetrica - Face segmentation helpers
// Extracts simple face segments from possibly self-intersecting face loops.
//
// Scope of this file:
// - analyze a face loop into simple cells
// - derive geometry cut entries/segments
// - determine visible cells from face-local +Z
//
// It intentionally stops at data/iteration. 3D keep/cut volumes that consume
// this analysis live elsewhere.

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
        kinds = [for (hidx = hs) (len(hedges[hidx]) > 4 ? hedges[hidx][4] : "parent")],
        span_ids = [for (hidx = hs) hedges[hidx][2]],
        cyc = [pts2d, pts3d, edge_ids, kinds, ns, span_ids],
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

// Non-zero winding-number containment for self-intersecting loops.
// Returns true when point is inside under non-zero winding rule.
function _ps_seg_point_in_poly_nonzero(pt, poly, eps=1e-9) =
    let(
        x = pt[0],
        y = pt[1],
        n = len(poly),
        wn = sum([
            for (i = [0:1:n-1])
                let(
                    j = (i + 1) % n,
                    a = poly[i],
                    b = poly[j],
                    ay = a[1], by = b[1],
                    up = (ay <= y) && (by > y),
                    down = (ay > y) && (by <= y),
                    is_left = _ps_seg_orient2(a, b, [x, y])
                )
                up ? ((is_left > eps) ? 1 : 0) :
                down ? ((is_left < -eps) ? -1 : 0) :
                0
        ])
    )
    wn != 0;

function _ps_seg_reverse_keep_first(list) =
    (len(list) <= 1) ? list : concat([list[0]], [for (i = [len(list)-1:-1:1]) list[i]]);

function _ps_seg_orient_cell(cell, target_sign, eps=1e-9) =
    let(
        area = _ps_seg_poly_area2(cell[0]),
        same = (abs(area) <= eps) || (area * target_sign > 0)
    )
    same ? cell : [
        _ps_seg_reverse_keep_first(cell[0]),
        _ps_seg_reverse_keep_first(cell[1]),
        _ps_reverse(cell[2]),
        _ps_reverse(cell[3])
    ];

function _ps_seg_point_seg_dist(pt, a, b, eps=1e-12) =
    let(
        ab = [b[0] - a[0], b[1] - a[1]],
        ap = [pt[0] - a[0], pt[1] - a[1]],
        den = ab[0]*ab[0] + ab[1]*ab[1],
        t = (den <= eps) ? 0 : max(0, min(1, (ap[0]*ab[0] + ap[1]*ab[1]) / den)),
        q = [a[0] + t * ab[0], a[1] + t * ab[1]]
    )
    norm([pt[0] - q[0], pt[1] - q[1]]);

function _ps_seg_point_param_on_segment(pt, a, b, eps=1e-9) =
    let(
        dx = b[0] - a[0],
        dy = b[1] - a[1],
        den = dx * dx + dy * dy,
        apx = pt[0] - a[0],
        apy = pt[1] - a[1],
        cross = abs(apx * dy - apy * dx),
        t = (den <= eps) ? undef : ((apx * dx + apy * dy) / den),
        L = sqrt(max(den, 0))
    )
    (is_undef(t) || cross > eps * max(1, L) || t < -eps || t > 1 + eps) ? undef : t;

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

function _ps_seg_node_kind(node3d, face_pts3d_local, eps=1e-9) =
    max([
        for (p = face_pts3d_local)
            (norm(node3d - p) <= eps) ? 1 : 0
    ]) == 1 ? "source_vertex" : "crossing";

function _ps_seg_keep_arr_cell(cell, outer2d, mode="nonzero", eps=1e-8) =
    let(
        input_area = _ps_seg_poly_area2(outer2d),
        input_sign = (input_area >= 0) ? 1 : -1,
        pts2d = cell[0],
        area_c = cell[4],
        probe = _ps_seg_cycle_probe_point(pts2d, eps),
        cycle_sign_ok = (abs(input_area) <= eps) ? true : (area_c * input_sign > eps),
        in_evenodd = _ps_seg_point_in_poly_evenodd(probe, outer2d, eps),
        in_nonzero = _ps_seg_point_in_poly_nonzero(probe, outer2d, eps)
    )
    (mode == "all") ? true :
    (mode == "nonzero") ? (in_nonzero && cycle_sign_ok) :
    in_evenodd;

function _ps_seg_orient_arr_cell(cell, target_sign, eps=1e-9) =
    let(
        area = cell[4],
        same = (abs(area) <= eps) || (area * target_sign > 0)
    )
    same ? cell : [
        _ps_seg_reverse_keep_first(cell[0]),
        _ps_seg_reverse_keep_first(cell[1]),
        _ps_seg_reverse_keep_first(cell[2]),
        _ps_reverse(cell[3]),
        -area
    ];

function _ps_seg_fill_cell_ids_from_arr(arr, mode="nonzero", eps=1e-8) =
    let(
        outer2d = arr[0],
        cells = arr[4]
    )
    [
        for (ci = [0:1:len(cells)-1])
            if (_ps_seg_keep_arr_cell(cells[ci], outer2d, mode, eps))
                ci
    ];

function _ps_seg_arr_occurrences(cells, cell_ids, target_sign, eps=1e-9) =
    [
        for (ci = cell_ids)
            let(
                cell = _ps_seg_orient_arr_cell(cells[ci], target_sign, eps),
                node_ids = cell[2],
                span_ids = cell[3],
                m = len(node_ids)
            )
            for (k = [0:1:m-1])
                [span_ids[k], ci, node_ids[k], node_ids[(k+1)%m]]
    ];

function _ps_seg_occ_count_for_span(occs, span_id) =
    len([for (o = occs) if (o[0] == span_id) 1]);

function _ps_seg_other_cell_for_span(occs, span_id, cell_idx) =
    let(
        hits = [for (o = occs) if (o[0] == span_id && o[1] != cell_idx) o[1]]
    )
    (len(hits) == 0) ? undef : hits[0];

function _ps_seg_boundary_occurrences(arr, filled_cell_ids, eps=1e-9) =
    let(
        outer2d = arr[0],
        cells = arr[4],
        input_area = _ps_seg_poly_area2(outer2d),
        input_sign = (input_area >= 0) ? 1 : -1,
        filled_occs = _ps_seg_arr_occurrences(cells, filled_cell_ids, input_sign, eps),
        all_occs = _ps_seg_arr_occurrences(cells, [for (i = [0:1:len(cells)-1]) i], input_sign, eps)
    )
    [
        for (o = filled_occs)
            if (_ps_seg_occ_count_for_span(filled_occs, o[0]) == 1)
                [o[0], o[1], o[2], o[3], _ps_seg_other_cell_for_span(all_occs, o[0], o[1])]
    ];

function _ps_seg_next_boundary_occ(boundary_occs, nodes, cur_idx, used=[]) =
    let(
        cur = boundary_occs[cur_idx],
        u = cur[2],
        v = cur[3],
        pu = nodes[u][0],
        pv = nodes[v][0],
        rev_ang = atan2(pu[1] - pv[1], pu[0] - pv[0]),
        outs = [
            for (i = [0:1:len(boundary_occs)-1])
                if (!_ps_list_contains(used, i) && boundary_occs[i][2] == v)
                    i
        ],
        pairs = [
            for (i = outs)
                let(
                    q = nodes[boundary_occs[i][3]][0],
                    ang = atan2(q[1] - pv[1], q[0] - pv[0]),
                    delta_raw = rev_ang - ang,
                    delta = (delta_raw < 0) ? (delta_raw + 360) : delta_raw
                )
                [delta, i]
        ]
    )
    (len(pairs) == 0) ? undef : pairs[_ps_seg_min_pair_idx(pairs)][1];

function _ps_seg_trace_boundary_loop(boundary_occs, nodes, start_idx, cur_idx, acc=[], depth=0) =
    let(
        max_depth = len(boundary_occs) + 2,
        cur = boundary_occs[cur_idx],
        acc2 = concat(acc, [cur_idx]),
        next_idx = _ps_seg_next_boundary_occ(boundary_occs, nodes, cur_idx, acc2)
    )
    (depth > max_depth) ? [] :
    ((cur[3] == boundary_occs[start_idx][2]) && (len(acc) > 0)) ? acc2 :
    is_undef(next_idx) ? [] :
    _ps_seg_trace_boundary_loop(boundary_occs, nodes, start_idx, next_idx, acc2, depth + 1);

function _ps_seg_extract_boundary_loops(boundary_occs, nodes, idx=0, used=[], loops=[]) =
    (idx >= len(boundary_occs)) ? loops :
    _ps_list_contains(used, idx) ? _ps_seg_extract_boundary_loops(boundary_occs, nodes, idx + 1, used, loops) :
    let(
        loop = _ps_seg_trace_boundary_loop(boundary_occs, nodes, idx, idx),
        used2 = concat(used, loop),
        loops2 = (len(loop) == 0) ? loops : concat(loops, [loop])
    )
    _ps_seg_extract_boundary_loops(boundary_occs, nodes, idx + 1, used2, loops2);

function _ps_seg_boundary_span_ts(span, a, b) =
    (a == span[1] && b == span[2]) ? [span[4], span[5]] :
    (a == span[2] && b == span[1]) ? [span[5], span[4]] :
    [span[4], span[5]];

function _ps_seg_boundary_span_adj_face_normal_local(adj_face_idx, poly_faces_idx, poly_verts_local) =
    is_undef(adj_face_idx) ? undef :
    ps_face_frame_normal(poly_verts_local, poly_faces_idx[adj_face_idx]);

function _ps_seg_boundary_span_adj_face_dir_span_local(source_ex, ex, ey, ez, adj_face_normal_local, eps=1e-9) =
    is_undef(adj_face_normal_local) ? undef :
    let(
        dir0_world = v_norm(v_cross(adj_face_normal_local, source_ex)),
        dir1_world = v_scale(dir0_world, -1),
        dir_world =
            (norm(dir0_world) <= eps)
                ? undef
                : (v_dot(dir0_world, ez) >= v_dot(dir1_world, ez))
                    ? dir0_world
                    : dir1_world
    )
    is_undef(dir_world) ? undef : [v_dot(dir_world, ex), v_dot(dir_world, ey), v_dot(dir_world, ez)];

function _ps_seg_boundary_span_filled_side(seg2d, filled_cell, eps=1e-9) =
    let(
        probe = is_undef(filled_cell) ? undef : _ps_seg_cycle_probe_point(filled_cell[0], eps),
        orient = is_undef(probe) ? 0 : _ps_seg_orient2(seg2d[0], seg2d[1], probe)
    )
    (orient > eps) ? 1 : (orient < -eps) ? -1 : 0;

/**
 * Function: Build internal dihedral-aware boundary-span site records for the current face.
 * Params: face_pts3d_local (loop in face-local 3D), face_idx (current face index), poly_faces_idx/poly_verts_local (full poly in current face-local coordinates), face_neighbors_idx/face_dihedrals (current face-edge metadata), mode (`"nonzero"`, `"evenodd"`, or `"all"`), eps (geometric tolerance)
 * Returns: list of boundary-span site records `[span_idx, center, ex, ey, ez, span_len, seg2d, loop_idx, source_edge_idx, source_t0, source_t1, kind, filled_cell_idx, other_cell_idx, adj_face_idx, dihedral, adj_face_normal_local, filled_side, adj_face_dir_span_local]`
 * Limitations/Gotchas: internal helper for `place_on_face_boundary_spans(...)`; `filled_side` is `+1` when the filled region lies on the left side of `seg2d`, `-1` on the right, and `0` only for degenerate/ambiguous cases
 */
function _ps_face_boundary_span_sites(face_pts3d_local, face_idx, poly_faces_idx, poly_verts_local, face_neighbors_idx, face_dihedrals, mode="nonzero", eps=1e-8) =
    let(
        arr = ps_face_arrangement(face_pts3d_local, eps),
        cells = arr[4],
        bm = ps_face_boundary_model(face_pts3d_local, mode, eps),
        spans = bm[3]
    )
    [
        for (si = [0:1:len(spans)-1])
            let(
                span = spans[si],
                seg2d = span[0],
                dx = seg2d[1][0] - seg2d[0][0],
                dy = seg2d[1][1] - seg2d[0][1],
                span_len = norm([dx, dy]),
                ex = (span_len <= eps) ? [1, 0, 0] : [dx / span_len, dy / span_len, 0],
                ey = [-ex[1], ex[0], 0],
                ez = [0, 0, 1],
                center = [(seg2d[0][0] + seg2d[1][0]) / 2, (seg2d[0][1] + seg2d[1][1]) / 2, 0],
                source_edge_idx = span[2],
                source_a = is_undef(source_edge_idx) ? undef : face_pts3d_local[source_edge_idx],
                source_b = is_undef(source_edge_idx) ? undef : face_pts3d_local[(source_edge_idx + 1) % len(face_pts3d_local)],
                source_ex = (is_undef(source_a) || is_undef(source_b) || norm(source_b - source_a) <= eps)
                    ? ex
                    : v_norm(source_b - source_a),
                filled_cell_idx = span[6],
                other_cell_idx = span[7],
                adj_face_idx =
                    (span[5] == "source" && !is_undef(source_edge_idx) && !is_undef(face_neighbors_idx) && source_edge_idx < len(face_neighbors_idx))
                        ? face_neighbors_idx[source_edge_idx]
                        : undef,
                dihedral =
                    (span[5] == "source" && !is_undef(source_edge_idx) && !is_undef(face_dihedrals) && source_edge_idx < len(face_dihedrals))
                        ? face_dihedrals[source_edge_idx]
                        : undef,
                adj_face_normal_local = _ps_seg_boundary_span_adj_face_normal_local(adj_face_idx, poly_faces_idx, poly_verts_local),
                filled_cell = is_undef(filled_cell_idx) ? undef : cells[filled_cell_idx],
                filled_side = _ps_seg_boundary_span_filled_side(seg2d, filled_cell, eps),
                adj_face_dir_span_local = _ps_seg_boundary_span_adj_face_dir_span_local(source_ex, ex, ey, ez, adj_face_normal_local, eps)
            )
            [
                si,
                center,
                ex,
                ey,
                ez,
                span_len,
                seg2d,
                span[1],
                source_edge_idx,
                span[3],
                span[4],
                span[5],
                filled_cell_idx,
                other_cell_idx,
                adj_face_idx,
                dihedral,
                adj_face_normal_local,
                filled_side,
                adj_face_dir_span_local
            ]
    ];

/**
 * Function: Build the planar arrangement induced by one face loop.
 * Params: face_pts3d_local (loop in face-local 3D), eps (geometric tolerance)
 * Returns: `[face_pts2d, crossings, nodes, spans, cells]`
 * Limitations/Gotchas: cells are all traced simple arrangement cells before any fill-rule filtering; they are not guaranteed convex
 */
function ps_face_arrangement(face_pts3d_local, eps=1e-8) =
    let(
        n = len(face_pts3d_local),
        face_pts2d = [for (p = face_pts3d_local) [p[0], p[1]]]
    )
    (n < 3) ? [face_pts2d, [], [], [], []] :
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
        nodes3d = _ps_seg_unique_nodes(raw_nodes, eps),
        segs = _ps_seg_segments_from_ts(face_pts3d_local, ts, nodes3d, eps),
        hedges = concat(
            [
                for (si = [0:1:len(segs)-1])
                    let(s = segs[si])
                    [s[0], s[1], si, s[2], "parent", s[3], s[4]]
            ],
            [
                for (si = [0:1:len(segs)-1])
                    let(s = segs[si])
                    [s[1], s[0], si, s[2], "parent", s[4], s[3]]
            ]
        ),
        input_area = _ps_seg_poly_area2(face_pts2d),
        input_sign = (input_area >= 0) ? 1 : -1,
        cycles = _ps_seg_extract_cycles(nodes3d, hedges, input_sign, 0, undef, [], eps),
        cycles_u = _ps_seg_dedupe_cycles(cycles),
        crossings = [
            for (h = hits)
                let(node_idx = _ps_seg_node_index(nodes3d, h[4], eps))
                [h[0], h[1], h[2], h[3], [h[4][0], h[4][1]], node_idx]
        ],
        nodes = [
            for (p = nodes3d)
                [[p[0], p[1]], _ps_seg_node_kind(p, face_pts3d_local, eps)]
        ],
        spans = [
            for (s = segs)
                [
                    [[nodes3d[s[0]][0], nodes3d[s[0]][1]], [nodes3d[s[1]][0], nodes3d[s[1]][1]]],
                    s[0],
                    s[1],
                    s[2],
                    s[3],
                    s[4],
                    "source"
                ]
        ],
        cells = [
            for (c = cycles_u)
                [c[0], c[1], c[4], c[5], _ps_seg_poly_area2(c[0])]
        ]
    )
    [face_pts2d, crossings, nodes, spans, cells];

/**
 * Function: Derive the true filled boundary from a face arrangement.
 * Params: face_pts3d_local (loop in face-local 3D), mode (`"nonzero"`, `"evenodd"`, or `"all"`), eps (geometric tolerance)
 * Returns: `[mode, filled_cell_ids, boundary_loops, boundary_spans]`
 * Limitations/Gotchas: boundary spans are currently sourced from arrangement spans only; later synthetic/cut spans can build on this record shape
 */
function ps_face_boundary_model(face_pts3d_local, mode="nonzero", eps=1e-8) =
    let(
        arr = ps_face_arrangement(face_pts3d_local, eps),
        nodes = arr[2],
        spans = arr[3],
        filled_cell_ids = _ps_seg_fill_cell_ids_from_arr(arr, mode, eps),
        boundary_occs = _ps_seg_boundary_occurrences(arr, filled_cell_ids, eps),
        occ_loops = _ps_seg_extract_boundary_loops(boundary_occs, nodes),
        loop_offsets = _ps_prefix_offsets([for (loop = occ_loops) len(loop)], [0]),
        boundary_loops = [
            for (li = [0:1:len(occ_loops)-1])
                let(
                    loop = occ_loops[li],
                    start_node_idx = boundary_occs[loop[0]][2],
                    pts2d = concat(
                        [[nodes[start_node_idx][0][0], nodes[start_node_idx][0][1]]],
                        [for (k = [0:1:len(loop)-2])
                            let(end_idx = boundary_occs[loop[k]][3])
                            [nodes[end_idx][0][0], nodes[end_idx][0][1]]
                        ]
                    ),
                    span_ids = [for (k = [0:1:len(loop)-1]) loop_offsets[li] + k]
                )
                [pts2d, span_ids]
        ],
        boundary_spans = [
            for (li = [0:1:len(occ_loops)-1])
                let(loop = occ_loops[li])
                for (k = [0:1:len(loop)-1])
                    let(
                        occ = boundary_occs[loop[k]],
                        span = spans[occ[0]],
                        a = occ[2],
                        b = occ[3],
                        seg2d = [[nodes[a][0][0], nodes[a][0][1]], [nodes[b][0][0], nodes[b][0][1]]],
                        ts = _ps_seg_boundary_span_ts(span, a, b)
                    )
                    [seg2d, li, span[3], ts[0], ts[1], span[6], occ[1], occ[4]]
        ]
    )
    [mode, filled_cell_ids, boundary_loops, boundary_spans];

/**
 * Function: Split a face loop into simple face-local cells.
 * Params: face_pts3d_local (loop in face-local 3D), mode (`"nonzero"`, `"evenodd"`, or `"all"`), eps (geometric tolerance)
 * Returns: `[[seg_pts2d, seg_pts3d_local, seg_parent_edge_ids, seg_edge_kinds], ...]`
 * Limitations/Gotchas: `mode="nonzero"` is the intended default for solid self-crossing faces; use `"evenodd"` only when parity fill is genuinely wanted
 */
function ps_face_segments(face_pts3d_local, mode="nonzero", eps=1e-8) =
    let(
        n = len(face_pts3d_local),
        arr = ps_face_arrangement(face_pts3d_local, eps),
        outer2d = arr[0],
        spans = arr[3],
        cells = arr[4]
    )
    (n < 3) ? [] :
    let(
        filtered = [
            for (c = cells)
                let(
                    pts2d = c[0],
                    pts3d = c[1],
                    span_ids = c[3],
                    edge_ids = [for (si = span_ids) spans[si][3]],
                    kinds = [for (si = span_ids) spans[si][6]]
                )
                if (_ps_seg_keep_arr_cell(c, outer2d, mode, eps))
                    [pts2d, pts3d, edge_ids, kinds]
        ]
    )
    (len(filtered) == 0)
        ? [[[for (p = face_pts3d_local) [p[0], p[1]]], face_pts3d_local, [for (i = [0:1:n-1]) i], [for (_i = [0:1:n-1]) "parent"]]]
        : filtered;

/**
 * Module: Iterate simple face-local cells for the current placed face.
 * Params: mode (`"nonzero"`, `"evenodd"`, or `"all"`), eps (geometric tolerance)
 * Returns: none; exposes `$ps_seg_*` metadata and calls children once per cell
 */
module place_on_face_segments(mode="nonzero", eps=1e-8) {
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

/**
 * Module: Iterate dihedral-aware boundary spans for the current placed face.
 * Params: mode (`"nonzero"`, `"evenodd"`, or `"all"`), eps (geometric tolerance)
 * Returns: none; exposes `$ps_boundary_span_*` metadata and places children in a face-local span frame
 */
module place_on_face_boundary_spans(mode="nonzero", eps=1e-8) {
    assert(!is_undef($ps_face_pts3d_local), "place_on_face_boundary_spans: requires place_on_faces context ($ps_face_pts3d_local)");
    assert(!is_undef($ps_face_idx), "place_on_face_boundary_spans: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "place_on_face_boundary_spans: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "place_on_face_boundary_spans: requires place_on_faces context ($ps_poly_verts_local)");
    assert(!is_undef($ps_face_neighbors_idx), "place_on_face_boundary_spans: requires place_on_faces context ($ps_face_neighbors_idx)");
    assert(!is_undef($ps_face_dihedrals), "place_on_face_boundary_spans: requires place_on_faces context ($ps_face_dihedrals)");
    sites = _ps_face_boundary_span_sites(
        $ps_face_pts3d_local,
        $ps_face_idx,
        $ps_poly_faces_idx,
        $ps_poly_verts_local,
        $ps_face_neighbors_idx,
        $ps_face_dihedrals,
        mode,
        eps
    );
    for (site = sites) {
        $ps_boundary_span_idx = site[0];
        $ps_boundary_span_count = len(sites);
        $ps_boundary_span_len = site[5];
        $ps_boundary_span_segment2d_local = site[6];
        $ps_boundary_span_loop_idx = site[7];
        $ps_boundary_span_source_edge_idx = site[8];
        $ps_boundary_span_source_t0 = site[9];
        $ps_boundary_span_source_t1 = site[10];
        $ps_boundary_span_kind = site[11];
        $ps_boundary_span_filled_cell_idx = site[12];
        $ps_boundary_span_other_cell_idx = site[13];
        $ps_boundary_span_adj_face_idx = site[14];
        $ps_boundary_span_dihedral = site[15];
        $ps_boundary_span_adj_face_normal_local = site[16];
        $ps_boundary_span_filled_side = site[17];
        $ps_boundary_span_adj_face_dir_span_local = site[18];

        multmatrix(ps_frame_matrix(site[1], site[2], site[3], site[4]))
            children();
    }
}

/**
 * Module: Render a safe 2D polygon fill for concave or self-intersecting loops.
 * Params: points (2D loop), mode (`"nonzero"`, `"evenodd"`, or `"all"`), eps (geometric tolerance)
 * Returns: none; emits 2D filled geometry via the segmented cells
 * Limitations/Gotchas: intended as a drop-in replacement for raw `polygon(points=...)` when backend behavior on crossing loops is unreliable
 */
module ps_polygon(points, mode="nonzero", eps=1e-8) {
    assert(!is_undef(points), "ps_polygon: points must be defined");
    assert(len(points) >= 3, "ps_polygon: need at least 3 points");
    pts3d = [for (p = points) [p[0], p[1], 0]];
    segs = ps_face_segments(pts3d, mode, eps);
    union() {
        for (s = segs)
            polygon(points = s[0]);
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
// 1) splitting into simple segments via ps_face_segments(..., mode),
// 2) ear-clipping each simple segment.
function _ps_seg_face_tris3(face_idx_loop, poly_verts_local, eps=1e-8, mode="evenodd") =
    let(
        face_pts3d = [for (vi = face_idx_loop) poly_verts_local[vi]],
        segs = ps_face_segments(face_pts3d, mode, eps)
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

function _ps_seg_interp2(a, b, t) =
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];

function _ps_seg_source_segments2d(poly2d, cut_segs2d) =
    concat(
        [
            for (i = [0:1:len(poly2d)-1])
                [poly2d[i], poly2d[(i+1)%len(poly2d)], true, i]
        ],
        [
            for (i = [0:1:len(cut_segs2d)-1])
                [cut_segs2d[i][0], cut_segs2d[i][1], false, i]
        ]
    );

function _ps_seg_hits_on_segments2d(src_segs, eps=1e-9) =
    [
        for (i = [0:1:len(src_segs)-1])
            for (j = [i+1:1:len(src_segs)-1])
                let(
                    s0 = src_segs[i],
                    s1 = src_segs[j],
                    a0 = s0[0], b0 = s0[1],
                    a1 = s1[0], b1 = s1[1],
                    hit = _ps_seg_proper_intersection(
                        [a0[0], a0[1], 0], [b0[0], b0[1], 0],
                        [a1[0], a1[1], 0], [b1[0], b1[1], 0],
                        eps
                    ),
                    t_a0_on_1 = _ps_seg_point_param_on_segment(a0, a1, b1, eps),
                    t_b0_on_1 = _ps_seg_point_param_on_segment(b0, a1, b1, eps),
                    t_a1_on_0 = _ps_seg_point_param_on_segment(a1, a0, b0, eps),
                    t_b1_on_0 = _ps_seg_point_param_on_segment(b1, a0, b0, eps)
                ) each concat(
                    !is_undef(hit) ? [[i, hit[0], j, hit[1], hit[2]]] : [],
                    (!is_undef(t_a0_on_1) && t_a0_on_1 > eps && t_a0_on_1 < 1 - eps) ? [[i, 0, j, t_a0_on_1, [a0[0], a0[1], 0]]] : [],
                    (!is_undef(t_b0_on_1) && t_b0_on_1 > eps && t_b0_on_1 < 1 - eps) ? [[i, 1, j, t_b0_on_1, [b0[0], b0[1], 0]]] : [],
                    (!is_undef(t_a1_on_0) && t_a1_on_0 > eps && t_a1_on_0 < 1 - eps) ? [[i, t_a1_on_0, j, 0, [a1[0], a1[1], 0]]] : [],
                    (!is_undef(t_b1_on_0) && t_b1_on_0 > eps && t_b1_on_0 < 1 - eps) ? [[i, t_b1_on_0, j, 1, [b1[0], b1[1], 0]]] : []
                )
    ];

function _ps_seg_src_ts(src_segs, hits, eps=1e-9) =
    [
        for (i = [0:1:len(src_segs)-1])
            _ps_seg_sort_uniq(
                concat(
                    [0, 1],
                    [for (h = hits) if (h[0] == i) h[1]],
                    [for (h = hits) if (h[2] == i) h[3]]
                ),
                eps
            )
    ];

function _ps_seg_nodes_from_src(src_segs, src_ts) =
    [
        for (si = [0:1:len(src_segs)-1])
            let(
                a = src_segs[si][0],
                b = src_segs[si][1],
                ts = src_ts[si]
            )
            for (t = ts)
                [_ps_seg_interp2(a, b, t)[0], _ps_seg_interp2(a, b, t)[1], 0]
    ];

function _ps_seg_subsegments_from_src(src_segs, src_ts, nodes, eps=1e-9) =
    [
        for (si = [0:1:len(src_segs)-1])
            let(
                s = src_segs[si],
                a = s[0],
                b = s[1],
                ts = src_ts[si]
            )
            for (k = [0:1:len(ts)-2])
                let(
                    t0 = ts[k],
                    t1 = ts[k+1],
                    p0 = [_ps_seg_interp2(a, b, t0)[0], _ps_seg_interp2(a, b, t0)[1], 0],
                    p1 = [_ps_seg_interp2(a, b, t1)[0], _ps_seg_interp2(a, b, t1)[1], 0],
                    u = _ps_seg_node_index(nodes, p0, eps),
                    v = _ps_seg_node_index(nodes, p1, eps)
                )
                if (!is_undef(u) && !is_undef(v) && u != v && (t1 - t0) > eps)
                    [u, v, si, s[2], s[3], t0, t1]
    ];

function _ps_seg_split_simple_loop(poly2d, cut_segs2d, eps=1e-8) =
    (len(poly2d) < 3) ? [] :
    let(
        src = _ps_seg_source_segments2d(poly2d, cut_segs2d),
        hits = _ps_seg_hits_on_segments2d(src, eps),
        src_ts = _ps_seg_src_ts(src, hits, eps),
        nodes = _ps_seg_unique_nodes(_ps_seg_nodes_from_src(src, src_ts), eps),
        segs = _ps_seg_subsegments_from_src(src, src_ts, nodes, eps),
        hedges = concat(
            [for (si = [0:1:len(segs)-1]) let(s = segs[si]) [s[0], s[1], si, s[4], s[3] ? "parent" : "cut", s[5], s[6]]],
            [for (si = [0:1:len(segs)-1]) let(s = segs[si]) [s[1], s[0], si, s[4], s[3] ? "parent" : "cut", s[6], s[5]]]
        ),
        input_area = _ps_seg_poly_area2(poly2d),
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
                    inside = _ps_seg_point_in_poly_evenodd(probe, poly2d, eps)
                )
                if (inside) [pts2d, pts3d, edge_ids, kinds]
        ]
    )
    (len(filtered) == 0) ? [[poly2d, [for (p = poly2d) [p[0], p[1], 0]], [for (i = [0:1:len(poly2d)-1]) i], [for (_i = [0:1:len(poly2d)-1]) "parent"]]] : filtered;

function _ps_seg_tri_z_at_xy(pt2d, tri, eps=1e-9) =
    let(
        a = tri[0],
        b = tri[1],
        c = tri[2],
        ax = a[0], ay = a[1], az = a[2],
        bx = b[0], by = b[1], bz = b[2],
        cx = c[0], cy = c[1], cz = c[2],
        px = pt2d[0], py = pt2d[1],
        den = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy)
    )
    (abs(den) <= eps) ? undef :
    let(
        u = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) / den,
        v = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) / den,
        w = 1 - u - v
    )
    ((u < -eps) || (v < -eps) || (w < -eps)) ? undef : (u * az + v * bz + w * cz);

function _ps_seg_pt_occluded(pt2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero") =
    len([
        for (fj = [0:1:len(poly_faces_idx)-1])
            if (fj != face_idx)
                let(
                    tris3 = _ps_seg_face_tris3(poly_faces_idx[fj], poly_verts_local, eps, mode)
                )
                for (tri = tris3)
                    let(z = _ps_seg_tri_z_at_xy(pt2d, tri, eps))
                    if (!is_undef(z) && z > eps) 1
    ]) > 0;

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

function _ps_seg_cut_dihed_from_tri(tri, eps=1e-9) =
    let(
        a = tri[0],
        b = tri[1],
        c = tri[2],
        n = v_norm(v_cross(b - a, c - a)),
        nz = min(1, max(-1, abs(n[2])))
    )
    180 - acos(nz);

function _ps_seg_cut_entries_dedupe(entries, eps=1e-8, i=0, acc=[]) =
    (i >= len(entries)) ? acc :
    let(
        e = entries[i],
        hit = len([for (a = acc) if (_ps_seg2_eq(a[0], e[0], eps)) 1]) > 0,
        acc2 = hit ? acc : concat(acc, [e])
    )
    _ps_seg_cut_entries_dedupe(entries, eps, i + 1, acc2);

/**
 * Function: Derive geometry cut entries where other faces cross the current face plane.
 * Params: face_pts2d (current face loop), face_idx (current face index), poly_faces_idx/poly_verts_local (full poly in current face-local coordinates), eps (tolerance), mode (cutter triangulation fill rule), filter_parent (drop cuts that coincide with parent edges)
 * Returns: `[[seg2d, cutter_face_idx, cut_dihed], ...]`
 */
function ps_face_geom_cut_entries(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero", filter_parent=true) =
    (is_undef(face_pts2d) || is_undef(poly_faces_idx) || is_undef(poly_verts_local)) ? [] :
    let(
        raw = [
            for (fj = [0:1:len(poly_faces_idx)-1])
                if (fj != face_idx)
                    let(
                        f = poly_faces_idx[fj],
                        tris3 = _ps_seg_face_tris3(f, poly_verts_local, eps, mode)
                    )
                    for (tri = tris3)
                        let(
                            seg = _ps_seg_tri_plane_segment(tri[0], tri[1], tri[2], eps),
                            dihed = _ps_seg_cut_dihed_from_tri(tri, eps)
                        )
                        if (!is_undef(seg) && norm([seg[1][0]-seg[0][0], seg[1][1]-seg[0][1]]) > eps)
                            [seg, fj, dihed]
        ],
        uniq = _ps_seg_cut_entries_dedupe(raw, eps),
        out = [
            for (e = uniq)
                if (!(filter_parent && _ps_seg_is_parent_edge(e[0], face_pts2d, eps)))
                    e
        ]
    )
    out;

/**
 * Function: Return only the 2D cut segments from `ps_face_geom_cut_entries(...)`.
 * Params: face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps, mode, filter_parent
 * Returns: `[seg2d, ...]`
 */
function ps_face_geom_cut_segments(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero", filter_parent=true) =
    [for (e = ps_face_geom_cut_entries(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps, mode, filter_parent)) e[0]];

/**
 * Function: Split the current face by geometry cuts and keep only the cells visible from local `+Z`.
 * Params: face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps, mode, filter_parent
 * Returns: `[[cell_pts2d, cell_pts3d_local, cell_edge_ids, cell_edge_kinds], ...]`
 */
function ps_face_visible_segments(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero", filter_parent=true) =
    let(
        base_segs = ps_face_segments([for (p = face_pts2d) [p[0], p[1], 0]], mode, eps),
        cut_segs = ps_face_geom_cut_segments(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps, mode, filter_parent),
        target_sign = (_ps_seg_poly_area2(face_pts2d) >= 0) ? 1 : -1
    )
    [
        for (base = base_segs)
            let(
                cells = _ps_seg_split_simple_loop(base[0], cut_segs, eps)
            )
            for (cell = cells)
                let(
                    probe = _ps_seg_cycle_probe_point(cell[0], eps),
                    hidden = _ps_seg_pt_occluded(probe, face_idx, poly_faces_idx, poly_verts_local, eps, mode)
                )
                if (!hidden) _ps_seg_orient_cell(cell, target_sign, eps)
    ];

/**
 * Module: Iterate geometry-derived cut segments for the current placed face.
 * Params: mode (cutter triangulation fill rule), eps (tolerance), filter_parent (drop cuts that coincide with parent edges)
 * Returns: none; exposes `$ps_face_cut_*` metadata and calls children once per cut segment
 */
module place_on_face_geom_cut_segments(mode="nonzero", eps=1e-8, filter_parent=true) {
    assert(!is_undef($ps_face_pts2d), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_idx), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "place_on_face_geom_cut_segments: requires place_on_faces context ($ps_poly_verts_local)");
    segs = ps_face_geom_cut_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, mode, filter_parent);
    for (si = [0:1:len(segs)-1]) {
        $ps_face_cut_idx = si;
        $ps_face_cut_count = len(segs);
        $ps_face_cut_segment2d_local = segs[si];
        $ps_face_cut_segments2d_local = segs;
        children();
    }
}

/**
 * Module: Iterate the retained visible cells for the current placed face.
 * Params: mode (cell/cutter fill rule), eps (tolerance), filter_parent (drop cuts that coincide with parent edges)
 * Returns: none; exposes `$ps_vis_seg_*` metadata and calls children once per visible cell
 */
module place_on_face_visible_segments(mode="nonzero", eps=1e-8, filter_parent=true) {
    assert(!is_undef($ps_face_pts2d), "place_on_face_visible_segments: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_idx), "place_on_face_visible_segments: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "place_on_face_visible_segments: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "place_on_face_visible_segments: requires place_on_faces context ($ps_poly_verts_local)");

    segs = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, mode, filter_parent);
    for (si = [0:1:len(segs)-1]) {
        s = segs[si];
        $ps_vis_seg_idx = si;
        $ps_vis_seg_count = len(segs);
        $ps_vis_seg_vertex_count = len(s[0]);
        $ps_vis_seg_pts2d = s[0];
        $ps_vis_seg_pts3d_local = s[1];
        $ps_vis_seg_edge_ids = s[2];
        $ps_vis_seg_edge_kinds = s[3];
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
module face_cut_stencil(face_thk, kerf=0.2, extend=0.5, z_pad=0.2, mode="nonzero", eps=1e-8, filter_parent=true) {
    assert(face_thk > 0, "face_cut_stencil: face_thk must be > 0");
    assert(!is_undef($ps_face_pts2d), "face_cut_stencil: requires place_on_faces context");
    assert(!is_undef($ps_face_idx), "face_cut_stencil: requires place_on_faces context");
    assert(!is_undef($ps_poly_faces_idx), "face_cut_stencil: requires place_on_faces context");
    assert(!is_undef($ps_poly_verts_local), "face_cut_stencil: requires place_on_faces context");

    segs = ps_face_geom_cut_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, mode, filter_parent);
    union() {
        for (s = segs)
            _ps_face_cut_strip(s, face_thk, kerf, extend, z_pad, eps);
    }
}

// Extruded keep-mask built from the visible cells of the current face.
// Use this in intersection() with a finished face body to keep only printable
// face pieces instead of subtracting an oversized stencil hull.
module face_visible_mask(face_thk, z_pad=0.2, mode="nonzero", eps=1e-8, filter_parent=true) {
    assert(face_thk > 0, "face_visible_mask: face_thk must be > 0");
    assert(!is_undef($ps_face_pts2d), "face_visible_mask: requires place_on_faces context");
    assert(!is_undef($ps_face_idx), "face_visible_mask: requires place_on_faces context");
    assert(!is_undef($ps_poly_faces_idx), "face_visible_mask: requires place_on_faces context");
    assert(!is_undef($ps_poly_verts_local), "face_visible_mask: requires place_on_faces context");

    segs = ps_face_visible_segments($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, mode, filter_parent);
    linear_extrude(height = face_thk + 2 * z_pad, center = true)
        union() {
            for (s = segs)
                polygon(points = s[0]);
        }
}
