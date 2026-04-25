// ---------------------------------------------------------------------------
// PolySymmetrica - Face-region volume helpers
// Builds positive face-local volumes from filled face boundary spans.

use <funcs.scad>
use <segments.scad>

/**
 * Function: Signed 2D triangle orientation.
 * Params: a/b/c (2D points)
 * Returns: positive for left turn, negative for right turn, zero for colinear
 */
function _ps_fr_orient2(a, b, c) =
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);

/**
 * Function: Intersect 2D lines represented as `n dot p = d`.
 * Params: n0/d0, n1/d1 (line equations), eps (parallel tolerance)
 * Returns: 2D intersection point, or `undef` for near-parallel lines
 */
function _ps_fr_line2_intersect(n0, d0, n1, d1, eps=1e-12) =
    let(det = n0[0]*n1[1] - n0[1]*n1[0])
    (abs(det) < eps) ? undef
  : [
        (d0*n1[1] - n0[1]*d1) / det,
        (n0[0]*d1 - d0*n1[0]) / det
    ];

/**
 * Function: Build a left-hand normal for a 2D point+direction line.
 * Params: line (`[point2d, dir2d, ...]`)
 * Returns: 2D normal vector
 */
function _ps_fr_line_normal(line) =
    [-line[1][1], line[1][0]];

/**
 * Function: Intersect two 2D point+direction lines.
 * Params: line0/line1 (`[point2d, dir2d, ...]`), eps (parallel tolerance)
 * Returns: 2D intersection point, or `undef` for near-parallel lines
 */
function _ps_fr_line_intersection(line0, line1, eps=1e-9) =
    let(
        n0 = _ps_fr_line_normal(line0),
        n1 = _ps_fr_line_normal(line1),
        d0 = v_dot(n0, line0[0]),
        d1 = v_dot(n1, line1[0])
    )
    _ps_fr_line2_intersect(n0, d0, n1, d1, eps);

/**
 * Function: Reconstruct a boundary-span point on its source edge.
 * Params: face_pts3d_local (source face loop), site (boundary-span site), t (source-edge parameter)
 * Returns: face-local 3D point, or `undef` when source edge metadata is missing
 */
function _ps_fr_span_source_point(face_pts3d_local, site, t) =
    let(
        source_edge_idx = site[8],
        n = len(face_pts3d_local),
        a = is_undef(source_edge_idx) ? undef : face_pts3d_local[source_edge_idx],
        b = is_undef(source_edge_idx) ? undef : face_pts3d_local[(source_edge_idx + 1) % n]
    )
    (is_undef(a) || is_undef(b)) ? undef : a + (b - a) * t;

/**
 * Function: Reconstruct the source-edge 3D segment for one boundary-span site.
 * Params: face_pts3d_local (source face loop), site (boundary-span site)
 * Returns: `[p0, p1]` in face-local 3D, falling back to the planar span when source data is absent
 */
function _ps_fr_span_seg3d(face_pts3d_local, site) =
    let(
        p0 = _ps_fr_span_source_point(face_pts3d_local, site, site[9]),
        p1 = _ps_fr_span_source_point(face_pts3d_local, site, site[10]),
        seg2d = site[6]
    )
    (is_undef(p0) || is_undef(p1))
        ? [[seg2d[0][0], seg2d[0][1], 0], [seg2d[1][0], seg2d[1][1], 0]]
        : [p0, p1];

/**
 * Function: Return the current-face in-plane ray out of the filled side of a span.
 * Params: site (boundary-span site)
 * Returns: span-local unit ray `[0,+/-1,0]` pointing outside the filled region
 */
function _ps_fr_span_exterior_ray(site) =
    [0, (site[17] < 0) ? 1 : -1, 0];

/**
 * Function: Return the current-face in-plane ray into the filled side of a span.
 * Params: site (boundary-span site)
 * Returns: span-local unit ray `[0,+/-1,0]` pointing inside the filled region
 */
function _ps_fr_span_filled_ray(site) =
    [0, (site[17] < 0) ? -1 : 1, 0];

/**
 * Function: Compute integer winding number of a 2D loop around a point.
 * Params: pt (2D point), poly (2D loop), eps (orientation tolerance)
 * Returns: integer winding number
 */
function _ps_fr_winding_number(pt, poly, eps=1e-9) =
    let(
        x = pt[0],
        y = pt[1],
        n = len(poly)
    )
    sum([
        for (i = [0:1:n-1])
            let(
                j = (i + 1) % n,
                a = poly[i],
                b = poly[j],
                is_left = _ps_fr_orient2(a, b, pt)
            )
            (a[1] <= y && b[1] > y && is_left > eps) ? 1 :
            (a[1] > y && b[1] <= y && is_left < -eps) ? -1 :
            0
    ]);

/**
 * Function: Build winding signs for arrangement cells under the source face loop.
 * Params: face_pts3d_local (source face loop), cells (face arrangement cells), eps (tolerance)
 * Returns: list of cell winding signs (`+1`, `-1`, or `0`)
 */
function _ps_fr_cell_winding_signs(face_pts3d_local, cells, eps=1e-8) =
    let(face_pts2d = [for (p = face_pts3d_local) [p[0], p[1]]])
    [
        for (cell = cells)
            let(
                probe = _ps_seg_cycle_probe_point(cell[0], eps),
                wn = _ps_fr_winding_number(probe, face_pts2d, eps)
            )
            (wn > 0) ? 1 : (wn < 0) ? -1 : 0
    ];

/**
 * Function: Select the face-plane ray used by anti-interference projection.
 * Params: site (boundary-span site), input_sign (source face signed-area sign), cell_winding_signs (per-arrangement-cell winding signs)
 * Returns: exterior ray for same-winding cells, filled ray for opposite-winding cells
 */
function _ps_fr_span_face_plane_ray(site, input_sign, cell_winding_signs) =
    let(
        cell_idx = site[12],
        cell_sign =
            (is_undef(cell_idx) || cell_idx < 0 || cell_idx >= len(cell_winding_signs))
                ? 0
                : cell_winding_signs[cell_idx],
        same_winding = (cell_sign == 0) || (cell_sign == input_sign)
    )
    same_winding ? _ps_fr_span_exterior_ray(site) : _ps_fr_span_filled_ray(site);

/**
 * Function: Build the anti-interference bisector direction in span-local coords.
 * Params: site (boundary-span site), input_sign (source face signed-area sign), cell_winding_signs (per-arrangement-cell winding signs), eps (zero-length tolerance)
 * Returns: span-local unit direction between the selected face-plane ray and adjacent-face +Z branch
 */
function _ps_fr_span_bisector_dir_span_local(site, input_sign, cell_winding_signs, eps=1e-8) =
    let(
        face_ray = _ps_fr_span_face_plane_ray(site, input_sign, cell_winding_signs),
        adj0 = site[18],
        adj_unit = (is_undef(adj0) || norm(adj0) <= eps) ? undef : v_norm(adj0),
        raw = is_undef(adj_unit) ? [0, 0, 1] : face_ray + adj_unit
    )
    (norm(raw) <= eps) ? [0, 0, 1] : v_norm(raw);

/**
 * Function: Transform a span-local vector into current face-local coordinates.
 * Params: site (boundary-span site), v_span (span-local vector)
 * Returns: face-local vector
 */
function _ps_fr_span_to_face_local(site, v_span) =
    site[2] * v_span[0] + site[3] * v_span[1] + site[4] * v_span[2];

/**
 * Function: Build the anti-interference bisector direction in face-local coords.
 * Params: site (boundary-span site), input_sign (source face signed-area sign), cell_winding_signs (per-arrangement-cell winding signs), eps (zero-length tolerance)
 * Returns: face-local unit direction
 */
function _ps_fr_span_bisector_dir_local(site, input_sign, cell_winding_signs, eps=1e-8) =
    v_norm(_ps_fr_span_to_face_local(site, _ps_fr_span_bisector_dir_span_local(site, input_sign, cell_winding_signs, eps)));

/**
 * Function: Compute scalar projection distance needed to reach a target Z plane.
 * Params: dz (target minus source Z), dir_z (projection direction Z), max_project (optional cap), eps (near-flat tolerance)
 * Returns: scalar offset along the projection direction
 */
function _ps_fr_project_offset(dz, dir_z, max_project=undef, eps=1e-8) =
    let(
        _ok = assert(
            abs(dz) <= eps || abs(dir_z) > eps || !is_undef(max_project),
            "ps_face_anti_interference: projection direction is too close to parallel to target Z plane; set max_project"
        )
    )
    (abs(dz) <= eps) ? 0 :
    (abs(dir_z) <= eps)
        ? (dz >= 0 ? 1 : -1) * abs(max_project)
        : let(
            raw = dz / dir_z,
            cap = is_undef(max_project) ? undef : abs(max_project)
        )
        is_undef(cap) ? raw : ps_clamp(raw, -cap, cap);

/**
 * Function: Report whether a projection offset would be capped.
 * Params: dz (target minus source Z), dir_z (projection direction Z), max_project (optional cap), eps (near-flat tolerance)
 * Returns: boolean
 */
function _ps_fr_project_was_capped(dz, dir_z, max_project=undef, eps=1e-8) =
    is_undef(max_project) ? false :
    (abs(dz) <= eps) ? false :
    (abs(dir_z) <= eps) ? true :
    abs(dz / dir_z) > abs(max_project) + eps;

/**
 * Function: Project one boundary span to a target Z plane as a 2D line.
 * Params: face_pts3d_local (source face loop), site (boundary-span site), z (target face-local Z), max_project (optional cap), eps (tolerance)
 * Returns: `[point2d, dir2d, was_capped, span_idx, source_edge_idx]`
 */
function _ps_fr_project_span_line(face_pts3d_local, site, z, input_sign, cell_winding_signs, max_project=undef, eps=1e-8) =
    let(
        seg3d = _ps_fr_span_seg3d(face_pts3d_local, site),
        mid = (seg3d[0] + seg3d[1]) / 2,
        dir = _ps_fr_span_bisector_dir_local(site, input_sign, cell_winding_signs, eps),
        offset = _ps_fr_project_offset(z - mid[2], dir[2], max_project, eps),
        p = mid + dir * offset,
        ex2d = [site[2][0], site[2][1]],
        line_dir = (norm(ex2d) <= eps) ? [1, 0] : v_norm(ex2d)
    )
    [[p[0], p[1]], line_dir, _ps_fr_project_was_capped(z - mid[2], dir[2], max_project, eps), site[0], site[8]];

/**
 * Function: Convert a circular list of projected boundary lines into loop vertices.
 * Params: lines (projected line records), eps (parallel tolerance)
 * Returns: 2D loop from intersections of adjacent lines
 */
function _ps_fr_projected_loop(lines, eps=1e-8) =
    let(n = len(lines))
    (n < 3) ? [] :
    [
        for (i = [0:1:n-1])
            let(hit = _ps_fr_line_intersection(lines[(i - 1 + n) % n], lines[i], eps))
            is_undef(hit) ? lines[i][0] : hit
    ];

/**
 * Function: Collect distinct loop ids from boundary-span sites preserving first-seen order.
 * Params: sites (boundary-span site records), i/acc (recursion state)
 * Returns: list of loop ids
 */
function _ps_fr_unique_loop_ids(sites, i=0, acc=[]) =
    (i >= len(sites)) ? acc :
    let(loop_idx = sites[i][7])
    _ps_fr_unique_loop_ids(
        sites,
        i + 1,
        _ps_list_contains(acc, loop_idx) ? acc : concat(acc, [loop_idx])
    );

/**
 * Function: Filter boundary-span sites to one boundary loop.
 * Params: sites (boundary-span site records), loop_idx (target loop id)
 * Returns: site records for that loop, in source boundary order
 */
function _ps_fr_sites_for_loop(sites, loop_idx) =
    [for (site = sites) if (site[7] == loop_idx) site];

/**
 * Function: Triangulate one projected cap loop into indexed polyhedron faces.
 * Params: loop2d (cap loop), offset (point-index offset), target_area_sign (desired triangle orientation), eps (tolerance)
 * Returns: list of triangle index faces
 */
function _ps_fr_cap_faces(loop2d, offset, target_area_sign, eps=1e-8) =
    [
        for (t = _ps_seg_triangulate_simple_poly_idx(loop2d, eps))
            let(
                area = _ps_fr_orient2(loop2d[t[0]], loop2d[t[1]], loop2d[t[2]]),
                oriented = (area * target_area_sign >= 0) ? t : [t[0], t[2], t[1]]
            )
            [for (idx = oriented) idx + offset]
    ];

/**
 * Function: Build side quad faces joining bottom and top loops.
 * Params: n (loop arity), loop_area_sign (bottom-loop signed area sign)
 * Returns: list of quad index faces
 */
function _ps_fr_side_faces(n, loop_area_sign) =
    [
        for (i = [0:1:n-1])
            let(j = (i + 1) % n)
            (loop_area_sign >= 0)
                ? [i, n + i, n + j, j]
                : [i, j, n + j, n + i]
    ];

/**
 * Function: Build one anti-interference shell record for one boundary loop.
 * Params: face_pts3d_local (source face loop), loop_sites (sites for one boundary loop), loop_idx (loop id), z0/z1 (target Z planes), input_sign (source loop winding sign), cell_winding_signs (per-cell winding signs), max_project (optional cap), eps (tolerance)
 * Returns: shell record `[points, faces, loop_idx, capped_count, bottom_loop2d, top_loop2d]`
 */
function _ps_fr_loop_shell(face_pts3d_local, loop_sites, loop_idx, z0, z1, input_sign, cell_winding_signs, max_project=undef, eps=1e-8) =
    let(
        lines0 = [for (site = loop_sites) _ps_fr_project_span_line(face_pts3d_local, site, z0, input_sign, cell_winding_signs, max_project, eps)],
        lines1 = [for (site = loop_sites) _ps_fr_project_span_line(face_pts3d_local, site, z1, input_sign, cell_winding_signs, max_project, eps)],
        loop0 = _ps_fr_projected_loop(lines0, eps),
        loop1 = _ps_fr_projected_loop(lines1, eps),
        n = min(len(loop0), len(loop1)),
        loop0n = [for (i = [0:1:n-1]) loop0[i]],
        loop1n = [for (i = [0:1:n-1]) loop1[i]],
        verts = concat(
            [for (p = loop0n) [p[0], p[1], z0]],
            [for (p = loop1n) [p[0], p[1], z1]]
        ),
        area = _ps_seg_poly_area2(loop0n),
        area_sign = (area >= 0) ? 1 : -1,
        faces = concat(
            _ps_fr_cap_faces(loop0n, 0, 1, eps),
            _ps_fr_cap_faces(loop1n, n, -1, eps),
            _ps_fr_side_faces(n, area_sign)
        ),
        capped_count = sum(concat(
            [for (line = lines0) line[2] ? 1 : 0],
            [for (line = lines1) line[2] ? 1 : 0]
        ))
    )
    [verts, faces, loop_idx, capped_count, loop0n, loop1n];

/**
 * Function: Build positive anti-interference shell meshes for one face.
 * Params: face_pts3d_local (current face loop in face-local 3D), face_idx (current face index), poly_faces_idx/poly_verts_local (full poly in current face-local coordinates), face_neighbors_idx/face_dihedrals (current face-edge metadata), z0/z1 (target local Z planes), mode (`"nonzero"`, `"evenodd"`, or `"all"`), max_project (optional projection-distance cap), eps (geometric tolerance)
 * Returns: list of shell records `[points, faces, loop_idx, capped_count, bottom_loop2d, top_loop2d]`
 * Limitations/Gotchas: emits one shell per filled boundary loop; holes/proxy punch-through volumes are intentionally outside this first primitive
 */
function ps_face_anti_interference_shells(
    face_pts3d_local,
    face_idx,
    poly_faces_idx,
    poly_verts_local,
    face_neighbors_idx,
    face_dihedrals,
    z0,
    z1,
    mode="nonzero",
    max_project=undef,
    eps=1e-8
) =
    let(
        _z0 = assert(!is_undef(z0), "ps_face_anti_interference_shells: z0 must be defined"),
        _z1 = assert(!is_undef(z1), "ps_face_anti_interference_shells: z1 must be defined"),
        arr = ps_face_arrangement(face_pts3d_local, eps),
        input_area = _ps_seg_poly_area2([for (p = face_pts3d_local) [p[0], p[1]]]),
        input_sign = (input_area >= 0) ? 1 : -1,
        cell_winding_signs = _ps_fr_cell_winding_signs(face_pts3d_local, arr[4], eps),
        sites = _ps_face_boundary_span_sites(
            face_pts3d_local,
            face_idx,
            poly_faces_idx,
            poly_verts_local,
            face_neighbors_idx,
            face_dihedrals,
            mode,
            eps
        ),
        loop_ids = _ps_fr_unique_loop_ids(sites)
    )
    [
        for (loop_idx = loop_ids)
            let(
                loop_sites = _ps_fr_sites_for_loop(sites, loop_idx),
                shell = _ps_fr_loop_shell(face_pts3d_local, loop_sites, loop_idx, z0, z1, input_sign, cell_winding_signs, max_project, eps)
            )
            if (len(shell[0]) >= 6 && len(shell[1]) >= 4)
                shell
    ];

/**
 * Module: Emit the current face's positive anti-interference volume.
 * Params: z0/z1 (target local Z planes), mode (`"nonzero"`, `"evenodd"`, or `"all"`), max_project (optional projection-distance cap), eps (geometric tolerance), convexity (OpenSCAD polyhedron convexity hint)
 * Returns: none; intended for use inside `place_on_faces(...)`, usually inside `intersection()`
 * Limitations/Gotchas: this is only the boundary-span volume primitive; it does not yet subtract or union proxy punch-through voids
 */
module ps_face_anti_interference_volume(z0, z1, mode="nonzero", max_project=undef, eps=1e-8, convexity=6) {
    assert(!is_undef($ps_face_pts3d_local), "ps_face_anti_interference_volume: requires place_on_faces context ($ps_face_pts3d_local)");
    assert(!is_undef($ps_face_idx), "ps_face_anti_interference_volume: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "ps_face_anti_interference_volume: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "ps_face_anti_interference_volume: requires place_on_faces context ($ps_poly_verts_local)");
    assert(!is_undef($ps_face_neighbors_idx), "ps_face_anti_interference_volume: requires place_on_faces context ($ps_face_neighbors_idx)");
    assert(!is_undef($ps_face_dihedrals), "ps_face_anti_interference_volume: requires place_on_faces context ($ps_face_dihedrals)");

    shells = ps_face_anti_interference_shells(
        $ps_face_pts3d_local,
        $ps_face_idx,
        $ps_poly_faces_idx,
        $ps_poly_verts_local,
        $ps_face_neighbors_idx,
        $ps_face_dihedrals,
        z0,
        z1,
        mode,
        max_project,
        eps
    );

    union() {
        for (shell = shells) {
            if (shell[3] > 0)
                echo(str("ps_face_anti_interference_volume: capped ", shell[3], " projection(s) on face ", $ps_face_idx, " loop ", shell[2]));

            polyhedron(points = shell[0], faces = shell[1], convexity = convexity);
        }
    }
}
