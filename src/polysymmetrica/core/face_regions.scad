use <funcs.scad>
use <render.scad>
use <segments.scad>

function _ps_fr_loop_from_adjacent_lines(lines) =
    let(m = len(lines))
    [
        for (k = [0:1:m-1])
            let(
                l0 = lines[(k - 1 + m) % m],
                l1 = lines[k]
            )
            _ps_line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
    ];

function _ps_fr_line_keep(line, p, eps=1e-9) =
    let(
        n = line[0],
        d = line[1],
        keep_ge = line[2],
        val = v_dot(n, p) - d
    )
    keep_ge ? (val >= -eps) : (val <= eps);

function _ps_fr_line_seg_intersection(line, a, b, eps=1e-9) =
    let(
        n = line[0],
        d = line[1],
        da = v_dot(n, a) - d,
        db = v_dot(n, b) - d,
        denom = da - db,
        t = (abs(denom) <= eps) ? 0 : da / denom
    )
    a + (b - a) * t;

function _ps_fr_clip_poly_line(poly, line, eps=1e-9) =
    let(
        m = len(poly)
    )
    (m == 0) ? [] :
    [
        for (i = [0:1:m-1])
            let(
                a = poly[i],
                b = poly[(i + 1) % m],
                a_keep = _ps_fr_line_keep(line, a, eps),
                b_keep = _ps_fr_line_keep(line, b, eps),
                x = _ps_fr_line_seg_intersection(line, a, b, eps)
            )
            each (
                a_keep && b_keep ? [b] :
                a_keep && !b_keep ? [x] :
                !a_keep && b_keep ? [x, b] :
                []
            )
    ];

function _ps_fr_clip_poly_lines(poly, lines, eps=1e-9, i=0) =
    (i >= len(lines) || len(poly) == 0)
        ? poly
        : _ps_fr_clip_poly_lines(_ps_fr_clip_poly_line(poly, lines[i], eps), lines, eps, i + 1);

function _ps_fr_lines_seed_square(lines, probe, eps=1e-9) =
    let(
        ds = [for (l = lines) abs(l[1])],
        max_d = (len(ds) == 0) ? 1 : max(ds),
        r = 4 * (norm(probe) + max_d + 1 + eps)
    )
    [
        probe + [-r, -r],
        probe + [ r, -r],
        probe + [ r,  r],
        probe + [-r,  r]
    ];

function _ps_fr_loop_from_halfplanes(lines, probe, eps=1e-9) =
    _ps_fr_clip_poly_lines(_ps_fr_lines_seed_square(lines, probe, eps), lines, eps);

function _ps_fr_offset_pts2d_inset_edges(pts, insets) =
    let(
        n = len(pts),
        centroid = v_scale([
            sum([for (p = pts) p[0]]),
            sum([for (p = pts) p[1]])
        ], 1 / n),
        lines = [
            for (k = [0:1:n-1])
                let(
                    p0 = pts[k],
                    p1 = pts[(k+1)%n],
                    e = v_norm(p1 - p0),
                    perp0 = v_norm([-e[1], e[0]]),
                    mid = (p0 + p1) / 2,
                    to_center = centroid - mid,
                    perp = (v_dot(perp0, to_center) < 0) ? -perp0 : perp0,
                    d = v_dot(perp, p0) + insets[k]
                )
                [perp, d]
        ]
    )
    [
        for (k = [0:1:n-1])
            let(
                l0 = lines[(k-1+n)%n],
                l1 = lines[k]
            )
            _ps_line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
    ];

function ps_face_region_inset_at_z(dihed, z) =
    let(angle = (180 - dihed) / 2)
    -z * tan(angle);

function _ps_fr_loop_from_seg_at_z(seg, face_diheds, z) =
    let(
        s_pts = seg[0],
        s_edge_ids = seg[2],
        s_kinds = seg[3],
        m = len(s_pts),
        insets = [
            for (k = [0:1:m-1])
                (_ps_seg_safe_at(s_kinds, k, "cut") == "parent")
                    ? ps_face_region_inset_at_z(_ps_seg_safe_at(face_diheds, _ps_seg_safe_at(s_edge_ids, k, 0), 180), z)
                    : 0
        ]
    )
    _ps_fr_offset_pts2d_inset_edges(s_pts, insets);

function ps_face_region_loops_at_z(face_pts2d, face_diheds, z, mode="nonzero", eps=1e-8) =
    let(
        segs = ps_face_segments([for (p = face_pts2d) [p[0], p[1], 0]], mode, eps)
    )
    [
        for (seg = segs)
            _ps_fr_loop_from_seg_at_z(seg, face_diheds, z)
    ];

module _ps_fr_loft_loops(loop_top2d, loop_bottom2d, z_top, z_bottom, eps=1e-8, cap_top=true, cap_bottom=true) {
    m = len(loop_top2d);
    assert(m >= 3, "_ps_fr_loft_loops: need at least 3 points");
    assert(len(loop_bottom2d) == m, "_ps_fr_loft_loops: loop sizes must match");

    top = [for (p = loop_top2d) [p[0], p[1], z_top]];
    bottom = [for (p = loop_bottom2d) [p[0], p[1], z_bottom]];
    top_tris = _ps_seg_triangulate_simple_poly_idx(loop_top2d, eps);
    bottom_tris = _ps_seg_triangulate_simple_poly_idx(loop_bottom2d, eps);

    polyhedron(
        points = concat(top, bottom),
        faces = concat(
            cap_top ? [for (t = top_tris) [t[0], t[1], t[2]]] : [],
            cap_bottom ? [for (t = bottom_tris) [2 * m - 1 - t[0], 2 * m - 1 - t[1], 2 * m - 1 - t[2]]] : [],
            [for (i = [0:1:m-1]) [(i+1)%m, i, m + i]],
            [for (i = [0:1:m-1]) [(i+1)%m, m + i, m + (i+1)%m]]
        ),
        convexity = 2
    );
}

// Volume around a face defined by the face plane at z=0 and the real face
// side planes split by dihedral/2. Positive z is outward from the face.
module ps_face_region_volume(face_pts2d, face_diheds, z0, z1, mode="nonzero", eps=1e-8) {
    z_min = min(z0, z1);
    z_max = max(z0, z1);
    segs = ps_face_segments([for (p = face_pts2d) [p[0], p[1], 0]], mode, eps);
    for (seg = segs) {
        loop0 = _ps_fr_loop_from_seg_at_z(seg, face_diheds, z_min);
        loop1 = _ps_fr_loop_from_seg_at_z(seg, face_diheds, z_max);
        if (len(loop0) >= 3 && len(loop1) >= 3 && len(loop0) == len(loop1))
            _ps_fr_loft_loops(loop1, loop0, z_max, z_min, eps);
    }
}

module ps_face_region_volume_ctx(z0, z1, mode="nonzero", eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_region_volume_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_dihedrals), "ps_face_region_volume_ctx: requires place_on_faces context ($ps_face_dihedrals)");
    ps_face_region_volume($ps_face_pts2d, $ps_face_dihedrals, z0, z1, mode, eps);
}

// Intersect arbitrary child geometry with the local admissible face region.
module ps_clip_to_face_region(face_pts2d, face_diheds, z0, z1, mode="nonzero", eps=1e-8) {
    intersection() {
        children();
        ps_face_region_volume(face_pts2d, face_diheds, z0, z1, mode, eps);
    }
}

// Context wrapper for ps_clip_to_face_region(...) inside place_on_faces(...).
module ps_clip_to_face_region_ctx(z0, z1, mode="nonzero", eps=1e-8) {
    intersection() {
        children();
        ps_face_region_volume_ctx(z0, z1, mode, eps);
    }
}

function ps_face_cut_join_dihed(cut_dihed) = 360 - cut_dihed;

function ps_face_cut_profile2d_from_cutter_normal(z0, z1, inward_n, cutter_n3, cut_clearance=0, cut_dihed=undef, eps=1e-9) =
    let(
        z_min = min(z0, z1),
        z_max = max(z0, z1),
        n_xy = [cutter_n3[0], cutter_n3[1]],
        n_use = (v_dot(n_xy, inward_n) > 0) ? -cutter_n3 : cutter_n3,
        b = v_norm([n_use[0], n_use[1], 1 + n_use[2]]),
        b_use = (v_dot([b[0], b[1]], inward_n) < 0) ? -b : b,
        b_u = v_dot([b_use[0], b_use[1]], inward_n),
        b_z = b_use[2],
        slope_geom = (abs(b_u) <= eps) ? 0 : (-b_z / b_u),
        slope_mag = is_undef(cut_dihed) ? abs(slope_geom) : tan((180 - cut_dihed) / 2),
        slope_signed = ((slope_geom < 0) ? -1 : 1) * slope_mag,
        u0_raw = cut_clearance / 2 + slope_signed * z_min,
        u1_raw = cut_clearance / 2 + slope_signed * z_max,
        shift = max(0, cut_clearance / 2 - min(u0_raw, u1_raw))
    )
    [
        [u0_raw + shift, z_min],
        [u1_raw + shift, z_max]
    ];

function ps_face_visible_cell_mask_loop(cell, cut_clearance=0) =
    let(
        s_pts = cell[0],
        s_kinds = cell[3],
        m = len(s_pts),
        insets = [
            for (k = [0:1:m-1])
                (_ps_seg_safe_at(s_kinds, k, "cut") == "parent") ? 0 : cut_clearance
                
        ]
    )
    _ps_fr_offset_pts2d_inset_edges(s_pts, insets);

function _ps_fr_profile_u_at_z(profile2d, z, eps=1e-9) =
    (len(profile2d) == 0) ? 0 :
    (len(profile2d) == 1) ? profile2d[0][0] :
    let(
        zs = [for (pz = profile2d) pz[1]],
        z_min = min(zs),
        z_max = max(zs)
    )
    (z <= z_min + eps)
        ? profile2d[0][0]
        : (z >= z_max - eps)
            ? profile2d[len(profile2d) - 1][0]
            : let(
                idxs = [
                    for (i = [0:1:len(profile2d)-2])
                        if (profile2d[i][1] <= z + eps && z <= profile2d[i + 1][1] + eps) i
                ],
                i = (len(idxs) == 0) ? 0 : idxs[0],
                p0 = profile2d[i],
                p1 = profile2d[i + 1],
                dz = p1[1] - p0[1],
                t = (abs(dz) <= eps) ? 0 : (z - p0[1]) / dz
            )
            p0[0] + (p1[0] - p0[0]) * t;

function _ps_fr_edge_plane_from_profile(pts, edge_idx, probe2d, profile2d, z_mid=0, eps=1e-9) =
    let(
        m = len(pts),
        p0 = pts[edge_idx],
        inward_n = _ps_fr_cell_edge_inward_n(pts, edge_idx, probe2d, eps),
        u0 = profile2d[0][0],
        z0 = profile2d[0][1],
        u1 = profile2d[len(profile2d) - 1][0],
        z1 = profile2d[len(profile2d) - 1][1],
        dz = z1 - z0,
        slope = (abs(dz) <= eps) ? 0 : (u1 - u0) / dz,
        intercept = u0 - slope * z0,
        n0 = [inward_n[0], inward_n[1], -slope],
        d0 = v_dot(inward_n, p0) + intercept,
        probe3 = [probe2d[0], probe2d[1], z_mid],
        keep_ok = v_dot(n0, probe3) >= d0 - 10 * eps
    )
    keep_ok ? [n0, d0] : [-n0, -d0];

function ps_face_visible_cell_region_planes(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8) =
    let(
        z_min = min(z0, z1),
        z_max = max(z0, z1),
        cut_lo = is_undef(band_z0) ? z_min : band_z0,
        cut_hi = is_undef(band_z1) ? z_max : band_z1,
        pts = cell[0],
        m = len(pts),
        probe = _ps_seg_cycle_probe_point(pts, eps),
        z_mid = (z_min + z_max) / 2
    )
    concat(
        [
            [[0, 0, 1], z_min],
            [[0, 0, -1], -z_max]
        ],
        [
            for (k = [0:1:m-1])
                let(
                    profile = [
                        [_ps_fr_visible_cell_edge_u_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, k, z_min, cut_lo, cut_hi, cut_clearance, eps), z_min],
                        [_ps_fr_visible_cell_edge_u_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, k, z_max, cut_lo, cut_hi, cut_clearance, eps), z_max]
                    ]
                )
                _ps_fr_edge_plane_from_profile(pts, k, probe, profile, z_mid, eps)
        ]
    );

function _ps_fr_plane_slice_line(plane, probe2d, z, eps=1e-9) =
    let(
        n3 = plane[0],
        n2 = [n3[0], n3[1]],
        d2 = plane[1] - n3[2] * z
    )
    (norm(n2) <= eps)
        ? undef
        : [n2, d2, v_dot(n2, probe2d) >= d2];

function ps_face_visible_cell_loop_at_z_from_region_planes(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8) =
    let(
        z_min = min(z0, z1),
        z_max = max(z0, z1),
        probe = _ps_seg_cycle_probe_point(cell[0], eps),
        planes = ps_face_visible_cell_region_planes(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0, band_z1, cut_clearance, eps),
        lines = [
            for (pi = [2:1:len(planes)-1])
                let(line = _ps_fr_plane_slice_line(planes[pi], probe, z, eps))
                if (!is_undef(line)) line
        ]
    )
    (z < z_min - eps || z > z_max + eps || len(lines) < 3)
        ? []
        : _ps_fr_loop_from_halfplanes(lines, probe, eps);

function _ps_fr_atom_edge_meta(cell, ia, ib, eps=1e-9) =
    let(
        pts = cell[0],
        edge_ids = cell[2],
        kinds = cell[3],
        cut_ids = cell[4],
        n = len(pts),
        fwd = [
            for (k = [0:1:n-1])
                if (k == ia && ((ia + 1) % n) == ib)
                    [_ps_seg_safe_at(edge_ids, k, undef), _ps_seg_safe_at(kinds, k, "inner"), _ps_seg_safe_at(cut_ids, k, undef)]
        ],
        rev = [
            for (k = [0:1:n-1])
                if (k == ib && ((ib + 1) % n) == ia)
                    [_ps_seg_safe_at(edge_ids, k, undef), _ps_seg_safe_at(kinds, k, "inner"), _ps_seg_safe_at(cut_ids, k, undef)]
        ]
    )
    (len(fwd) > 0) ? fwd[0] :
    (len(rev) > 0) ? rev[0] :
    [undef, "inner", undef];

function _ps_fr_visible_cell_atoms(cell, eps=1e-9) =
    let(
        pts = cell[0],
        tris = _ps_seg_triangulate_simple_poly_idx(pts, eps)
    )
    _ps_poly_is_convex2(pts, eps)
        ? [cell]
        : [
            for (tri = tris)
                let(
                    ia = tri[0],
                    ib = tri[1],
                    ic = tri[2],
                    m01 = _ps_fr_atom_edge_meta(cell, ia, ib, eps),
                    m12 = _ps_fr_atom_edge_meta(cell, ib, ic, eps),
                    m20 = _ps_fr_atom_edge_meta(cell, ic, ia, eps)
                )
                [
                    [pts[ia], pts[ib], pts[ic]],
                    undef,
                    [m01[0], m12[0], m20[0]],
                    [m01[1], m12[1], m20[1]],
                    [m01[2], m12[2], m20[2]]
                ]
        ];

function _ps_fr_cell_edge_inward_n(cell_pts2d, edge_idx, probe, eps=1e-9) =
    let(
        m = len(cell_pts2d),
        a = cell_pts2d[edge_idx],
        b = cell_pts2d[(edge_idx + 1) % m],
        e = v_norm(b - a),
        perp0 = v_norm([-e[1], e[0]]),
        mid = (a + b) / 2
    )
    (v_dot(perp0, probe - mid) < 0) ? -perp0 : perp0;

function _ps_fr_visible_cell_edge_u_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, edge_idx, z, band_z0, band_z1, cut_clearance=0, eps=1e-9) =
    let(
        kinds = cell[3],
        edge_ids = cell[2],
        cut_ids = cell[4],
        kind = _ps_seg_safe_at(kinds, edge_idx, "cut")
    )
    (kind == "inner")
        ? 0
        : (kind == "parent")
        ? ps_face_region_inset_at_z(_ps_seg_safe_at(face_diheds, _ps_seg_safe_at(edge_ids, edge_idx, 0), 180), z)
            : let(
                cid = _ps_seg_safe_at(cut_ids, edge_idx, undef),
                probe = _ps_seg_cycle_probe_point(cell[0], eps),
                inward_n = _ps_fr_cell_edge_inward_n(cell[0], edge_idx, probe, eps),
                valid_cid = !is_undef(cid) && cid >= 0 && cid < len(cut_entries)
            )
        !valid_cid
            ? cut_clearance
            : let(
                cutter_f = poly_faces_idx[cut_entries[cid][1]],
                cutter_n3 = ps_face_frame_normal(poly_verts_local, cutter_f),
                profile = ps_face_cut_profile2d_from_cutter_normal(band_z0, band_z1, inward_n, cutter_n3, cut_clearance, cut_entries[cid][2], eps)
            )
            _ps_fr_profile_u_at_z(profile, z, eps);

function _ps_fr_visible_cell_lines_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, band_z0, band_z1, cut_clearance=0, eps=1e-9) =
    let(
        pts = cell[0],
        m = len(pts),
        probe = _ps_seg_cycle_probe_point(pts, eps)
    )
    [
        for (k = [0:1:m-1])
            let(
                p0 = pts[k],
                inward_n = _ps_fr_cell_edge_inward_n(pts, k, probe, eps),
                u = _ps_fr_visible_cell_edge_u_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, k, z, band_z0, band_z1, cut_clearance, eps),
                d = v_dot(inward_n, p0) + u,
                keep_ge = v_dot(inward_n, probe) >= d
            )
            [inward_n, d, keep_ge]
    ];

function _ps_fr_visible_cell_loop_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, band_z0, band_z1, cut_clearance=0, eps=1e-9) =
    let(
        lines = _ps_fr_visible_cell_lines_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, band_z0, band_z1, cut_clearance, eps),
        seed = _ps_fr_loop_from_adjacent_lines(lines)
    )
    seed;

function ps_face_visible_cell_loop_at_z_clipped(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, band_z0, band_z1, cut_clearance=0, eps=1e-9) =
    let(
        probe = _ps_seg_cycle_probe_point(cell[0], eps),
        lines = _ps_fr_visible_cell_lines_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, band_z0, band_z1, cut_clearance, eps)
    )
    _ps_fr_loop_from_halfplanes(lines, probe, eps);

module ps_face_visible_cell_volume(cell, z0, z1, cut_clearance=0, eps=1e-8) {
    loop = ps_face_visible_cell_mask_loop(cell, cut_clearance);
    if (len(loop) >= 3)
        translate([0, 0, z0])
            linear_extrude(height = z1 - z0)
                ps_polygon(points = loop, mode = "nonzero");
}

module _ps_fr_visible_convex_atom_region_volume(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8) {
    z_min = min(z0, z1);
    z_max = max(z0, z1);
    cut_lo = is_undef(band_z0) ? z_min : band_z0;
    cut_hi = is_undef(band_z1) ? z_max : band_z1;
    levels0 = [z_min, cut_lo, cut_hi, z_max];
    levels = [
        for (i = [0:1:len(levels0)-1])
            if ((i == 0) || abs(levels0[i] - levels0[i - 1]) > eps) levels0[i]
    ];
    loops = [
        for (z = levels)
            [_ps_fr_visible_cell_loop_at_z(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, cut_lo, cut_hi, cut_clearance, eps), z]
    ];
    if (min([for (lz = loops) len(lz[0])]) >= 3)
        _ps_fr_stack_same_arity_loops(loops, eps);
}

module ps_face_visible_cell_region_volume(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8) {
    atoms = _ps_fr_visible_cell_atoms(cell, eps);
    union() {
        for (atom = atoms)
            _ps_fr_visible_convex_atom_region_volume(atom, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0, band_z1, cut_clearance, eps);
    }
}

module ps_face_visible_cell_region_volume_ctx(z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, mode="nonzero", eps=1e-8, filter_parent=true) {
    assert(!is_undef($ps_face_pts2d), "ps_face_visible_cell_region_volume_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_idx), "ps_face_visible_cell_region_volume_ctx: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "ps_face_visible_cell_region_volume_ctx: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "ps_face_visible_cell_region_volume_ctx: requires place_on_faces context ($ps_poly_verts_local)");
    assert(!is_undef($ps_face_dihedrals), "ps_face_visible_cell_region_volume_ctx: requires place_on_faces context ($ps_face_dihedrals)");
    assert(!is_undef($ps_vis_seg_pts2d), "ps_face_visible_cell_region_volume_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_pts2d)");
    assert(!is_undef($ps_vis_seg_edge_kinds), "ps_face_visible_cell_region_volume_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_edge_kinds)");
    cell = [$ps_vis_seg_pts2d, undef, $ps_vis_seg_edge_ids, $ps_vis_seg_edge_kinds, is_undef($ps_vis_seg_cut_entry_ids) ? [] : $ps_vis_seg_cut_entry_ids];
    cut_entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, mode, filter_parent);
    ps_face_visible_cell_region_volume(cell, $ps_face_dihedrals, cut_entries, $ps_poly_faces_idx, $ps_poly_verts_local, z0, z1, band_z0, band_z1, cut_clearance, eps);
}

module ps_face_visible_cell_volume_ctx(z0, z1, cut_clearance=0, eps=1e-8) {
    assert(!is_undef($ps_vis_seg_pts2d), "ps_face_visible_cell_volume_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_pts2d)");
    assert(!is_undef($ps_vis_seg_edge_kinds), "ps_face_visible_cell_volume_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_edge_kinds)");
    cell = [$ps_vis_seg_pts2d, undef, undef, $ps_vis_seg_edge_kinds, is_undef($ps_vis_seg_cut_entry_ids) ? [] : $ps_vis_seg_cut_entry_ids];
    ps_face_visible_cell_volume(cell, z0, z1, cut_clearance, eps);
}

// Intersect arbitrary child geometry with one visible face cell.
// Baseline mode clips to the simple visible-cell prism; `apply_cut_bands=true`
// switches to the direct visible-cell region built from all edge side-lines
// together.
module ps_clip_to_visible_face_cell_ctx(z0, z1, cut_clearance=0, mode="nonzero", eps=1e-8, apply_cut_bands=false, band_z0=undef, band_z1=undef) {
    if (apply_cut_bands) {
        intersection() {
            children();
            ps_face_visible_cell_region_volume_ctx(z0, z1, band_z0, band_z1, cut_clearance, mode, eps);
        }
    } else {
        difference() {
            intersection() {
                children();
                ps_face_visible_cell_volume_ctx(z0, z1, cut_clearance, eps);
            }
        }
    }
}

// Apply visible-cell clipping to arbitrary child geometry for all visible cells
// of the current face. This is the generic segmented-face consumer used by
// examples such as face_plate_visible(...).
module ps_clip_to_visible_face_segments_ctx(z0, z1, cut_clearance=0, mode="nonzero", eps=1e-8, filter_parent=true, apply_cut_bands=false, band_z0=undef, band_z1=undef) {
    union() {
        place_on_face_visible_segments(mode, eps, filter_parent) {
            ps_clip_to_visible_face_cell_ctx(z0, z1, cut_clearance, mode, eps, apply_cut_bands, band_z0, band_z1)
                children();
        }
    }
}

module _ps_fr_stack_same_arity_loops(loop_levels, eps=1e-8) {
    level_count = len(loop_levels);
    assert(level_count >= 2, "_ps_fr_stack_same_arity_loops: need at least 2 levels");
    m = len(loop_levels[0][0]);
    assert(m >= 3, "_ps_fr_stack_same_arity_loops: need at least 3 points per loop");
    for (li = [1:1:level_count-1])
        assert(len(loop_levels[li][0]) == m, "_ps_fr_stack_same_arity_loops: loop sizes must match");

    points = [
        for (li = [0:1:level_count-1])
            let(loop2d = loop_levels[li][0], z = loop_levels[li][1])
            for (p = loop2d) [p[0], p[1], z]
    ];
    bottom_tris = _ps_seg_triangulate_simple_poly_idx(loop_levels[0][0], eps);
    top_tris = _ps_seg_triangulate_simple_poly_idx(loop_levels[level_count - 1][0], eps);
    faces = concat(
        [for (t = bottom_tris) [m - 1 - t[0], m - 1 - t[1], m - 1 - t[2]]],
        [for (t = top_tris) [(level_count - 1) * m + t[0], (level_count - 1) * m + t[1], (level_count - 1) * m + t[2]]],
        [
            for (li = [0:1:level_count-2])
                for (i = [0:1:m-1])
                    let(
                        lo = li * m,
                        up = (li + 1) * m,
                        i1 = (i + 1) % m
                    )
                    each [
                        [up + i1, up + i, lo + i],
                        [up + i1, lo + i, lo + i1]
                    ]
        ]
    );
    polyhedron(points = points, faces = faces, convexity = 2);
}
