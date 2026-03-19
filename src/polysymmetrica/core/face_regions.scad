use <funcs.scad>
use <render.scad>
use <segments.scad>

function _ps_fr_safe_at(v, i, dflt) =
    (is_undef(i) || i < 0 || i >= len(v)) ? dflt : v[i];

function _ps_fr_line_intersect_2d(n1, d1, n2, d2) =
    let(det = n1[0] * n2[1] - n1[1] * n2[0])
    (abs(det) < 1e-9)
        ? [0, 0]
        : [
            (d1 * n2[1] - n1[1] * d2) / det,
            (n1[0] * d2 - d1 * n2[0]) / det
        ];

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
            _ps_fr_line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
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
                (_ps_fr_safe_at(s_kinds, k, "cut") == "parent")
                    ? ps_face_region_inset_at_z(_ps_fr_safe_at(face_diheds, _ps_fr_safe_at(s_edge_ids, k, 0), 180), z)
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

function ps_face_cut_profile2d_from_cutter_normal(z0, z1, inward_n, cutter_n3, cut_clearance=0, eps=1e-9) =
    let(
        z_min = min(z0, z1),
        z_max = max(z0, z1),
        spans_zero = (z_min < 0) && (z_max > 0),
        n_xy = [cutter_n3[0], cutter_n3[1]],
        n_use = (v_dot(n_xy, inward_n) > 0) ? -cutter_n3 : cutter_n3,
        b = v_norm([n_use[0], n_use[1], 1 + n_use[2]]),
        b_u = v_dot([b[0], b[1]], inward_n),
        b_z = b[2],
        slope_mag = (abs(b_u) <= eps) ? 0 : abs(-b_z / b_u)
    )
    spans_zero
        ? [
            [cut_clearance / 2 + slope_mag * abs(z_min), z_min],
            [cut_clearance / 2, 0],
            [cut_clearance / 2 + slope_mag * abs(z_max), z_max]
        ]
        : [
            [cut_clearance / 2 + slope_mag * abs(z_min), z_min],
            [cut_clearance / 2 + slope_mag * abs(z_max), z_max]
        ];

function ps_face_visible_cell_mask_loop(cell, cut_clearance=0) =
    let(
        s_pts = cell[0],
        s_kinds = cell[3],
        m = len(s_pts),
        insets = [
            for (k = [0:1:m-1])
                (_ps_fr_safe_at(s_kinds, k, "cut") == "parent") ? 0 : cut_clearance
                
        ]
    )
    _ps_fr_offset_pts2d_inset_edges(s_pts, insets);

module ps_face_visible_cell_volume(cell, z0, z1, cut_clearance=0, eps=1e-8) {
    loop = ps_face_visible_cell_mask_loop(cell, cut_clearance);
    if (len(loop) >= 3)
        translate([0, 0, z0])
            linear_extrude(height = z1 - z0)
                ps_polygon(points = loop, mode = "nonzero");
}

module ps_face_visible_cell_volume_ctx(z0, z1, cut_clearance=0, eps=1e-8) {
    assert(!is_undef($ps_vis_seg_pts2d), "ps_face_visible_cell_volume_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_pts2d)");
    assert(!is_undef($ps_vis_seg_edge_kinds), "ps_face_visible_cell_volume_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_edge_kinds)");
    cell = [$ps_vis_seg_pts2d, undef, undef, $ps_vis_seg_edge_kinds, is_undef($ps_vis_seg_cut_entry_ids) ? [] : $ps_vis_seg_cut_entry_ids];
    ps_face_visible_cell_volume(cell, z0, z1, cut_clearance, eps);
}

// Intersect arbitrary child geometry with one visible face cell. Optional
// cut-band subtraction is available for future segmented join relief work,
// but is disabled by default while the stable baseline remains plain cell
// clipping.
module ps_clip_to_visible_face_cell_ctx(z0, z1, cut_clearance=0, along_pad=0, mode="nonzero", eps=1e-8, apply_cut_bands=false, band_z0=undef, band_z1=undef, band_overcut=1e-3) {
    cell_cut_clearance = apply_cut_bands ? 0 : cut_clearance;
    cut_z0 = is_undef(band_z0) ? z0 : band_z0;
    cut_z1 = is_undef(band_z1) ? z1 : band_z1;
    cut_z0_eff = cut_z0 - band_overcut;
    cut_z1_eff = cut_z1 + band_overcut;
    difference() {
        intersection() {
            children();
            ps_face_visible_cell_volume_ctx(z0, z1, cell_cut_clearance, eps);
        }
        if (apply_cut_bands)
            ps_face_visible_segment_cut_bands_ctx(cut_z0_eff, cut_z1_eff, cut_clearance, along_pad, mode, eps, band_overcut);
    }
}

// Apply visible-cell clipping to arbitrary child geometry for all visible cells
// of the current face. This is the generic segmented-face consumer used by
// examples such as face_plate_visible(...).
module ps_clip_to_visible_face_segments_ctx(z0, z1, cut_clearance=0, along_pad=0, mode="nonzero", eps=1e-8, filter_parent=true, apply_cut_bands=false, band_z0=undef, band_z1=undef, band_overcut=1e-3) {
    union() {
        place_on_face_visible_segments(mode, eps, filter_parent) {
            ps_clip_to_visible_face_cell_ctx(z0, z1, cut_clearance, along_pad, mode, eps, apply_cut_bands, band_z0, band_z1, band_overcut)
                children();
        }
    }
}

function _ps_fr_cut_band_loop_u_overcut(seg2d, inward_n, u, along_pad=0, overcut=0) =
    let(
        a = seg2d[0],
        b = seg2d[1],
        e = v_norm(b - a),
        t = e * (along_pad + overcut),
        du0 = -inward_n * overcut,
        du1 = inward_n * (u + overcut)
    )
    [
        a - t + du0,
        b + t + du0,
        b + t + du1,
        a - t + du1
    ];

module _ps_fr_stack_quad_loops(loop_levels, eps=1e-8) {
    level_count = len(loop_levels);
    assert(level_count >= 2, "_ps_fr_stack_quad_loops: need at least 2 levels");
    points = [
        for (li = [0:1:level_count-1])
            let(loop2d = loop_levels[li][0], z = loop_levels[li][1])
            for (p = loop2d) [p[0], p[1], z]
    ];
    faces = concat(
        [[3, 2, 1, 0]],
        [[for (i = [0:1:3]) (level_count - 1) * 4 + i]],
        [
            for (li = [0:1:level_count-2])
                for (i = [0:1:3])
                    let(
                        lo = li * 4,
                        up = (li + 1) * 4,
                        i1 = (i + 1) % 4
                    )
                    each [
                        [up + i1, up + i, lo + i],
                        [up + i1, lo + i, lo + i1]
                    ]
        ]
    );
    polyhedron(points = points, faces = faces, convexity = 2);
}

module ps_face_cut_band_volume_profiled(seg2d, inward_n, profile2d, along_pad=0, eps=1e-8, overcut=0) {
    assert(len(profile2d) >= 2, "ps_face_cut_band_volume_profiled: need at least 2 profile points");
    loops = [
        for (pz = profile2d)
            [_ps_fr_cut_band_loop_u_overcut(seg2d, inward_n, pz[0], along_pad, overcut), pz[1]]
    ];
    _ps_fr_stack_quad_loops(loops, eps);
}

module ps_face_visible_segment_cut_bands_ctx(z0, z1, cut_clearance=0, along_pad=0, mode="nonzero", eps=1e-8, band_overcut=1e-3) {
    assert(!is_undef($ps_face_pts2d), "ps_face_visible_segment_cut_bands_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_idx), "ps_face_visible_segment_cut_bands_ctx: requires place_on_faces context ($ps_face_idx)");
    assert(!is_undef($ps_poly_faces_idx), "ps_face_visible_segment_cut_bands_ctx: requires place_on_faces context ($ps_poly_faces_idx)");
    assert(!is_undef($ps_poly_verts_local), "ps_face_visible_segment_cut_bands_ctx: requires place_on_faces context ($ps_poly_verts_local)");
    assert(!is_undef($ps_vis_seg_pts2d), "ps_face_visible_segment_cut_bands_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_pts2d)");
    assert(!is_undef($ps_vis_seg_edge_kinds), "ps_face_visible_segment_cut_bands_ctx: requires place_on_face_visible_segments context ($ps_vis_seg_edge_kinds)");

    entries = ps_face_geom_cut_entries($ps_face_pts2d, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local, eps, mode, true);
    probe = _ps_seg_cycle_probe_point($ps_vis_seg_pts2d, eps);
    m = len($ps_vis_seg_pts2d);
    for (k = [0:1:m-1]) {
        kind = _ps_fr_safe_at($ps_vis_seg_edge_kinds, k, "parent");
        cid = _ps_fr_safe_at($ps_vis_seg_cut_entry_ids, k, undef);
        if (kind == "cut" && !is_undef(cid) && cid >= 0 && cid < len(entries)) {
            a = $ps_vis_seg_pts2d[k];
            b = $ps_vis_seg_pts2d[(k + 1) % m];
            e = v_norm(b - a);
            perp0 = v_norm([-e[1], e[0]]);
            mid = (a + b) / 2;
            inward_n = (v_dot(perp0, probe - mid) < 0) ? -perp0 : perp0;
            cutter_f = $ps_poly_faces_idx[entries[cid][1]];
            cutter_n3 = ps_face_frame_normal($ps_poly_verts_local, cutter_f);
            ps_face_cut_band_volume_profiled([a, b], inward_n, ps_face_cut_profile2d_from_cutter_normal(z0, z1, inward_n, cutter_n3, cut_clearance, eps), along_pad, eps, band_overcut);
        }
    }
}
