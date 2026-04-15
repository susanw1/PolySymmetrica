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

function _ps_fr_probe_point(poly, eps=1e-9) =
    let(
        probe = _ps_seg_cycle_probe_point(poly, eps),
        n = len(poly),
        centroid = (n == 0) ? [0, 0] : v_scale([
            sum([for (p = poly) p[0]]),
            sum([for (p = poly) p[1]])
        ], 1 / n)
    )
    is_undef(probe) ? centroid : probe;

function _ps_fr_norm2(v, eps=1e-9) =
    let(L = norm(v))
    (L <= eps) ? [0, 0] : v / L;

function _ps_fr_face_span2(face_pts2d) =
    let(
        xs = [for (p = face_pts2d) p[0]],
        ys = [for (p = face_pts2d) p[1]]
    )
    norm([max(xs) - min(xs), max(ys) - min(ys)]);

function _ps_fr_line_keep(line, p, eps=1e-9) =
    (is_undef(p) || is_undef(line) || len(line) < 3 || is_undef(line[0]) || len(line[0]) != 2 || is_undef(line[1]) || is_undef(line[2])) ? false :
    let(
        n = line[0],
        d = line[1],
        keep_ge = line[2],
        val = v_dot(n, p) - d
    )
    keep_ge ? (val >= -eps) : (val <= eps);

function _ps_fr_line_seg_intersection(line, a, b, eps=1e-9) =
    (is_undef(a) || is_undef(b) || is_undef(line) || len(line) < 3 || is_undef(line[0]) || len(line[0]) != 2 || is_undef(line[1])) ? undef :
    let(
        n = line[0],
        d = line[1],
        da = v_dot(n, a) - d,
        db = v_dot(n, b) - d,
        denom = da - db,
        t = (abs(denom) <= eps) ? 0 : da / denom
    )
    (abs(denom) <= eps) ? undef : a + (b - a) * t;

function _ps_fr_clip_poly_line(poly, line, eps=1e-9) =
    let(
        clean = [for (p = poly) if (!is_undef(p) && len(p) == 2) p],
        m = len(clean)
    )
    (m == 0) ? [] :
    [
        for (i = [0:1:m-1])
            let(
                a = clean[i],
                b = clean[(i + 1) % m],
                a_keep = _ps_fr_line_keep(line, a, eps),
                b_keep = _ps_fr_line_keep(line, b, eps),
                x = _ps_fr_line_seg_intersection(line, a, b, eps)
            )
            each (
                a_keep && b_keep ? [b] :
                a_keep && !b_keep ? (is_undef(x) ? [] : [x]) :
                !a_keep && b_keep ? (is_undef(x) ? [b] : [x, b]) :
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
                    e = _ps_fr_norm2(p1 - p0),
                    perp0 = _ps_fr_norm2([-e[1], e[0]]),
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

module _ps_fr_raw_filled_face_volume(face_pts3d_local, z0, z1, eps=1e-8) {
    z_min = min(z0, z1);
    z_max = max(z0, z1);

    translate([0, 0, z_min])
        linear_extrude(height = z_max - z_min)
            union() {
                for (cell = ps_face_filled_cells(face_pts3d_local, eps))
                    ps_polygon(points = cell[0]);
            }
}

module _ps_fr_atom_prism(atom, z0, z1) {
    z_min = min(z0, z1);
    z_max = max(z0, z1);

    translate([0, 0, z_min])
        linear_extrude(height = z_max - z_min)
            ps_polygon(points = atom[0]);
}

module _ps_fr_filled_cell_prism(cell, z0, z1) {
    z_min = min(z0, z1);
    z_max = max(z0, z1);

    translate([0, 0, z_min])
        linear_extrude(height = z_max - z_min)
            ps_polygon(points = cell[0]);
}

module _ps_fr_interference_boundary_cutter(face_span, z0, z1, inward_sign=1, eps=1e-8) {
    along_far = max(1, 4 * face_span, 4 * abs(z1 - z0));
    outward_far = max(1, 2 * face_span, 4 * abs(z1 - z0));
    z_far = max(1, 4 * face_span, 4 * abs(z1 - z0));
    cut_eps = max(eps, 1e-6);

    // Keep-side depends on the segment winding.
    // If inward_sign > 0 then local +Y is the keep side; otherwise local -Y is.
    // So the cutter must remove the opposite half-space.
    translate([0, -inward_sign * (outward_far - cut_eps) / 2, 0])
        cube([2 * along_far + 2 * cut_eps, outward_far + cut_eps, z_far + 2 * cut_eps], center = true);
}

function _ps_fr_atom_edge_inward_sign(atom_pts2d, seg2d, eps=1e-9) =
    let(
        mid = _ps_seg_boundary_midpoint(seg2d),
        left_n = _ps_seg_boundary_left_normal(seg2d, eps),
        probe = _ps_seg_cycle_probe_point(atom_pts2d, eps)
    )
    (v_dot(probe - mid, left_n) >= 0) ? 1 : -1;

module _ps_fr_atom_boundary_interference_cutters(atom, face_diheds, face_span, z0, z1, eps=1e-8) {
    pts2d = atom[0];
    source_edge_ids = atom[2];
    edge_kinds = atom[3];
    n = len(pts2d);

    for (k = [0:1:n-1]) {
        edge_kind = _ps_seg_safe_at(edge_kinds, k, "inner");

        if (edge_kind == "boundary") {
            p0 = pts2d[k];
            p1 = pts2d[(k + 1) % n];
            seg2d = [p0, p1];
            seg_len = norm(p1 - p0);
            ex2 = (seg_len <= eps) ? [1, 0] : (p1 - p0) / seg_len;
            ey2 = _ps_seg_boundary_left_normal(seg2d, eps);
            inward_sign = _ps_fr_atom_edge_inward_sign(pts2d, seg2d, eps);
            inward2 = ey2 * inward_sign;
            source_edge_idx = _ps_seg_safe_at(source_edge_ids, k, undef);
            dihed = _ps_seg_safe_at(face_diheds, source_edge_idx, 180);
            mid = _ps_seg_boundary_midpoint(seg2d);
            ex3 = [ex2[0], ex2[1], 0];
            ez3 = _ps_seg_boundary_edge_ez3(inward2, dihed, eps);
            ey3 = v_norm(v_cross(ez3, ex3));

            multmatrix(ps_frame_matrix([mid[0], mid[1], 0], ex3, ey3, ez3))
                _ps_fr_interference_boundary_cutter(face_span, z0, z1, inward_sign, eps);
        }
    }
}

module ps_face_interference_atom_boundary_cutter_ctx(atom_idx, edge_idx, z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_interference_atom_boundary_cutter_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_pts3d_local), "ps_face_interference_atom_boundary_cutter_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);
    face_span = _ps_fr_face_span2($ps_face_pts2d);
    assert(atom_idx >= 0 && atom_idx < len(atoms), str("ps_face_interference_atom_boundary_cutter_ctx: atom_idx out of range ", atom_idx));

    atom = atoms[atom_idx];
    pts2d = atom[0];
    source_edge_ids = atom[2];
    edge_kinds = atom[3];
    n = len(pts2d);

    assert(edge_idx >= 0 && edge_idx < n, str("ps_face_interference_atom_boundary_cutter_ctx: edge_idx out of range ", edge_idx));
    assert(_ps_seg_safe_at(edge_kinds, edge_idx, "inner") == "boundary", str("ps_face_interference_atom_boundary_cutter_ctx: edge ", edge_idx, " is not a boundary edge"));

    p0 = pts2d[edge_idx];
    p1 = pts2d[(edge_idx + 1) % n];
    seg2d = [p0, p1];
    seg_len = norm(p1 - p0);
    ex2 = (seg_len <= eps) ? [1, 0] : (p1 - p0) / seg_len;
    ey2 = _ps_seg_boundary_left_normal(seg2d, eps);
    inward_sign = _ps_fr_atom_edge_inward_sign(pts2d, seg2d, eps);
    inward2 = ey2 * inward_sign;
    source_edge_idx = _ps_seg_safe_at(source_edge_ids, edge_idx, undef);
    dihed = _ps_seg_safe_at($ps_face_dihedrals, source_edge_idx, 180);
    mid = _ps_seg_boundary_midpoint(seg2d);
    ex3 = [ex2[0], ex2[1], 0];
    ez3 = _ps_seg_boundary_edge_ez3(inward2, dihed, eps);
    ey3 = v_norm(v_cross(ez3, ex3));

    multmatrix(ps_frame_matrix([mid[0], mid[1], 0], ex3, ey3, ez3))
        _ps_fr_interference_boundary_cutter(face_span, z0, z1, inward_sign, eps);
}

// Default interference cutter for the current true filled-boundary segment,
// assuming we are already in that segment's edge frame.
module ps_face_boundary_seg_default_interference_cutter_local_ctx(z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_boundary_seg_default_interference_cutter_local_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_boundary_seg_inward_is_positive_ey), "ps_face_boundary_seg_default_interference_cutter_local_ctx: requires place_on_face_filled_boundary_* context");

    inward_sign = $ps_face_boundary_seg_inward_is_positive_ey ? 1 : -1;
    face_span = _ps_fr_face_span2($ps_face_pts2d);

    _ps_fr_interference_boundary_cutter(face_span, z0, z1, -inward_sign, eps);
}

// Boundary cutter orientation for material bodies derived from the face region
// itself, for example crossing lobes used in phase-2 subtraction.
module ps_face_boundary_seg_body_interference_cutter_local_ctx(z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_boundary_seg_body_interference_cutter_local_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_boundary_seg_inward_is_positive_ey), "ps_face_boundary_seg_body_interference_cutter_local_ctx: requires place_on_face_filled_boundary_* context");

    inward_sign = $ps_face_boundary_seg_inward_is_positive_ey ? 1 : -1;
    face_span = _ps_fr_face_span2($ps_face_pts2d);

    _ps_fr_interference_boundary_cutter(face_span, z0, z1, inward_sign, eps);
}

// Default interference cutter for one true filled-boundary segment.
// Call inside place_on_face_filled_boundary_segments(...).
module ps_face_boundary_seg_default_interference_cutter_ctx(z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_boundary_seg_default_interference_cutter_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_boundary_seg2d), "ps_face_boundary_seg_default_interference_cutter_ctx: requires place_on_face_filled_boundary_segments context");
    assert(!is_undef($ps_face_boundary_seg_inward2d), "ps_face_boundary_seg_default_interference_cutter_ctx: requires place_on_face_filled_boundary_segments context");

    seg2d = $ps_face_boundary_seg2d;
    p0 = seg2d[0];
    p1 = seg2d[1];
    mid = _ps_seg_boundary_midpoint(seg2d);
    seg_len = _ps_seg_boundary_seg_len(seg2d);
    ex2 = (seg_len <= eps) ? [1, 0] : (p1 - p0) / seg_len;
    inward2 = $ps_face_boundary_seg_inward2d;
    dihed = $ps_face_boundary_seg_dihedral;
    ex3 = [ex2[0], ex2[1], 0];
    ez3 = _ps_seg_boundary_edge_ez3(inward2, dihed, eps);
    ey3 = v_norm(v_cross(ez3, ex3));

    multmatrix(ps_frame_matrix([mid[0], mid[1], 0], ex3, ey3, ez3))
        ps_face_boundary_seg_default_interference_cutter_local_ctx(z0, z1, eps);
}

// First-order face body for one source edge:
// owning filled-cell prism minus that edge's own default interference cutter.
// This stays in face-local coordinates and is intended as the phase-1 body
// used to trim other crossing edge cutters.
module ps_face_source_edge_phase1_body_ctx(source_edge_idx, z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts3d_local), "ps_face_source_edge_phase1_body_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    cells = ps_face_filled_cells($ps_face_pts3d_local, eps);
    segs = ps_face_filled_boundary_segments($ps_face_pts3d_local, eps);
    matches = [for (s = segs) if (s[1] == source_edge_idx) s];

    assert(len(matches) > 0, str("ps_face_source_edge_phase1_body_ctx: no boundary record for source edge ", source_edge_idx));

    union()
        for (rec = matches) {
            cell_idx = rec[4];
            assert(cell_idx >= 0 && cell_idx < len(cells), str("ps_face_source_edge_phase1_body_ctx: cell_idx out of range ", cell_idx));

            difference() {
                _ps_fr_filled_cell_prism(cells[cell_idx], z0, z1);

                place_on_face_filled_boundary_source_edges(source_edge_indices = [source_edge_idx], eps = eps)
                    ps_face_boundary_seg_default_interference_cutter_local_ctx(z0, z1, eps);
            }
        }
}

// Bounded lobe body used for phase-2 target-edge trimming:
// a lobe polygon extruded through the face slab, with phase-1 interference
// applied only to the real boundary source edges that bound that lobe.
module ps_face_crossing_lobe_body_ctx(poly2d, source_edge_ids, z0, z1, eps=1e-8) {
    difference() {
        translate([0, 0, z0])
            linear_extrude(height = z1 - z0)
                ps_polygon(points = poly2d);

        place_on_face_filled_boundary_source_edges(source_edge_indices = source_edge_ids, eps = eps)
            ps_face_boundary_seg_body_interference_cutter_local_ctx(z0, z1, eps);
    }
}

module place_on_face_interference_atom_boundary_cutters_ctx(atom_idx, eps=1e-8) {
    assert(!is_undef($ps_face_pts3d_local), "place_on_face_interference_atom_boundary_cutters_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);
    assert(atom_idx >= 0 && atom_idx < len(atoms), str("place_on_face_interference_atom_boundary_cutters_ctx: atom_idx out of range ", atom_idx));

    atom = atoms[atom_idx];
    edge_kinds = atom[3];
    n = len(atom[0]);

    for (ei = [0 : 1 : n - 1]) {
        if (_ps_seg_safe_at(edge_kinds, ei, "inner") == "boundary") {
            $ps_face_interference_atom_boundary_edge_idx = ei;
            children();
        }
    }
}

// Debug/helper surface for the actual interference cutters applied to a face.
// This exposes the same cutter union used by ps_face_interference_volume_ctx(...),
// but without differencing it from the face body.
module ps_face_interference_cutters_ctx(z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_interference_cutters_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_pts3d_local), "ps_face_interference_cutters_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);
    face_span = _ps_fr_face_span2($ps_face_pts2d);

    union() {
        for (atom = atoms)
            _ps_fr_atom_boundary_interference_cutters(atom, $ps_face_dihedrals, face_span, z0, z1, eps);
    }
}

module ps_face_interference_atom_prisms_ctx(z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts3d_local), "ps_face_interference_atom_prisms_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);

    union() {
        for (atom = atoms)
            _ps_fr_atom_prism(atom, z0, z1);
    }
}

module ps_face_interference_atom_prism_ctx(atom_idx, z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts3d_local), "ps_face_interference_atom_prism_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);
    assert(atom_idx >= 0 && atom_idx < len(atoms), str("ps_face_interference_atom_prism_ctx: atom_idx out of range ", atom_idx));

    _ps_fr_atom_prism(atoms[atom_idx], z0, z1);
}

module ps_face_interference_atom_cutters_ctx(atom_idx, z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_interference_atom_cutters_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_pts3d_local), "ps_face_interference_atom_cutters_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);
    face_span = _ps_fr_face_span2($ps_face_pts2d);
    assert(atom_idx >= 0 && atom_idx < len(atoms), str("ps_face_interference_atom_cutters_ctx: atom_idx out of range ", atom_idx));

    _ps_fr_atom_boundary_interference_cutters(atoms[atom_idx], $ps_face_dihedrals, face_span, z0, z1, eps);
}

module ps_face_interference_atom_volume_ctx(atom_idx, z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_interference_atom_volume_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_pts3d_local), "ps_face_interference_atom_volume_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);
    face_span = _ps_fr_face_span2($ps_face_pts2d);
    assert(atom_idx >= 0 && atom_idx < len(atoms), str("ps_face_interference_atom_volume_ctx: atom_idx out of range ", atom_idx));

    difference() {
        _ps_fr_atom_prism(atoms[atom_idx], z0, z1);
        _ps_fr_atom_boundary_interference_cutters(atoms[atom_idx], $ps_face_dihedrals, face_span, z0, z1, eps);
    }
}

// Interference-only face volume:
// the raw filled face slab, minus one-sided dihedral/2 cutters placed on the
// true filled-boundary subsegments. This intentionally excludes foreign cutters
// and local seat/clearance geometry.
module ps_face_interference_volume_ctx(z0, z1, eps=1e-8) {
    assert(!is_undef($ps_face_pts2d), "ps_face_interference_volume_ctx: requires place_on_faces context ($ps_face_pts2d)");
    assert(!is_undef($ps_face_pts3d_local), "ps_face_interference_volume_ctx: requires place_on_faces context ($ps_face_pts3d_local)");

    atoms = ps_face_filled_atoms($ps_face_pts3d_local, eps);
    face_span = _ps_fr_face_span2($ps_face_pts2d);

    union() {
        for (atom = atoms)
            difference() {
                _ps_fr_atom_prism(atom, z0, z1);
                _ps_fr_atom_boundary_interference_cutters(atom, $ps_face_dihedrals, face_span, z0, z1, eps);
            }
    }
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

function _ps_fr_run_end_plane(entry) =
    let(
        line = entry[0],
        n2 = line[0],
        d = line[1],
        keep_ge = line[2],
        n3 = [n2[0], n2[1], 0]
    )
    keep_ge ? [n3, d] : [-n3, -d];

function ps_face_visible_cell_region_planes(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8, include_run_ends=false) =
    let(
        z_min = min(z0, z1),
        z_max = max(z0, z1),
        cut_lo = is_undef(band_z0) ? z_min : band_z0,
        cut_hi = is_undef(band_z1) ? z_max : band_z1,
        pts = cell[0],
        m = len(pts),
        probe = _ps_fr_probe_point(pts, eps),
        z_mid = (z_min + z_max) / 2,
        run_end_entries = include_run_ends ? ps_face_visible_cell_cut_run_end_entries(cell, eps) : []
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
        ],
        [for (entry = run_end_entries) _ps_fr_run_end_plane(entry)]
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

function ps_face_visible_cell_loop_at_z_from_region_planes(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8, include_run_ends=false) =
    let(
        z_min = min(z0, z1),
        z_max = max(z0, z1),
        probe = _ps_fr_probe_point(cell[0], eps),
        planes = ps_face_visible_cell_region_planes(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0, band_z1, cut_clearance, eps, include_run_ends),
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
        cut_run_ids = (len(cell) > 5) ? cell[5] : [],
        n = len(pts),
        fwd = [
            for (k = [0:1:n-1])
                if (k == ia && ((ia + 1) % n) == ib)
                    [_ps_seg_safe_at(edge_ids, k, undef), _ps_seg_safe_at(kinds, k, "inner"), _ps_seg_safe_at(cut_ids, k, undef), _ps_seg_safe_at(cut_run_ids, k, undef)]
        ],
        rev = [
            for (k = [0:1:n-1])
                if (k == ib && ((ib + 1) % n) == ia)
                    [_ps_seg_safe_at(edge_ids, k, undef), _ps_seg_safe_at(kinds, k, "inner"), _ps_seg_safe_at(cut_ids, k, undef), _ps_seg_safe_at(cut_run_ids, k, undef)]
        ]
    )
    (len(fwd) > 0) ? fwd[0] :
    (len(rev) > 0) ? rev[0] :
    [undef, "inner", undef, undef];

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
                    [m01[2], m12[2], m20[2]],
                    [m01[3], m12[3], m20[3]]
                ]
        ];

function _ps_fr_atom_has_inner_edges(cell) =
    len([for (kind = cell[3]) if (kind == "inner") 1]) > 0;

function _ps_fr_atom_can_use_clipped_loops(cell, clip_loops, eps=1e-9) =
    let(
        clip_sizes = [for (lz = clip_loops) len(lz[0])]
    )
    !_ps_fr_atom_has_inner_edges(cell) &&
    min(clip_sizes) >= 3 &&
    len([for (n = clip_sizes) if (n != clip_sizes[0]) 1]) == 0;

function _ps_fr_cut_run_starts_at_edge(cell, edge_idx) =
    let(
        kinds = cell[3],
        run_ids = (len(cell) > 5) ? cell[5] : [],
        n = len(kinds),
        prev_e = (edge_idx - 1 + n) % n,
        run = _ps_seg_safe_at(run_ids, edge_idx, undef)
    )
    (_ps_seg_safe_at(kinds, edge_idx, undef) == "cut") &&
    (_ps_seg_safe_at(run_ids, prev_e, undef) != run);

function _ps_fr_cut_run_ends_at_edge(cell, edge_idx) =
    let(
        kinds = cell[3],
        run_ids = (len(cell) > 5) ? cell[5] : [],
        n = len(kinds),
        next_e = (edge_idx + 1) % n,
        run = _ps_seg_safe_at(run_ids, edge_idx, undef)
    )
    (_ps_seg_safe_at(kinds, edge_idx, undef) == "cut") &&
    (_ps_seg_safe_at(run_ids, next_e, undef) != run);

// Full-depth boundary planes for finite cut spans. These are local 2D line
// constraints; in 3D they lift to planes constant across z0..z1.
function ps_face_visible_cell_cut_run_end_entries(cell, eps=1e-9) =
    let(
        pts = cell[0],
        run_ids = (len(cell) > 5) ? cell[5] : [],
        n = len(pts)
    )
    [
        for (ei = [0:1:n-1])
            let(
                a = pts[ei],
                b = pts[(ei + 1) % n],
                dir = _ps_fr_norm2(b - a),
                run_id = _ps_seg_safe_at(run_ids, ei, undef)
            )
            each concat(
                _ps_fr_cut_run_starts_at_edge(cell, ei)
                    ? [[[dir, v_dot(dir, a), true], ei, true, run_id]]
                    : [],
                _ps_fr_cut_run_ends_at_edge(cell, ei)
                    ? [[[dir, v_dot(dir, b), false], ei, false, run_id]]
                    : []
            )
    ];

function ps_face_visible_cell_cut_run_end_lines(cell, eps=1e-9) =
    [for (e = ps_face_visible_cell_cut_run_end_entries(cell, eps)) e[0]];

function _ps_fr_cell_edge_inward_n(cell_pts2d, edge_idx, probe, eps=1e-9) =
    let(
        m = len(cell_pts2d),
        a = cell_pts2d[edge_idx],
        b = cell_pts2d[(edge_idx + 1) % m],
        e = _ps_fr_norm2(b - a),
        perp0 = _ps_fr_norm2([-e[1], e[0]]),
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
                probe = _ps_fr_probe_point(cell[0], eps),
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
        probe = _ps_fr_probe_point(pts, eps)
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
        probe = _ps_fr_probe_point(cell[0], eps),
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
    clip_loops = [
        for (z = levels)
            [ps_face_visible_cell_loop_at_z_clipped(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z, cut_lo, cut_hi, cut_clearance, eps), z]
    ];
    use_clipped = _ps_fr_atom_can_use_clipped_loops(cell, clip_loops, eps);
    active_loops = use_clipped ? clip_loops : loops;
    if (min([for (lz = active_loops) len(lz[0])]) >= 3)
        _ps_fr_stack_same_arity_loops(active_loops, eps);
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
    cell = [
        $ps_vis_seg_pts2d,
        undef,
        $ps_vis_seg_edge_ids,
        $ps_vis_seg_edge_kinds,
        is_undef($ps_vis_seg_cut_entry_ids) ? [] : $ps_vis_seg_cut_entry_ids,
        is_undef($ps_vis_seg_cut_run_ids) ? [] : $ps_vis_seg_cut_run_ids
    ];
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
