use <../../core/funcs.scad>
use <../../core/render.scad>
use <../../core/segments.scad>
use <../../core/validate.scad>

// Base in-plane inset applied to each face edge (mm).
FACE_PLATE_EDGE_INSET = 0.2;
// Minimum face radius before adding the pillow (mm).
FACE_PLATE_PILLOW_MIN_RAD = 5;
// Pillow inset at the face surface (mm).
FACE_PLATE_PILLOW_INSET = 2;
// Additional pillow inset at the raised height (mm).
FACE_PLATE_PILLOW_RAMP = 1;
// Pillow thickness above the face (mm).
FACE_PLATE_PILLOW_THK = 0.4;

// Thickness of the top polygon above the ramped edges.
FACE_PLATE_TOP_THK = 0.3;
// Clearance height for face sockets (mm) - make larger if face is far inset into a face.
FACE_PLATE_CLEAR_HEIGHT = 10;

function _line_intersect_2d(n1, d1, n2, d2) =
    let(det = n1[0] * n2[1] - n1[1] * n2[0])
    (abs(det) < 1e-9)
        ? [0, 0]
        : [
            (d1 * n2[1] - n1[1] * d2) / det,
            (n1[0] * d2 - d1 * n2[0]) / det
        ];


function ps_rot_axis(v, axis, ang) =
    let(
        a = v_norm(axis),
        c = cos(ang),
        s = sin(ang),
        term1 = v_scale(v, c),
        term2 = v_scale(v_cross(a, v), s),
        term3 = v_scale(a, v_dot(a, v) * (1 - c))
    )
    v_add(v_add(term1, term2), term3);


function _offset_pts2d(pts, diheds, ht) =
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
                    perp0 = v_norm([ -e[1], e[0] ]),
                    mid = (p0 + p1) / 2,
                    to_center = centroid - mid,
                    perp = (v_dot(perp0, to_center) < 0) ? -perp0 : perp0,
                    inset = -(ht/2) * tan((180 - diheds[k]) / 2),
                    d = v_dot(perp, p0) - inset
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
            _line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
    ];


function _offset_pts2d_inset(pts, inset) =
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
                    perp0 = v_norm([ -e[1], e[0] ]),
                    mid = (p0 + p1) / 2,
                    to_center = centroid - mid,
                    perp = (v_dot(perp0, to_center) < 0) ? -perp0 : perp0,
                    d = v_dot(perp, p0) + inset
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
            _line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
    ];

function _offset_pts2d_inset_edges(pts, insets) =
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
                    perp0 = v_norm([ -e[1], e[0] ]),
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
            _line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
    ];
    
function _gap_inset(dihed, gap) =
    let(s = sin(dihed / 2))
    (abs(s) < 1e-6) ? gap : (gap / (2 * s));


function _bottom_pts2d_from_bevel(top_pts, diheds, ht) =
    let(
        n = len(top_pts),
        centroid = v_scale([
            sum([for (p = top_pts) p[0]]),
            sum([for (p = top_pts) p[1]])
        ], 1 / n),
        lines = [
            for (k = [0:1:n-1])
                let(
                    p0 = top_pts[k],
                    p1 = top_pts[(k+1)%n],
                    e = v_norm([p1[0]-p0[0], p1[1]-p0[1], 0]),
                    n0 = [0,0,1],
                    // LHR: rotate outward normal toward adjacent face by +dihedral.
                    n1 = ps_rot_axis(n0, e, diheds[k]),
                    n_side = v_norm(n0 + n1),
                    n_xy = [n_side[0], n_side[1]],
                    d2_raw = v_dot(n_xy, p0) + n_side[2] * ht,
                    flip = (v_dot(n_xy, centroid) > d2_raw),
                    n_xy2 = flip ? [-n_xy[0], -n_xy[1]] : n_xy,
                    d2 = flip ? -d2_raw : d2_raw
                )
                [n_xy2, d2]
        ]
    )
    [
        for (k = [0:1:n-1])
            let(
                l0 = lines[(k-1+n)%n],
                l1 = lines[k]
            )
            _line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
    ];




function _bottom_pts2d_from_planes(pts, diheds, ht) =
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
                    perp0 = v_norm([ e[1], -e[0] ]),
                    mid = (p0 + p1) / 2,
                    to_center = centroid - mid,
                    perp = (v_dot(perp0, to_center) < 0) ? -perp0 : perp0,
                    alpha = (180 - diheds[k]) / 2,
                    n_xy = v_scale(perp, sin(alpha)),
                    n_z = cos(alpha),
                    d2 = v_dot(n_xy, p0) + n_z * ht
                )
                [n_xy, d2]
        ]
    )
    [
        for (k = [0:1:n-1])
            let(
                l0 = lines[(k-1+n)%n],
                l1 = lines[k]
            )
            _line_intersect_2d(l0[0], l0[1], l1[0], l1[1])
    ];


function _pts_centroid2d(pts) =
    let(n = len(pts))
    [sum([for (p = pts) p[0]]) / n, sum([for (p = pts) p[1]]) / n];

function _pts_radius2d(pts, centroid) =
    let(n = len(pts))
    sum([for (p = pts) norm([p[0] - centroid[0], p[1] - centroid[1]])]) / n;

function _safe_at(v, i, dflt) =
    (is_undef(i) || i < 0 || i >= len(v)) ? dflt : v[i];

function _segmented_inset_loops2d(pts, insets, eps=1e-8) =
    let(
        segs = ps_face_segments([for (p = pts) [p[0], p[1], 0]], "nonzero", eps)
    )
    [
        for (s = segs)
            let(
                s_pts = s[0],
                s_edge_ids = s[2],
                m = len(s_pts),
                s_insets = [
                    for (k = [0:1:m-1])
                        let(seg = [s_pts[k], s_pts[(k+1)%m]])
                        _ps_seg_is_parent_edge(seg, pts, eps)
                            ? _safe_at(insets, _safe_at(s_edge_ids, k, 0), 0)
                            : 0
                ]
            )
            _offset_pts2d_inset_edges(s_pts, s_insets)
    ];

function _segmented_body_loops(pts, diheds, insets, eps=1e-8) =
    let(
        segs = ps_face_segments([for (p = pts) [p[0], p[1], 0]], "nonzero", eps)
    )
    [
        for (s = segs)
            let(
                s_pts = s[0],
                s_edge_ids = s[2],
                m = len(s_pts),
                s_insets = [
                    for (k = [0:1:m-1])
                        let(seg = [s_pts[k], s_pts[(k+1)%m]])
                        _ps_seg_is_parent_edge(seg, pts, eps)
                            ? _safe_at(insets, _safe_at(s_edge_ids, k, 0), 0)
                            : 0
                ],
                s_diheds = [
                    for (k = [0:1:m-1])
                        let(seg = [s_pts[k], s_pts[(k+1)%m]])
                        _ps_seg_is_parent_edge(seg, pts, eps)
                            ? _safe_at(diheds, _safe_at(s_edge_ids, k, 0), 180)
                            : 180
                ],
                s_top2d = _offset_pts2d_inset_edges(s_pts, s_insets)
            )
            [s_top2d, s_diheds]
    ];


// Bevel by constructing an inset bottom polygon and lofting.
// Expects LHR/CW polygon order for pts and aligned dihedrals.
module face_plate(idx, pts, face_thk, diheds, insets_override, clear_space,
    edge_inset = FACE_PLATE_EDGE_INSET,
    pillow_min_rad = FACE_PLATE_PILLOW_MIN_RAD,
    pillow_inset = FACE_PLATE_PILLOW_INSET,
    pillow_ramp = FACE_PLATE_PILLOW_RAMP,
    pillow_thk = FACE_PLATE_PILLOW_THK,
    top_thk = FACE_PLATE_TOP_THK,
    base_z = undef,
    clear_height = FACE_PLATE_CLEAR_HEIGHT,
    eps = 1e-4
) {
    n = len(pts);
    insets = (is_undef(insets_override) || len(insets_override) != n)
        ? [for (k = [0:1:len(pts)-1]) _gap_inset(diheds[k], edge_inset)]
        : insets_override;
    body_loops = (edge_inset > 0)
        ? _segmented_body_loops(pts, diheds, insets, eps)
        : [[pts, diheds]];
    roof_loops = [for (bd = body_loops) bd[0]];
    pts_gap = roof_loops[0];
    centroid = _pts_centroid2d(pts_gap);
    rad = _pts_radius2d(pts_gap, centroid);

    base_z_eff = is_undef(base_z)? -face_thk / 2 : base_z; // base Z coord
    ramped_thk = face_thk - top_thk; // actual thickness of ramped part
    color(len(pts) == 3 ? "white" : "red") {
        // ramped part
        for (bd = body_loops) {
            loop2d = bd[0];
            loop_diheds = bd[1];
            m = len(loop2d);
            top = [for (p = loop2d) [p[0], p[1], base_z_eff + ramped_thk]];
            bottom2d = _bottom_pts2d_from_bevel(loop2d, loop_diheds, ramped_thk);
            bottom = [for (p = bottom2d) [p[0], p[1], base_z_eff]];
            points = concat(top, bottom);
            faces = concat(
                [[for (i = [0:1:m-1]) i]],
                [[for (i = [0:1:m-1]) (2*m-1-i)]],
                [for (i = [0:1:m-1]) [(i+1)%m, i, m + i]],
                [for (i = [0:1:m-1]) [(i+1)%m, m + i, m + (i+1)%m]]
            );
            polyhedron(points = points, faces = faces, convexity = 2);
        }

        // flat roof part: segment self-intersecting loops before insetting.
        translate([0, 0, base_z_eff + ramped_thk]) linear_extrude(top_thk)
            union() {
                for (loop = roof_loops)
                    ps_polygon(points = loop, mode = "nonzero");
            }

        // top pillow part
        if (rad > pillow_min_rad) {
            // Build the pillow from the same simple loops as the face body/roof.
            for (loop = roof_loops) {
                loop_centroid = _pts_centroid2d(loop);
                loop_rad = _pts_radius2d(loop, loop_centroid);
                if (loop_rad > pillow_inset + pillow_ramp + eps) {
                    s0 = max(0, 1 - pillow_inset / loop_rad);
                    s1 = max(0, 1 - (pillow_inset + pillow_ramp) / loop_rad);
                    p0 = [for (p = loop) [
                        loop_centroid[0] + (p[0] - loop_centroid[0]) * s0,
                        loop_centroid[1] + (p[1] - loop_centroid[1]) * s0
                    ]];
                    scale_xy = s0 <= eps ? 0 : (s1 / s0);
                    translate([0, 0, base_z_eff + face_thk])
                        linear_extrude(height = pillow_thk, scale = scale_xy)
                            ps_polygon(points = p0, mode = "nonzero");
                }
            }
        }
    }

    // Conditionally clears the airspace above the face, to remove material from the face-mount above the face
    if (clear_space) {
        color("magenta") translate([0, 0, base_z_eff + face_thk - eps]) linear_extrude(height = clear_height)
            union() {
                for (loop = roof_loops)
                    ps_polygon(points = loop, mode = "nonzero");
            }
    }
}



// Please keep this:
//difference() {
//    translate([0,0,-5]) cube(20);
//
//    !face_plate(0, [[10,0],[0,-10],[-10,0],[0,10]], face_thk=1.2, diheds=[140,80,140,80], insets_override=undef, clear_space=false);
    !face_plate(0, [for (i=[0:4]) [10*cos(-144*i), 10*sin(-144*i)]], face_thk=1.2, diheds=[60,80,100,120,140], insets_override=undef, clear_space=false);
//}
