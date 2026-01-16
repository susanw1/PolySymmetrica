use <../../core/funcs.scad>

// Base in-plane inset applied to each face edge (mm).
FACE_PLATE_EDGE_INSET = 0.2;
// Minimum face radius before adding the pillow (mm).
FACE_PLATE_PILLOW_MIN_RAD = 5;
// Pillow inset at the face surface (mm).
FACE_PLATE_PILLOW_INSET = 2;
// Additional pillow inset at the raised height (mm).
FACE_PLATE_PILLOW_RAMP = 1;
// Pillow height above the face (mm).
FACE_PLATE_PILLOW_H = 0.4;
// Clearance height for facet sockets (mm) - make larger if face is far inset into a face.
FACE_PLATE_CLEAR_HEIGHT = 10;

function _line_intersect_2d(n1, d1, n2, d2) =
    let(det = n1[0] * n2[1] - n1[1] * n2[0])
    (abs(det) < 1e-9)
        ? [0, 0]
        : [
            (d1 * n2[1] - n1[1] * d2) / det,
            (n1[0] * d2 - d1 * n2[0]) / det
        ];


function _rot_axis(v, axis, ang) =
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
                    n1 = _rot_axis(n0, e, -diheds[k]),
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


// Bevel by constructing an inset bottom polygon and lofting.
module face_plate(idx, pts, ht, diheds, insets_override, clear_space,
    edge_inset = FACE_PLATE_EDGE_INSET,
    pillow_min_rad = FACE_PLATE_PILLOW_MIN_RAD,
    pillow_inset = FACE_PLATE_PILLOW_INSET,
    pillow_ramp = FACE_PLATE_PILLOW_RAMP,
    pillow_h = FACE_PLATE_PILLOW_H,
    clear_height = FACE_PLATE_CLEAR_HEIGHT
) {
    n = len(pts);
    insets = (is_undef(insets_override) || len(insets_override) != n)
        ? [for (k = [0:1:len(pts)-1]) _gap_inset(diheds[k], edge_inset)]
        : insets_override;
    pts_gap = (edge_inset > 0) ? _offset_pts2d_inset_edges(pts, insets) : pts;
    centroid = _pts_centroid2d(pts_gap);
    rad = _pts_radius2d(pts_gap, centroid);
    top = [for (p = pts_gap) [p[0], p[1], ht/2]];
    bottom2d = _bottom_pts2d_from_bevel(pts_gap, diheds, ht);
    bottom = [for (p = bottom2d) [p[0], p[1], -ht/2]];
    points = concat(top, bottom);
    faces = concat(
        [ [for (i = [0:1:n-1]) i] ],
        [ [for (i = [0:1:n-1]) (2*n-1-i)] ],
        [for (i = [0:1:n-1]) [i, n + (i+1)%n, (i+1)%n]],
        [for (i = [0:1:n-1]) [i, n + i, n + (i+1)%n]]
    );

    echo("points", points, "faces", faces);
    color(len(pts) == 3 ? "white" : "red") {
        polyhedron(points = points, faces = faces);
        if (rad > pillow_min_rad) {
            // Cheap convex "pillow" by hulling two inset offsets at different heights.
            hull() {
                translate([0, 0, ht/2])
                    linear_extrude(height = 0.01)
                        offset(delta = -pillow_inset)
                            polygon(points = pts_gap);
                translate([0, 0, ht/2 + pillow_h])
                    linear_extrude(height = 0.01)
                        offset(delta = -(pillow_inset + pillow_ramp))
                            polygon(points = pts_gap);
            }
        }
    }

    if (clear_space) {
        translate([0, 0, ht/2]) linear_extrude(height = clear_height) polygon(points = pts_gap);
    }
}

