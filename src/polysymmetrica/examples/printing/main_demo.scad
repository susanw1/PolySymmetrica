use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>


t = undef;
IR = 20;

//p = (tetrahedron());
//p = poly_truncate(octahedron());
//p = (dodecahedron());
p = poly_truncate(dodecahedron());

EDGE_T = 3;
FACE_T = 0.8;
FIN_T = 0.8;

SINGLE_FACE = 0;

PILLOW_MIN_RAD = 5;
PILLOW_INSET = 2;
PILLOW_RAMP = 1;
PILLOW_H = 0.4;

module edge_seg(len, pc, r_poly) {
    translate([0, 0, 0]) {
        // builds the edge bar itself
        hull() {
            translate([len/2, 0, 0]) sphere(d = EDGE_T, $fn = 40);
            translate([-len/2, 0, 0]) sphere(d = EDGE_T, $fn = 40);
//            translate([0, 0, -EDGE_T/4]) cube([len, EDGE_T, EDGE_T/2], center = true);
        }
        
        // creates the support fin under the edge
        hull() {
            translate([len/2,0,0]) sphere(d = FIN_T, $fn = 10);
            translate([-len/2,0,0]) sphere(d = FIN_T, $fn = 10);
            translate(pc) sphere(d = FIN_T, $fn = 10);
        }
    }    
} 


function _line_intersect_2d(n1, d1, n2, d2) =
    let(det = n1[0] * n2[1] - n1[1] * n2[0])
    (abs(det) < 1e-9)
        ? [0, 0]
        : [
            (d1 * n2[1] - n1[1] * d2) / det,
            (n1[0] * d2 - d1 * n2[0]) / det
        ];

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

function _pts_centroid2d(pts) =
    let(n = len(pts))
    [sum([for (p = pts) p[0]]) / n, sum([for (p = pts) p[1]]) / n];

function _pts_radius2d(pts, centroid) =
    let(n = len(pts))
    sum([for (p = pts) norm([p[0] - centroid[0], p[1] - centroid[1]])]) / n;

// Bevel by constructing an inset bottom polygon and lofting.
module face_plate(idx, pts, ht, diheds, clear_space) {
    n = len(pts);
    centroid = _pts_centroid2d(pts);
    rad = _pts_radius2d(pts, centroid);
    top = [for (p = pts) [p[0], p[1], ht/2]];
    bottom2d = _offset_pts2d(pts, diheds, ht);
    bottom = [for (p = bottom2d) [p[0], p[1], -ht/2]];
    points = concat(top, bottom);
    faces = concat(
        [ [for (i = [0:1:n-1]) i] ],
        [ [for (i = [0:1:n-1]) (2*n-1-i)] ],
        [for (i = [0:1:n-1]) [i, n + (i+1)%n, (i+1)%n]],
        [for (i = [0:1:n-1]) [i, n + i, n + (i+1)%n]]
    );

    echo("points", points, "faces", faces);
    color(len(pts) == 3 ? "white" : "red")
        union() {
            polyhedron(points = points, faces = faces);
            if (rad > PILLOW_MIN_RAD) {
                // Cheap convex "pillow" by hulling two inset offsets at different heights.
                hull() {
                    translate([0, 0, ht/2])
                        linear_extrude(height = 0.01)
                            offset(delta = -PILLOW_INSET)
                                polygon(points = pts);
                    translate([0, 0, ht/2 + PILLOW_H])
                        linear_extrude(height = 0.01)
                            offset(delta = -(PILLOW_INSET + PILLOW_RAMP))
                                polygon(points = pts);
                }
            }
        }

    if (clear_space) {
        translate([0, 0, ht/2]) linear_extrude(height = EDGE_T) polygon(points = pts);
    }
}


difference() {
    if (is_undef(SINGLE_FACE)) {
        color("gray")
        place_on_edges(p, IR) {
            edge_seg($ps_edge_len, $ps_poly_center_local, $ps_edge_midradius);
        }
    }
    
    !place_on_faces(p, IR) {
        echo($ps_vertex_count, $ps_facet_radius, $ps_facet_dihedrals, $ps_face_pts2d);

        if (is_undef(SINGLE_FACE) || $ps_facet_idx == SINGLE_FACE || $ps_facet_idx == 1) {
//            face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_facet_dihedrals, is_undef(SINGLE_FACE));
            face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_facet_dihedrals, false);
        }
    }

    sphere(r = IR * 5/6);
}

*edge_seg(10, [0,0,-20], 20);
