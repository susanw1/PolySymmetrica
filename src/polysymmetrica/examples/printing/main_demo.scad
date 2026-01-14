use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>


t = undef;
IR = 20;

LAYER1 = -100;

//p = (tetrahedron());
//p = (octahedron());
p = (dodecahedron());
//p = poly_truncate(dodecahedron());

EDGE_T = 3;
FACE_T = 0.8;
FIN_T = 0.8;

SINGLE_FACE = 0;

module edge_seg(len, pc, r_poly) {
    translate([0, 0, 0]) {
        // builds the edge bar itself
        hull() {
            translate([len/2, 0, 0]) sphere(d = EDGE_T, $fn = 40);
            translate([-len/2, 0, 0]) sphere(d = EDGE_T, $fn = 40);
            translate([0, 0, -EDGE_T/4]) cube([len, EDGE_T, EDGE_T/2], center = true);
        }
        
        // creates the support fin under the edge
        hull() {
            translate([len/2,0,0]) cube(FIN_T, center = true);
            translate([-len/2,0,0]) cube(FIN_T, center = true);
            translate(pc) cube(FIN_T, center = true);
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

// Bevel by constructing an inset bottom polygon and lofting.
module face_plate(idx, pts, ht, diheds, clear_space) {
    n = len(pts);
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
        polyhedron(points = points, faces = faces);

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
    
    place_on_faces(p, IR) {
        echo($ps_vertex_count, $ps_facet_radius, $ps_facet_dihedrals, $ps_face_pts2d);

        if (is_undef(SINGLE_FACE) || $ps_facet_idx == SINGLE_FACE) {
//            face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_facet_dihedrals, is_undef(SINGLE_FACE));
            face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_facet_dihedrals, false);
        }
    }

    sphere(r = IR * 3/4);
}


*edge_seg(10, [0,0,-20], 20);

pts = [[12.9968, 2.22045e-16, 0.4], [4.01623, 12.3607, 0.4], [-10.5146, 7.63932, 0.4], [-10.5146, -7.63932, 0.4], [4.01623, -12.3607, 0.4], [12.6912, 0, -0.4], [3.9218, 12.0701, -0.4], [-10.2674, 7.45971, -0.4], [-10.2674, -7.45971, -0.4], [3.9218, -12.0701, -0.4]];
faces = [[0, 1, 2, 3, 4], [9, 8, 7, 6, 5], [0, 6, 1], [1, 7, 2], [2, 8, 3], [3, 9, 4], [4, 5, 0], [0, 5, 6], [1, 6, 7], [2, 7, 8], [3, 8, 9], [4, 9, 5]];

*difference() {
    polyhedron(points = pts, faces = faces);
    translate([0,0,30]) #sphere(r = IR * 3/4);
}

*difference() {
    union() {
        place_on_faces(p, IR) {
            echo($ps_vertex_count, $ps_facet_radius, $ps_facet_dihedrals, $ps_face_pts2d);

            if (is_undef(SINGLE_FACE) || $ps_facet_idx == SINGLE_FACE) {
                face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_facet_dihedrals, false);
            }
        }
    }
    #sphere(r = IR * 3/4);
}
