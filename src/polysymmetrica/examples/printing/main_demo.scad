use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>
use <../../core/classify.scad>
use <../../core/prisms.scad>
use <../../core/construction.scad>
use <../../core/segments.scad>
use <../../core/face_regions.scad>

use <../../models/platonics_all.scad>
use <../../models/archimedians_all.scad>
use <../../models/johnsons_all.scad>

use <edge_seg.scad>
use <face_plate.scad>

SC = 1;
IR = 20 * SC;

//p = (tetrahedron());
//p = poly_truncate(dodecahedron());
//p = (dodecahedron());
//base = icosahedron();
//sol = solve_cantitruncate_trig(base);
//s = poly_cantellate_norm(base, 0.5);
//p = poly_dual(great_rhombicuboctahedron());
//p = poly_attach(octahedron(), icosahedron(), f1=[0,7]);
//p = poly_attach(octahedron(), icosahedron(), f1=[0,1,2,3,4,5,6,7]);
//p = poly_attach(p1, icosahedron(), f1=0);

//AP_CANT_ROWS = [
//    ["face", "family", 0, ["df", 0.30]],
//    ["face", "family", 1, ["df", 0.30]]
//];

//p = poly_cantellate(poly_antiprism(6), params_overrides=AP_CANT_ROWS);
//p = poly_prism(5);
//p = poly_antiprism(5);
//p = poly_prism(n=5, p=2);
//p = poly_antiprism(n=5, p=2, angle = 0);
//p = poly_antiprism(n=7, p=2, angle = 0);
p = poly_antiprism(n=7, p=3, angle = 15);

//p = j1_square_pyramid();l
//p = poly_dual(j2_pentagonal_pyramid());

EDGE_T = 3.5 * SC; // 3.5
FACE_T = 1.6 * SC; // 1.6 * SC;
FIN_T = 1 * SC;
INSET = 1.1 * SC;
FACET_BASE_T = 1;
FACET_BASE_W = 2.2;
BASE_Z = -FACE_T / 4;

DEBUG_SHOW_MODEL = true;
DEBUG_SEG_LABELS = true;
DEBUG_SEG_EDGE_LABELS = true;
DEBUG_LOOP_COMPARE = false;
DEBUG_LOOP_COMPARE_FACES = [2,5];
DEBUG_LOOP_COMPARE_CELLS = ["2/0"];
DEBUG_LOOP_COMPARE_Z = BASE_Z + FACE_T + 0.6;
DEBUG_LOOP_COMPARE_DRAW_Z = BASE_Z + FACE_T + 2.5;
DEBUG_LIVE_FACE_SLAB = false;
DEBUG_LIVE_FACE_SLAB_FACES = [2,5];
DEBUG_LIVE_FACE_SLAB_CELLS = ["2/0", "5/1"];
DEBUG_LIVE_FACE_SLAB_Z0 = BASE_Z;
DEBUG_LIVE_FACE_SLAB_Z1 = BASE_Z + FACE_T;
DEBUG_LABEL_Z = BASE_Z + FACE_T + 0.15;
DEBUG_LABEL_T = 0.15;
DEBUG_LABEL_SIZE = 1.2;
DEBUG_EDGE_LABEL_SIZE = 0.8;
DEBUG_PLANE_STRIP_EXTEND = 50;
DEBUG_LOOP_STRIP_W = 0.35;

//// Or, experimental:
//IR = 12 * SC;
//EDGE_T = 3.0 * SC; // 3.5
//FACE_T = 1.2 * SC; // 1.6 * SC;

function _sum_pairs2(pts) =
    [for (axis = [0:1]) sum([for (p = pts) p[axis]])];

function _centroid2(pts) =
    let(n = len(pts))
        assert(n > 0, "_centroid2: empty point set")
        [for (axis = [0:1]) _sum_pairs2(pts)[axis] / n];

function _mid2(a, b) = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];

function _has_str(xs, target) =
    len([for (x = xs) if (str(x) == target) 1]) > 0;

module _debug_edge_strip(a, b, w = DEBUG_LOOP_STRIP_W, extend = DEBUG_PLANE_STRIP_EXTEND) {
    d = b - a;
    L = norm(d);
    if (L > 1e-8) {
        ux = d[0] / L;
        uy = d[1] / L;
        nx = -uy;
        ny = ux;
        a2 = [a[0] - ux * extend, a[1] - uy * extend];
        b2 = [b[0] + ux * extend, b[1] + uy * extend];
        polygon(points = [
            [a2[0] + nx * w / 2, a2[1] + ny * w / 2],
            [b2[0] + nx * w / 2, b2[1] + ny * w / 2],
            [b2[0] - nx * w / 2, b2[1] - ny * w / 2],
            [a2[0] - nx * w / 2, a2[1] - ny * w / 2]
        ]);
    }
}

module debug_seg_labels(show_faces = undef) {
    place_on_faces(p, IR) {
        if (is_undef(show_faces) || len(search($ps_face_idx, [for (i=show_faces) i])) > 0) {
            cut_entries = ps_face_geom_cut_entries(
                $ps_face_pts2d,
                $ps_face_idx,
                $ps_poly_faces_idx,
                $ps_poly_verts_local,
                1e-8,
                "nonzero",
                true
            );

            place_on_face_visible_segments() {
                color("black")
                translate([
                    _centroid2($ps_vis_seg_pts2d)[0],
                    _centroid2($ps_vis_seg_pts2d)[1],
                    DEBUG_LABEL_Z
                ])
                    linear_extrude(DEBUG_LABEL_T)
                        text(
                            str("f", $ps_face_idx, "/c", $ps_vis_seg_idx),
                            size = DEBUG_LABEL_SIZE,
                            halign = "center",
                            valign = "center"
                        );

                if (DEBUG_SEG_EDGE_LABELS) {
                    for (ei = [0:1:len($ps_vis_seg_pts2d)-1]) {
                        if ($ps_vis_seg_edge_kinds[ei] == "cut") {
                            cid = $ps_vis_seg_cut_entry_ids[ei];
                            cutter_face = cut_entries[cid][1];

                            color("purple")
                            translate([
                                _mid2(
                                    $ps_vis_seg_pts2d[ei],
                                    $ps_vis_seg_pts2d[(ei + 1) % len($ps_vis_seg_pts2d)]
                                )[0],
                                _mid2(
                                    $ps_vis_seg_pts2d[ei],
                                    $ps_vis_seg_pts2d[(ei + 1) % len($ps_vis_seg_pts2d)]
                                )[1],
                                DEBUG_LABEL_Z + DEBUG_LABEL_T
                            ])
                                linear_extrude(DEBUG_LABEL_T)
                                    text(
                                        str("e", ei, "/c", cid, "/f", cutter_face),
                                        size = DEBUG_EDGE_LABEL_SIZE,
                                        halign = "center",
                                        valign = "center"
                                    );
                        }
                    }
                }
            }
        }
    }
}

module _debug_loop_overlay(pts, col, label_txt, z = DEBUG_LOOP_COMPARE_DRAW_Z, w = DEBUG_LOOP_STRIP_W) {
    if (len(pts) >= 2) {
        color(col)
        translate([0, 0, z])
            linear_extrude(DEBUG_LABEL_T)
                union() {
                    for (i = [0:1:len(pts)-1]) {
                        _debug_edge_strip(pts[i], pts[(i + 1) % len(pts)], w, 0);
                    }
                }

        color(col)
        translate([
            _centroid2(pts)[0],
            _centroid2(pts)[1],
            z + DEBUG_LABEL_T
        ])
            linear_extrude(DEBUG_LABEL_T)
                text(
                    label_txt,
                    size = DEBUG_EDGE_LABEL_SIZE,
                    halign = "center",
                    valign = "center"
                );
    }
}

module debug_loop_compare(show_faces = undef, show_cells = undef, z_sample = DEBUG_LOOP_COMPARE_Z) {
    place_on_faces(p, IR) {
        if (is_undef(show_faces) || len(search($ps_face_idx, [for (i=show_faces) i])) > 0) {
            cut_entries = ps_face_geom_cut_entries(
                $ps_face_pts2d,
                $ps_face_idx,
                $ps_poly_faces_idx,
                $ps_poly_verts_local,
                1e-8,
                "nonzero",
                true
            );

            place_on_face_visible_segments() {
                if (
                    is_undef(show_cells) ||
                    _has_str(show_cells, str($ps_face_idx, "/", $ps_vis_seg_idx))
                ) {
                    cell = [
                        $ps_vis_seg_pts2d,
                        undef,
                        $ps_vis_seg_edge_ids,
                        $ps_vis_seg_edge_kinds,
                        $ps_vis_seg_cut_entry_ids
                    ];
                    adj_loop = _ps_fr_visible_cell_loop_at_z(
                        cell,
                        $ps_face_dihedrals,
                        cut_entries,
                        $ps_poly_faces_idx,
                        $ps_poly_verts_local,
                        z_sample,
                        BASE_Z,
                        BASE_Z + FACE_T + 0.6,
                        INSET,
                        1e-8
                    );
                    clip_loop = ps_face_visible_cell_loop_at_z_clipped(
                        cell,
                        $ps_face_dihedrals,
                        cut_entries,
                        $ps_poly_faces_idx,
                        $ps_poly_verts_local,
                        z_sample,
                        BASE_Z,
                        BASE_Z + FACE_T + 0.6,
                        INSET,
                        1e-8
                    );

                    _debug_loop_overlay(adj_loop, "orange", str("adj z=", z_sample));
                    _debug_loop_overlay(clip_loop, "cyan", str("clip z=", z_sample), DEBUG_LOOP_COMPARE_DRAW_Z + 0.25);
                }
            }
        }
    }
}

module debug_live_face_slab(show_faces = undef, show_cells = undef) {
    place_on_faces(p, IR) {
        if (is_undef(show_faces) || len(search($ps_face_idx, [for (i=show_faces) i])) > 0) {
            band_z0 = DEBUG_LIVE_FACE_SLAB_Z0;
            band_z1 = DEBUG_LIVE_FACE_SLAB_Z1;
            band_pad = max(0.6, band_z1 - band_z0);
            clip_z0 = band_z0 - band_pad;
            clip_z1 = band_z1 + band_pad;

            place_on_face_visible_segments("nonzero", 1e-8, true) {
                if (
                    is_undef(show_cells) ||
                    _has_str(show_cells, str($ps_face_idx, "/", $ps_vis_seg_idx))
                ) {
                    color([0.2, 1.0, 0.3, 0.45])
                    ps_clip_to_visible_face_cell_ctx(
                        clip_z0,
                        clip_z1,
                        cut_clearance = INSET,
                        mode = "nonzero",
                        eps = 1e-8,
                        apply_cut_bands = true,
                        band_z0 = band_z0,
                        band_z1 = band_z1
                    )
                        translate([0, 0, DEBUG_LIVE_FACE_SLAB_Z0])
                            linear_extrude(DEBUG_LIVE_FACE_SLAB_Z1 - DEBUG_LIVE_FACE_SLAB_Z0)
                                ps_polygon(points = $ps_face_pts2d, mode = "nonzero");
                }
            }
        }
    }
}


/**
* Generate skeletal frame:
    show_faces = undef;

* Generate single face, eg face#8:
    show_faces = [8];
    and place '!'

* Generate multiple faces, eg face#0,#9, and #23:
    show_faces = [0, 9, 23];
    and place '!'
*/
module model(show_faces = undef, clear_airspace = true) {
    let (inter_radius = IR)
    union() {
        union() {
            // Constructs the edge-based frame
            color("gray")
            place_on_edges(p, IR) {
                edge_seg($ps_edge_pts_local, $ps_poly_center_local, edge_t = EDGE_T);
            }

            // mounting plate - the edge frame isn't quite substantial enough
            color("blue")
            place_on_faces(p, IR) {
                translate([0,0, BASE_Z - FACET_BASE_T]) linear_extrude(FACET_BASE_T)
                    difference() {
                        ps_polygon(points = $ps_face_pts2d);
                        offset(-FACET_BASE_W) ps_polygon(points = $ps_face_pts2d);
                    }
            }
        }
        // Constructs faces, removes them from frame to create face-fitting sockets.
        place_on_faces(p, IR) {
            color($ps_face_idx < 2? "red" : ["yellow", "green", "blue", "white", "orange"][$ps_face_idx % 5], 1)
            // add '!' here to force faces-only:
            if (is_undef(show_faces) || len(search($ps_face_idx, [for (i=show_faces) i])) > 0) {
                face_plate_visible($ps_face_idx, $ps_face_pts2d, FACE_T, $ps_face_dihedrals, undef, clear_airspace,
                    edge_inset = INSET, base_z = BASE_Z, clear_height = 0.6,
                    seg_apply_cut_bands = true);
            }
        }
    }
}

if (DEBUG_SHOW_MODEL) {
    model(undef, true);
}

if (DEBUG_SEG_LABELS) {
    debug_seg_labels(DEBUG_LOOP_COMPARE_FACES);
}

if (DEBUG_LOOP_COMPARE) {
    debug_loop_compare(DEBUG_LOOP_COMPARE_FACES, DEBUG_LOOP_COMPARE_CELLS);
}

if (DEBUG_LIVE_FACE_SLAB) {
    debug_live_face_slab(DEBUG_LIVE_FACE_SLAB_FACES, DEBUG_LIVE_FACE_SLAB_CELLS);
}
