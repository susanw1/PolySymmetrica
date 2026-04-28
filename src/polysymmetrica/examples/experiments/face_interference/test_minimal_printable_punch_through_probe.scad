use <../../../core/face_regions.scad>
use <../../../core/funcs.scad>
use <../../../core/placement.scad>
use <../../../core/prisms.scad>
use <../../../core/segments.scad>
use <../../../core/render.scad>

// Minimal pre-punch-through printable integration probe for
// poly_antiprism(7,3,15). This deliberately stays narrow: it combines the
// stable face-local data APIs into one positive keep-body without reintroducing
// the old proxy/carve stack.
//
// The final panel in each row computes:
//   raw face slab ∩ visible-cell mask ∩ positive anti-interference volume
//
// Orange cut strips are drawn as inspection aids only. They are not used as the
// primary integration strategy here, because this probe is testing the positive
// intersection model before any broader printable/proxy machinery returns.
// It does not yet model explicit foreign intrusion/clearance volumes.

IR = 32;
MODE = "nonzero";
FILTER_PARENT_CUTS = true;

STAR_FACE_IDX = 1;
TRI_FACE_IDX = 12;
STAR_SOURCE_EDGE_IDX = 0;

FACE_THK = 0.62;
VOL_Z_MIN = -1;
VOL_Z_MAX = 2;
MAX_PROJECT = 45;
LINE_R = 0.45;
TXT_H = 0.30;
TXT_S = 2.4;
PANEL_X = 116;
PANEL_Y = 112;
PANEL_LABEL_Y = -50;
CUT_KERF = 1.0;

P = poly_antiprism(7, 3, angle = 15);

/**
 * Function: Pick a stable, readable named color.
 * Params: i (index)
 * Returns: named color string
 */
function example_color(i) =
    (i % 8 == 0) ? "tomato" :
    (i % 8 == 1) ? "deepskyblue" :
    (i % 8 == 2) ? "gold" :
    (i % 8 == 3) ? "mediumseagreen" :
    (i % 8 == 4) ? "orchid" :
    (i % 8 == 5) ? "darkorange" :
    (i % 8 == 6) ? "turquoise" :
    "sienna";

/**
 * Function: Compute a unit 2D direction from one point to another.
 * Params: a/b (2D points)
 * Returns: unit direction, or `[1,0]` for a degenerate segment
 */
function unit2d(a, b) =
    let(d = [b[0] - a[0], b[1] - a[1]])
    (norm(d) <= 1e-9) ? [1, 0] : d / norm(d);

/**
 * Module: Draw text in the current local XY plane.
 * Params: s (text), size (font size), z (local Z)
 * Returns: none
 */
module draw_text2d(s, size = TXT_S, z = 0) {
    translate([0, 0, z])
        linear_extrude(height = TXT_H)
            text(s, size = size, halign = "center", valign = "center");
}

/**
 * Module: Draw a world-space panel label.
 * Params: s (label text)
 * Returns: none
 */
module draw_panel_label(s) {
    translate([0, PANEL_LABEL_Y, -18])
        draw_text2d(s, size = TXT_S + 0.35);
}

/**
 * Module: Draw a local 2D segment as an extruded rounded ribbon.
 * Params: seg2d (2D segment), r (ribbon radius), z (local Z offset)
 * Returns: none
 */
module draw_segment_stroke(seg2d, r = LINE_R, z = 0) {
    translate([0, 0, z])
        linear_extrude(height = FACE_THK * 5, center = true)
            hull() {
                translate(seg2d[0]) circle(r = r, $fn = 18);
                translate(seg2d[1]) circle(r = r, $fn = 18);
            }
}

/**
 * Module: Draw a local 2D polygon loop as a thin face fill.
 * Params: pts2d (2D polygon loop), z (local Z offset)
 * Returns: none
 */
module draw_polygon_fill(pts2d, z = 0) {
    translate([0, 0, z])
        linear_extrude(height = FACE_THK * 0.25, center = true)
            polygon(points = pts2d);
}

/**
 * Module: Draw labels for source edges on a selected face.
 * Params: pts2d (source face loop), highlight_source_edge_idx (optional source edge)
 * Returns: none
 */
module draw_source_edge_labels(pts2d, highlight_source_edge_idx = undef) {
    centre = ps_centroid2d(pts2d);
    for (ei = [0:1:len(pts2d)-1]) {
        seg2d = [pts2d[ei], pts2d[(ei + 1) % len(pts2d)]];
        mid = ps_segment_midpoint2d(seg2d);
        offset_dir = unit2d(mid, centre);
        highlighted = !is_undef(highlight_source_edge_idx) && ei == highlight_source_edge_idx;

        color(highlighted ? "red" : "dimgray")
            translate([mid[0] + 3.0 * offset_dir[0], mid[1] + 3.0 * offset_dir[1], FACE_THK * 0.65])
                draw_text2d(str(highlighted ? "*se" : "se", ei), size = 1.25);
    }
}

/**
 * Module: Emit a raw face-local material slab from the filled face polygon.
 * Params: none
 * Returns: none
 */
module face_material_slab() {
    translate([0, 0, -FACE_THK / 2])
        linear_extrude(height = FACE_THK)
            ps_polygon($ps_face_pts2d, MODE);
}

/**
 * Module: Emit the positive anti-interference volume for the current face.
 * Params: none
 * Returns: none
 */
module anti_interference_volume() {
    ps_face_anti_interference_volume(
        VOL_Z_MIN,
        VOL_Z_MAX,
        mode = MODE,
        max_project = MAX_PROJECT
    );
}

/**
 * Module: Emit the minimal printable keep-body under test.
 * Params: none
 * Returns: none
 */
module printable_keep_body() {
    intersection() {
        face_material_slab();
        face_visible_mask(FACE_THK, z_pad = 0.04, mode = MODE, filter_parent = FILTER_PARENT_CUTS);
        anti_interference_volume();
    }
}

/**
 * Module: Draw orange geometry-cut strips as inspection aids.
 * Params: none
 * Returns: none
 */
module draw_cut_strips() {
    cuts = ps_face_geom_cut_entries(
        $ps_face_pts2d,
        $ps_face_idx,
        $ps_poly_faces_idx,
        $ps_poly_verts_local,
        mode = MODE,
        filter_parent = FILTER_PARENT_CUTS
    );

    for (ci = [0:1:len(cuts)-1]) {
        cut = cuts[ci];
        mid = ps_segment_midpoint2d(cut[0]);

        color("darkorange", 0.70)
            draw_segment_stroke(cut[0], r = CUT_KERF, z = FACE_THK * 0.55);

        color("black")
            translate([mid[0], mid[1], FACE_THK * 0.95])
                draw_text2d(str("c", ci, "/f", cut[1]), size = 1.05);
    }
}

/**
 * Module: Draw the surrounding poly and label the selected face.
 * Params: face_idx (selected face), label_s (panel label)
 * Returns: none
 */
module draw_context_panel(face_idx, label_s) {
    poly_render(P, IR);

    color("silver")
        place_on_edges(P, IR)
            cube([$ps_edge_len, 0.35, 0.85], center = true);

    color("gold")
        place_on_vertices(P, IR)
            sphere(0.95, $fn = 12);

    place_on_faces(P, IR) {
        if ($ps_face_idx == face_idx)
            color("tomato", 0.48)
                face_material_slab();

        color(($ps_face_idx == face_idx) ? "black" : "gray")
            translate([0, 0, 1.3])
                draw_text2d(str("f", $ps_face_idx), size = 1.45);
    }

    draw_panel_label(label_s);
}

/**
 * Module: Draw visible cells and geometry-cut provenance for one face.
 * Params: face_idx (selected face), source_edge_idx (optional highlighted edge), label_s (panel label)
 * Returns: none
 */
module draw_visible_data_panel(face_idx, source_edge_idx, label_s) {
    place_on_faces(P, IR) {
        if ($ps_face_idx == face_idx) {
            color("gainsboro", 0.18)
                face_material_slab();

            place_on_face_visible_segments(mode = MODE, filter_parent = FILTER_PARENT_CUTS) {
                color(example_color($ps_vis_seg_idx), 0.40)
                    draw_polygon_fill($ps_vis_seg_pts2d, z = -FACE_THK * 0.35);

                centre = ps_centroid2d($ps_vis_seg_pts2d);
                color("white")
                    translate([centre[0], centre[1], FACE_THK * 0.72])
                        draw_text2d(str("v", $ps_vis_seg_idx), size = 1.25);
            }

            draw_cut_strips();
            draw_source_edge_labels($ps_face_pts2d, source_edge_idx);
        }
    }

    draw_panel_label(label_s);
}

/**
 * Module: Draw anti-interference volume and its boundary span skeleton.
 * Params: face_idx (selected face), source_edge_idx (optional highlighted source edge), label_s (panel label)
 * Returns: none
 */
module draw_volume_data_panel(face_idx, source_edge_idx, label_s) {
    place_on_faces(P, IR) {
        if ($ps_face_idx == face_idx) {
            color("gainsboro", 0.14)
                face_material_slab();

            color("deepskyblue", 0.35)
                anti_interference_volume();

            place_on_face_boundary_spans(mode = MODE) {
                highlighted = !is_undef(source_edge_idx) && $ps_boundary_span_source_edge_idx == source_edge_idx;
                color(highlighted ? "red" : "black")
                    cube([$ps_boundary_span_len, highlighted ? 0.9 : 0.55, 0.55], center = true);
            }
        }
    }

    draw_panel_label(label_s);
}

/**
 * Module: Draw the minimal printable keep-body plus cut provenance overlay.
 * Params: face_idx (selected face), source_edge_idx (optional highlighted source edge), label_s (panel label)
 * Returns: none
 */
module draw_printable_result_panel(face_idx, source_edge_idx, label_s) {
    place_on_faces(P, IR) {
        if ($ps_face_idx == face_idx) {
            color("white")
                printable_keep_body();

            color("deepskyblue", 0.18)
                anti_interference_volume();

            draw_cut_strips();
            draw_source_edge_labels($ps_face_pts2d, source_edge_idx);
        }
    }

    draw_panel_label(label_s);
}

/**
 * Module: Echo summary counts for one row.
 * Params: label_s (row label), face_idx (selected face)
 * Returns: none
 */
module echo_row_summary(label_s, face_idx) {
    place_on_faces(P, IR) {
        if ($ps_face_idx == face_idx) {
            bm = ps_face_boundary_model($ps_face_pts3d_local, MODE);
            cuts = ps_face_geom_cut_entries(
                $ps_face_pts2d,
                $ps_face_idx,
                $ps_poly_faces_idx,
                $ps_poly_verts_local,
                mode = MODE,
                filter_parent = FILTER_PARENT_CUTS
            );
            visible = ps_face_visible_segments(
                $ps_face_pts2d,
                $ps_face_idx,
                $ps_poly_faces_idx,
                $ps_poly_verts_local,
                mode = MODE,
                filter_parent = FILTER_PARENT_CUTS
            );

            echo(str(
                "minimal printable pre-punch-through ", label_s, " f", face_idx,
                ": boundary_loops=", len(bm[2]),
                " boundary_spans=", len(bm[3]),
                " geom_cuts=", len(cuts),
                " visible_segments=", len(visible)
            ));
        }
    }
}

/**
 * Module: Draw one target face row.
 * Params: face_idx (selected face), source_edge_idx (optional highlighted source edge), label_s (row label), y (row offset)
 * Returns: none
 */
module draw_probe_row(face_idx, source_edge_idx, label_s, y) {
    translate([-1.5 * PANEL_X, y, 0])
        draw_context_panel(face_idx, str(label_s, " context"));

    translate([-0.5 * PANEL_X, y, 0])
        draw_visible_data_panel(face_idx, source_edge_idx, str(label_s, " visible/cuts"));

    translate([0.5 * PANEL_X, y, 0])
        draw_volume_data_panel(face_idx, source_edge_idx, str(label_s, " anti-volume"));

    translate([1.5 * PANEL_X, y, 0])
        draw_printable_result_panel(face_idx, source_edge_idx, str(label_s, " keep-body"));

    echo_row_summary(label_s, face_idx);
}

draw_probe_row(STAR_FACE_IDX, STAR_SOURCE_EDGE_IDX, "star", 0.5 * PANEL_Y);
draw_probe_row(TRI_FACE_IDX, undef, "tri", -0.5 * PANEL_Y);
