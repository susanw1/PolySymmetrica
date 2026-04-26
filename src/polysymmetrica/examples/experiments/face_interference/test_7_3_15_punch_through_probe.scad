use <../../../core/face_regions.scad>
use <../../../core/funcs.scad>
use <../../../core/placement.scad>
use <../../../core/prisms.scad>
use <../../../core/segments.scad>

// Dedicated self-crossing punch-through probe for poly_antiprism(7, 3, angle=15).
// This is a developer experiment, not a printable example. It keeps the target
// face data visible as arrangement cells, filled boundaries, cut segments, and
// the positive anti-interference volume so punch-through behavior can be
// inspected without reintroducing the full proxy/carve stack.

IR = 32;
STAR_FACE_IDX = 1;
STAR_SOURCE_EDGE_IDX = 0;
TRI_FACE_IDX = 12;
TRI_SOURCE_EDGE_IDX = undef;
MODE = "nonzero";
FILTER_PARENT_CUTS = true;

FACE_THK = 0.28;
NODE_R = 0.85;
LINE_R = 0.42;
TXT_H = 0.32;
TXT_S = 2.4;
PANEL_X = 112;
PANEL_Y = 112;
PANEL_LABEL_Y = -50;
VOL_Z_MIN = -2;
VOL_Z_MAX = 2;
MAX_PROJECT = 45;

P = poly_antiprism(7, 3, angle = 15);

/**
 * Function: Pick a stable, readable named color.
 * Params: i (index)
 * Returns: named color string
 */
function draw_color(i) =
    (i % 8 == 0) ? "tomato" :
    (i % 8 == 1) ? "deepskyblue" :
    (i % 8 == 2) ? "gold" :
    (i % 8 == 3) ? "mediumseagreen" :
    (i % 8 == 4) ? "orchid" :
    (i % 8 == 5) ? "darkorange" :
    (i % 8 == 6) ? "turquoise" :
    "sienna";

/**
 * Function: Compute a 2D segment midpoint.
 * Params: seg2d (`[[x0,y0],[x1,y1]]`)
 * Returns: midpoint `[x, y]`
 */
function segment_midpoint2d(seg2d) =
    [(seg2d[0][0] + seg2d[1][0]) / 2, (seg2d[0][1] + seg2d[1][1]) / 2];

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
        draw_text2d(s, size = TXT_S + 0.4);
}

/**
 * Module: Draw a local 2D segment as an extruded rounded ribbon.
 * Params: seg2d (2D segment), r (ribbon radius), z (local Z offset)
 * Returns: none
 */
module draw_segment_stroke(seg2d, r = LINE_R, z = 0) {
    translate([0, 0, z])
        linear_extrude(height = FACE_THK, center = true)
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
        linear_extrude(height = FACE_THK, center = true)
            polygon(points = pts2d);
}

/**
 * Module: Draw labels for the original source edges of the selected face.
 * Params: pts2d (source face loop)
 * Returns: none
 */
module draw_source_edge_labels(pts2d, source_edge_idx = undef) {
    centre = ps_centroid2d(pts2d);
    for (ei = [0:1:len(pts2d)-1]) {
        seg2d = [pts2d[ei], pts2d[(ei + 1) % len(pts2d)]];
        mid = segment_midpoint2d(seg2d);
        offset_dir = unit2d(mid, centre);
        is_target_source = !is_undef(source_edge_idx) && ei == source_edge_idx;
        label = is_target_source ? str("*se", ei) : str("se", ei);

        color(is_target_source ? "red" : "dimgray")
            translate([mid[0] + 3.0 * offset_dir[0], mid[1] + 3.0 * offset_dir[1], 1.1 * FACE_THK])
                draw_text2d(label, size = 1.35);
    }
}

/**
 * Module: Draw global poly context with face and edge labels.
 * Params: poly (poly descriptor), face_idx (selected face)
 * Returns: none
 */
module draw_poly_context(poly, face_idx) {
    color("silver")
        place_on_edges(poly, IR)
            cube([$ps_edge_len, 0.35, 0.85], center = true);

    color("gold")
        place_on_vertices(poly, IR)
            sphere(1.0, $fn = 12);

    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx)
            color("tomato", 0.46)
                linear_extrude(height = FACE_THK, center = true)
                    ps_polygon($ps_face_pts2d, MODE);

        color(($ps_face_idx == face_idx) ? "black" : "gray")
            translate([0, 0, 1.3])
                draw_text2d(str("f", $ps_face_idx), size = 1.55);
    }

    place_on_edges(poly, IR) {
        if (_ps_list_contains($ps_edge_adj_faces_idx, face_idx))
            color("black")
                translate([0, 0, 1.0])
                    draw_text2d(str("e", $ps_edge_idx), size = 1.35);
    }
}

/**
 * Module: Echo summary counts for the target face once during render.
 * Params: label_s (case label), face_pts3d_local, face_idx, poly_faces_idx, poly_verts_local
 * Returns: none
 */
module draw_echo_target_summary(label_s, face_pts3d_local, face_idx, poly_faces_idx, poly_verts_local) {
    arr = ps_face_arrangement(face_pts3d_local);
    bm = ps_face_boundary_model(face_pts3d_local, MODE);
    cuts = ps_face_geom_cut_entries(
        ps_xy(face_pts3d_local),
        face_idx,
        poly_faces_idx,
        poly_verts_local,
        mode = MODE,
        filter_parent = FILTER_PARENT_CUTS
    );
    visible = ps_face_visible_segments(
        ps_xy(face_pts3d_local),
        face_idx,
        poly_faces_idx,
        poly_verts_local,
        mode = MODE,
        filter_parent = FILTER_PARENT_CUTS
    );

    echo(str(
        "7/3+15 punch-through probe ", label_s, " f", face_idx,
        ": crossings=", len(arr[1]),
        " nodes=", len(arr[2]),
        " arrangement_spans=", len(arr[3]),
        " cells=", len(arr[4]),
        " filled_cells=", len(bm[1]),
        " boundary_loops=", len(bm[2]),
        " boundary_spans=", len(bm[3]),
        " geom_cuts=", len(cuts),
        " visible_segments=", len(visible)
    ));
}

/**
 * Module: Draw arrangement nodes, spans, and source-edge labels for the target face.
 * Params: poly (poly descriptor), face_idx (target face), source_edge_idx (optional highlighted source edge), label_s (panel label)
 * Returns: none
 */
module draw_arrangement_panel(poly, face_idx, source_edge_idx, label_s) {
    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            arr = ps_face_arrangement($ps_face_pts3d_local);
            nodes = arr[2];
            spans = arr[3];

            color("gainsboro", 0.28)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            for (si = [0:1:len(spans)-1]) {
                span = spans[si];
                is_target_source = !is_undef(source_edge_idx) && span[3] == source_edge_idx;
                color(is_target_source ? "red" : "slategray")
                    draw_segment_stroke(span[0], r = is_target_source ? LINE_R * 0.8 : LINE_R * 0.45);
            }

            for (ni = [0:1:len(nodes)-1]) {
                node = nodes[ni];
                color(node[1] == "crossing" ? "deepskyblue" : "black")
                    translate([node[0][0], node[0][1], FACE_THK])
                        sphere(r = NODE_R, $fn = 12);
            }

            draw_source_edge_labels($ps_face_pts2d, source_edge_idx);
            draw_echo_target_summary(label_s, $ps_face_pts3d_local, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local);
        }
    }

    draw_panel_label(str(label_s, " arrangement"));
}

/**
 * Module: Draw nonzero filled boundary loops and source-edge lineage.
 * Params: poly (poly descriptor), face_idx (target face), source_edge_idx (optional highlighted source edge), label_s (panel label)
 * Returns: none
 */
module draw_boundary_panel(poly, face_idx, source_edge_idx, label_s) {
    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            bm = ps_face_boundary_model($ps_face_pts3d_local, MODE);
            loops = bm[2];
            spans = bm[3];

            color("gainsboro", 0.18)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            for (li = [0:1:len(loops)-1])
                color(draw_color(li), 0.38)
                    draw_polygon_fill(loops[li][0], z = -0.15);

            for (si = [0:1:len(spans)-1]) {
                span = spans[si];
                is_target_source = !is_undef(source_edge_idx) && span[2] == source_edge_idx;
                mid = segment_midpoint2d(span[0]);

                color(is_target_source ? "red" : "black")
                    draw_segment_stroke(span[0], r = is_target_source ? LINE_R * 0.95 : LINE_R * 0.58, z = 0.25);

                color("white")
                    translate([mid[0], mid[1], 1.3 * FACE_THK])
                        draw_text2d(str("b", si, "/se", span[2]), size = 1.15);
            }

            draw_source_edge_labels($ps_face_pts2d, source_edge_idx);
        }
    }

    draw_panel_label(str(label_s, " boundary"));
}

/**
 * Module: Draw geometry cuts from other faces and the retained visible cells.
 * Params: poly (poly descriptor), face_idx (target face), source_edge_idx (optional highlighted source edge), label_s (panel label)
 * Returns: none
 */
module draw_cut_visibility_panel(poly, face_idx, source_edge_idx, label_s) {
    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            color("gainsboro", 0.16)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            place_on_face_visible_segments(mode = MODE, filter_parent = FILTER_PARENT_CUTS) {
                color(draw_color($ps_vis_seg_idx), 0.35)
                    draw_polygon_fill($ps_vis_seg_pts2d, z = -0.05);

                centre = ps_centroid2d($ps_vis_seg_pts2d);
                color("white")
                    translate([centre[0], centre[1], 1.2 * FACE_THK])
                        draw_text2d(str("v", $ps_vis_seg_idx), size = 1.3);
            }

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
                mid = segment_midpoint2d(cut[0]);

                color("darkorange")
                    draw_segment_stroke(cut[0], r = LINE_R * 0.5, z = 0.45);

                color("black")
                    translate([mid[0], mid[1], 1.5 * FACE_THK])
                        draw_text2d(str("c", ci, "/f", cut[1]), size = 1.1);
            }

            draw_source_edge_labels($ps_face_pts2d, source_edge_idx);
        }
    }

    draw_panel_label(str(label_s, " cuts"));
}

/**
 * Module: Draw dihedral-aware boundary spans and positive anti-interference volume.
 * Params: poly (poly descriptor), face_idx (target face), source_edge_idx (optional highlighted source edge), label_s (panel label)
 * Returns: none
 */
module draw_volume_panel(poly, face_idx, source_edge_idx, label_s) {
    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            color("gainsboro", 0.18)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            color("deepskyblue", 0.32)
                ps_face_anti_interference_volume(
                    VOL_Z_MIN,
                    VOL_Z_MAX,
                    mode = MODE,
                    max_project = MAX_PROJECT
                );

            place_on_face_boundary_spans(mode = MODE) {
                is_target_source = !is_undef(source_edge_idx) && $ps_boundary_span_source_edge_idx == source_edge_idx;
                color(is_target_source ? "red" : "black")
                    cube([$ps_boundary_span_len, is_target_source ? 0.95 : 0.55, 0.55], center = true);

                color("white")
                    translate([0, 0, 1.1])
                        draw_text2d(
                            str(
                                "s", $ps_boundary_span_idx,
                                "/se", $ps_boundary_span_source_edge_idx,
                                "/f", $ps_boundary_span_adj_face_idx
                            ),
                            size = 0.95
                        );
            }
        }
    }

    draw_panel_label(str(label_s, " volume"));
}

module draw_probe_row(face_idx, source_edge_idx, label_s) {
    translate([-2 * PANEL_X, 0, 0])
        draw_poly_context(P, face_idx);

    translate([-1 * PANEL_X, 0, 0])
        draw_arrangement_panel(P, face_idx, source_edge_idx, label_s);

    translate([0 * PANEL_X, 0, 0])
        draw_boundary_panel(P, face_idx, source_edge_idx, label_s);

    translate([1 * PANEL_X, 0, 0])
        draw_cut_visibility_panel(P, face_idx, source_edge_idx, label_s);

    translate([2 * PANEL_X, 0, 0])
        draw_volume_panel(P, face_idx, source_edge_idx, label_s);
}

translate([0, 0.5 * PANEL_Y, 0])
    draw_probe_row(STAR_FACE_IDX, STAR_SOURCE_EDGE_IDX, "star");

translate([0, -0.5 * PANEL_Y, 0])
    draw_probe_row(TRI_FACE_IDX, TRI_SOURCE_EDGE_IDX, "tri");
