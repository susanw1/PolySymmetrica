use <../../../core/face_regions.scad>
use <../../../core/funcs.scad>
use <../../../core/placement.scad>
use <../../../core/segments.scad>
use <../../../core/truncation.scad>
use <../../../core/render.scad>
use <../../../models/tetrahedron.scad>

// Canonical self-crossing stress case for the antitruncated tetrahedron.
// The selected hex face has short end spans and long crossing spans whose
// adjacent-face direction changes along the same source edge. This file keeps
// that data visible through arrangement, boundary, source-edge, span, and
// anti-interference-volume views.

IR = 34;
HEX_FACE_IDX = 0;
LONG_SOURCE_EDGE_IDX = 1;
SHORT_SOURCE_EDGE_IDX = 0;
MODE = "nonzero";

FACE_THK = 0.26;
LINE_R = 0.48;
NODE_R = 0.78;
TXT_H = 0.30;
TXT_S = 2.4;
PANEL_X = 118;
PANEL_LABEL_Y = -50;
VOL_Z_MIN = -2;
VOL_Z_MAX = 2;
MAX_PROJECT = 36;
DIH_RAY_LEN = 8;

P = poly_truncate(tetrahedron(), t = -0.5);

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
 * Module: Draw a 3D direction marker in the current local coordinate system.
 * Params: dir (3D direction), len (marker length), r (marker radius)
 * Returns: none
 */
module draw_direction3d(dir, len = DIH_RAY_LEN, r = 0.36) {
    if (!is_undef(dir) && norm(dir) > 1e-9) {
        d = v_norm(dir);
        axis = [-d[1], d[0], 0];
        angle = acos(d[2]);
        rot_axis = (norm(axis) <= 1e-9) ? [1, 0, 0] : axis;

        rotate(a = angle, v = rot_axis) {
            translate([0, 0, len / 2])
                cylinder(h = len, r = r, center = true, $fn = 14);
            translate([0, 0, len])
                sphere(r = r * 1.8, $fn = 12);
        }
    }
}

/**
 * Module: Draw labels for the selected face's original source edges.
 * Params: pts2d (source face loop)
 * Returns: none
 */
module draw_source_edge_labels(pts2d) {
    centre = ps_centroid2d(pts2d);
    for (ei = [0:1:len(pts2d)-1]) {
        seg2d = [pts2d[ei], pts2d[(ei + 1) % len(pts2d)]];
        mid = ps_segment_midpoint2d(seg2d);
        offset_dir = unit2d(mid, centre);
        is_long = ei == LONG_SOURCE_EDGE_IDX;
        is_short = ei == SHORT_SOURCE_EDGE_IDX;

        color(is_long ? "red" : (is_short ? "navy" : "dimgray"))
            translate([mid[0] + 3.0 * offset_dir[0], mid[1] + 3.0 * offset_dir[1], 1.1 * FACE_THK])
                draw_text2d(str(is_long ? "*se" : "se", ei), size = 1.35);
    }
}

/**
 * Module: Echo summary counts for the antitet hex face.
 * Params: face_pts3d_local, face_idx, poly_faces_idx, poly_verts_local
 * Returns: none
 */
module echo_hex_summary(face_pts3d_local, face_idx, poly_faces_idx, poly_verts_local) {
    arr = ps_face_arrangement(face_pts3d_local);
    bm = ps_face_boundary_model(face_pts3d_local, MODE);
    cuts = ps_face_geom_cut_entries(
        ps_xy(face_pts3d_local),
        face_idx,
        poly_faces_idx,
        poly_verts_local,
        mode = MODE
    );

    echo(str(
        "antitruncated tetrahedron hex f", face_idx,
        ": crossings=", len(arr[1]),
        " nodes=", len(arr[2]),
        " arrangement_spans=", len(arr[3]),
        " cells=", len(arr[4]),
        " filled_cells=", len(bm[1]),
        " boundary_loops=", len(bm[2]),
        " boundary_spans=", len(bm[3]),
        " geom_cuts=", len(cuts)
    ));
}

/**
 * Module: Draw global context for the antitruncated tetrahedron.
 * Params: none
 * Returns: none
 */
module draw_context_panel() {
    poly_render(P, IR);

    color("gold")
        place_on_vertices(P, IR)
            sphere(1.05, $fn = 12);

    color("silver")
        place_on_edges(P, IR)
            cube([$ps_edge_len, 0.4, 0.9], center = true);

    place_on_faces(P, IR) {
        color(($ps_face_idx == HEX_FACE_IDX) ? "black" : "dimgray")
            translate([0, 0, 1.5])
                draw_text2d(str("f", $ps_face_idx), size = 1.45);

        if ($ps_face_idx == HEX_FACE_IDX)
            color("tomato", 0.45)
                linear_extrude(height = FACE_THK, center = true)
                    ps_polygon($ps_face_pts2d, MODE);
    }

    draw_panel_label("anti-tet context");
}

/**
 * Module: Draw arrangement spans and crossing nodes for the hex face.
 * Params: none
 * Returns: none
 */
module draw_arrangement_panel() {
    place_on_faces(P, IR) {
        if ($ps_face_idx == HEX_FACE_IDX) {
            arr = ps_face_arrangement($ps_face_pts3d_local);
            nodes = arr[2];
            spans = arr[3];

            color("gainsboro", 0.22)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            for (si = [0:1:len(spans)-1]) {
                span = spans[si];
                is_long = span[3] == LONG_SOURCE_EDGE_IDX;
                color(is_long ? "red" : "slategray")
                    draw_segment_stroke(span[0], r = is_long ? LINE_R * 0.8 : LINE_R * 0.45);
            }

            for (ni = [0:1:len(nodes)-1]) {
                node = nodes[ni];
                color(node[1] == "crossing" ? "deepskyblue" : "black")
                    translate([node[0][0], node[0][1], FACE_THK])
                        sphere(r = NODE_R, $fn = 12);
            }

            draw_source_edge_labels($ps_face_pts2d);
            echo_hex_summary($ps_face_pts3d_local, $ps_face_idx, $ps_poly_faces_idx, $ps_poly_verts_local);
        }
    }

    draw_panel_label("arrangement");
}

/**
 * Module: Draw filled boundary source-edge groups in parent face coordinates.
 * Params: none
 * Returns: none
 */
module draw_source_edges_panel() {
    place_on_faces(P, IR) {
        if ($ps_face_idx == HEX_FACE_IDX) {
            color("gainsboro", 0.16)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            place_on_face_filled_boundary_source_edges(mode = MODE, coords = "parent") {
                is_long = $ps_boundary_source_edge_idx == LONG_SOURCE_EDGE_IDX;

                color(is_long ? "red" : example_color($ps_boundary_source_edge_idx), 0.34)
                    draw_segment_stroke($ps_boundary_source_edge_segment2d_local, r = is_long ? LINE_R * 1.0 : LINE_R * 0.7);

                for (si = [0:1:$ps_boundary_source_edge_span_count-1]) {
                    span_seg = $ps_boundary_source_edge_span_segments2d_local[si];
                    mid = ps_segment_midpoint2d(span_seg);
                    side = $ps_boundary_source_edge_sides[si];
                    d = span_seg[1] - span_seg[0];
                    len = norm(d);
                    left = (len <= 1e-9) ? [0, 1] : [-d[1] / len, d[0] / len];
                    tip = mid + ((side >= 0) ? 1 : -1) * 4.8 * left;

                    color(is_long ? "red" : example_color($ps_boundary_source_edge_idx))
                        draw_segment_stroke(span_seg, r = LINE_R * 0.48, z = 0.35);

                    color(side >= 0 ? "limegreen" : "crimson")
                        translate([0, 0, 1.2 * FACE_THK]) {
                            draw_segment_stroke([mid, tip], r = LINE_R * 0.28);
                            translate([tip[0], tip[1], 0])
                                sphere(r = 0.72, $fn = 12);
                        }
                }
            }

            draw_source_edge_labels($ps_face_pts2d);
        }
    }

    draw_panel_label("source edges");
}

/**
 * Module: Draw boundary-span frames and adjacent-face direction markers.
 * Params: none
 * Returns: none
 */
module draw_boundary_span_panel() {
    place_on_faces(P, IR) {
        if ($ps_face_idx == HEX_FACE_IDX) {
            color("gainsboro", 0.14)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            place_on_face_boundary_spans(mode = MODE) {
                is_long = $ps_boundary_span_source_edge_idx == LONG_SOURCE_EDGE_IDX;
                is_short = $ps_boundary_span_source_edge_idx == SHORT_SOURCE_EDGE_IDX;

                color(is_long ? "red" : (is_short ? "navy" : "black"))
                    cube([$ps_boundary_span_len, is_long ? 0.95 : 0.55, 0.55], center = true);

                color(is_long ? "deepskyblue" : "mediumseagreen")
                    draw_direction3d($ps_boundary_span_adj_face_dir_span_local);

                color("white")
                    translate([0, 0, 1.1])
                        draw_text2d(
                            str(
                                "s", $ps_boundary_span_idx,
                                "/se", $ps_boundary_span_source_edge_idx,
                                "/f", $ps_boundary_span_adj_face_idx
                            ),
                            size = 0.88
                        );
            }
        }
    }

    draw_panel_label("span directions");
}

/**
 * Module: Draw the positive anti-interference volume derived from the filled boundary.
 * Params: none
 * Returns: none
 */
module draw_volume_panel() {
    place_on_faces(P, IR) {
        if ($ps_face_idx == HEX_FACE_IDX) {
            color("gainsboro", 0.16)
                translate([0, 0, -FACE_THK])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, MODE);

            color("deepskyblue", 0.34)
                ps_face_anti_interference_volume(
                    VOL_Z_MIN,
                    VOL_Z_MAX,
                    mode = MODE,
                    max_project = MAX_PROJECT
                );

            color("black")
                place_on_face_boundary_spans(mode = MODE)
                    cube([$ps_boundary_span_len, LINE_R * 0.75, LINE_R * 0.75], center = true);
        }
    }

    draw_panel_label("positive volume");
}

translate([-2 * PANEL_X, 0, 0])
    draw_context_panel();

translate([-1 * PANEL_X, 0, 0])
    draw_arrangement_panel();

translate([0 * PANEL_X, 0, 0])
    draw_source_edges_panel();

translate([1 * PANEL_X, 0, 0])
    draw_boundary_span_panel();

translate([2 * PANEL_X, 0, 0])
    draw_volume_panel();
