use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/prisms.scad>
use <../../core/segments.scad>
use <../../core/truncation.scad>
use <../../models/dodecahedron.scad>
use <../../models/tetrahedron.scad>

// Boundary-model demo surface.
// Shows the reconstructed filled boundary loops/spans for:
// - one star face on a star antiprism
// - one self-crossing hex face on an antitruncated tetrahedron
// - one plain face on a dodecahedron as a control

IR = 30;
FACE_THK = 0.35;
LINE_R = 0.9;
TXT_H = 0.4;
TXT_S = 3.6;
PANEL_X = 135;
PANEL_Y = 118;
ROW_LABEL_X = -2.35 * PANEL_X;
RAW_FILL_Z = -1.2 * FACE_THK;
BOUNDARY_FILL_Z = 0.15 * FACE_THK;
SPAN_MARKER_Z = 1.2 * FACE_THK;
SPAN_ARROW_LEN = 7;
DIH_RAY_LEN = 10;
FACE_LABEL_Z = 0.8;
SOURCE_EDGE_MARKER_Z = 1.5 * FACE_THK;
SOURCE_SIDE_LEN = 5.2;

STAR_POLY = poly_antiprism(5, 2);
STAR_FACE_IDX = 1;

ANTI_POLY = poly_truncate(tetrahedron(), t = -0.5);
ANTI_FACE_IDX = 0;

DODECA_POLY = dodecahedron();
DODECA_FACE_IDX = 0;

/**
 * Function: Pick a readable named color for one boundary loop.
 * Params: i (loop index)
 * Returns: a named color string
 */
function boundary_loop_color(i) =
    (i % 6 == 0) ? "tomato" :
    (i % 6 == 1) ? "deepskyblue" :
    (i % 6 == 2) ? "gold" :
    (i % 6 == 3) ? "mediumseagreen" :
    (i % 6 == 4) ? "orchid" :
    "darkorange";

/**
 * Function: Pick a readable named color for one source edge.
 * Params: i (source edge index)
 * Returns: a named color string
 */
function source_edge_color(i) =
    (i % 8 == 0) ? "tomato" :
    (i % 8 == 1) ? "deepskyblue" :
    (i % 8 == 2) ? "gold" :
    (i % 8 == 3) ? "mediumseagreen" :
    (i % 8 == 4) ? "orchid" :
    (i % 8 == 5) ? "darkorange" :
    (i % 8 == 6) ? "turquoise" :
    "sienna";

/**
 * Module: Draw a world-space label below one panel.
 * Params: s (label string)
 * Returns: none
 */
module draw_panel_label(s) {
    translate([0, -58, -34])
        linear_extrude(height = TXT_H)
            text(s, size = TXT_S, halign = "center", valign = "center");
}

/**
 * Module: Draw a world-space label for one row of panels.
 * Params: s (row label string)
 * Returns: none
 */
module draw_row_label(s) {
    translate([ROW_LABEL_X, 0, -34])
        rotate([0, 0, 90])
            linear_extrude(height = TXT_H)
                text(s, size = TXT_S + 0.4, halign = "center", valign = "center");
}

/**
 * Module: Draw a thin local face fill from a point loop.
 * Params: pts2d (2D polygon loop)
 * Returns: none
 */
module draw_polygon(pts2d) {
    linear_extrude(height = FACE_THK, center = true)
        polygon(points = pts2d);
}

/**
 * Module: Draw a stroked local 2D segment as a thin extruded ribbon.
 * Params: seg2d (2D segment), r (stroke radius)
 * Returns: none
 */
module draw_local_segment_stroke(seg2d, r = LINE_R) {
    linear_extrude(height = FACE_THK * 0.9, center = true)
        hull() {
            translate(seg2d[0]) circle(r = r, $fn = 20);
            translate(seg2d[1]) circle(r = r, $fn = 20);
        }
}

/**
 * Module: Draw a local span-side marker in the current boundary-span frame.
 * Params: filled_side (`+1` or `-1`), len (marker length)
 * Returns: none
 */
module draw_span_side_marker(filled_side, len = SPAN_ARROW_LEN) {
    dir = (filled_side >= 0) ? 1 : -1;

    color((filled_side >= 0) ? "limegreen" : "crimson")
        union() {
            translate([0, dir * len / 2, SPAN_MARKER_Z])
                cube([0.55, len, 0.55], center = true);

            translate([0, dir * len, SPAN_MARKER_Z])
                sphere(r = 0.95, $fn = 12);
        }
}

/**
 * Module: Draw the adjacent-face plane direction in the current span frame.
 * Params: dir_span_local (adjacent-face direction `[x,y,z]` in span-local coords)
 * Returns: none
 */
module draw_dihedral_direction(dir_span_local) {
    dir = v_norm(dir_span_local);
    axis = [-dir[1], dir[0], 0];
    ang = acos(max(-1, min(1, dir[2])));

    if (norm(dir) > 1e-9)
        color("deepskyblue")
            rotate(a = ang, v = (norm(axis) <= 1e-9) ? [1, 0, 0] : axis)
                union() {
                    translate([0, 0, DIH_RAY_LEN / 2])
                        cylinder(h = DIH_RAY_LEN, r = 0.45, center = true, $fn = 14);
                    translate([0, 0, DIH_RAY_LEN])
                        sphere(r = 0.95, $fn = 12);
                }
}

/**
 * Module: Draw source-edge labels on the selected face.
 * Params: pts2d (face loop in local 2D)
 * Returns: none
 */
module draw_source_edge_labels(pts2d) {
    for (ei = [0:1:len(pts2d)-1]) {
        seg2d = [pts2d[ei], pts2d[(ei+1)%len(pts2d)]];
        mid = ps_segment_midpoint2d(seg2d);
        inward = ps_centroid2d(pts2d) - mid;
        offset = (norm(inward) <= 1e-9) ? [0, 0] : inward / norm(inward);

        color("dimgray")
            translate([mid[0] + 2.8 * offset[0], mid[1] + 2.8 * offset[1], FACE_THK / 2 + 0.02])
                linear_extrude(height = TXT_H)
                    text(str("se", ei), size = 1.5, halign = "center", valign = "center");
    }
}

/**
 * Module: Draw a light poly wireframe context for the demo poly.
 * Params: poly (poly descriptor), ir (placement scale)
 * Returns: none
 */
module draw_wireframe(poly, ir = IR, show_face_labels = false) {
    color("silver")
        place_on_edges(poly, ir)
            cube([$ps_edge_len, 0.8, 0.8], center = true);

    color("gold")
        place_on_vertices(poly, ir)
            sphere(1.2, $fn = 12);

    if (show_face_labels)
        color("black")
            place_on_faces(poly, ir)
                translate([0, 0, FACE_LABEL_Z])
                    linear_extrude(height = TXT_H)
                        text(str("f", $ps_face_idx), size = 1.9, halign = "center", valign = "center");
}

/**
 * Module: Draw one boundary-model panel for one selected face.
 * Params: poly (poly descriptor), face_idx (target face index), mode (`"nonzero"` or `"evenodd"`), label_s (panel label)
 * Returns: none
 */
module draw_panel_boundary_model(poly, face_idx, mode, label_s) {
    draw_wireframe(poly, IR);

    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            bm = ps_face_boundary_model($ps_face_pts3d_local, mode = mode);
            loops = bm[2];
            spans = bm[3];

            color("gainsboro", 0.30)
                translate([0, 0, RAW_FILL_Z])
                    draw_polygon($ps_face_pts2d);

            for (li = [0:1:len(loops)-1]) {
                loop = loops[li];

                color(boundary_loop_color(li), 0.55)
                    translate([0, 0, BOUNDARY_FILL_Z])
                        draw_polygon(loop[0]);

                color("white") {
                    centre = ps_centroid2d(loop[0]);
                    translate([centre[0], centre[1], FACE_THK / 2 + 0.02])
                        linear_extrude(height = TXT_H)
                            text(str("L", li), size = 2.2, halign = "center", valign = "center");
                }
            }

            for (si = [0:1:len(spans)-1]) {
                span = spans[si];

                color("black")
                    draw_local_segment_stroke(span[0], r = LINE_R * 0.6);

                color("black") {
                    mid = ps_segment_midpoint2d(span[0]);
                    translate([mid[0], mid[1], FACE_THK / 2 + 0.02])
                        linear_extrude(height = TXT_H)
                            text(str(si), size = 1.5, halign = "center", valign = "center");
                }
            }
        }
    }

    draw_panel_label(label_s);
}

/**
 * Module: Draw one span-site panel for one selected face.
 * Params: poly (poly descriptor), face_idx (target face index), mode (`"nonzero"` or `"evenodd"`), label_s (panel label)
 * Returns: none
 */
module draw_panel_boundary_spans(poly, face_idx, mode, label_s) {
    draw_wireframe(poly, IR, true);

    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            bm = ps_face_boundary_model($ps_face_pts3d_local, mode = mode);
            loops = bm[2];

            color("gainsboro", 0.25)
                translate([0, 0, RAW_FILL_Z])
                    draw_polygon($ps_face_pts2d);

            for (li = [0:1:len(loops)-1]) {
                loop = loops[li];

                color(boundary_loop_color(li), 0.30)
                    translate([0, 0, BOUNDARY_FILL_Z])
                        draw_polygon(loop[0]);
            }

            draw_source_edge_labels($ps_face_pts2d);

            place_on_face_boundary_spans(mode = mode) {
                color("black")
                    cube([$ps_boundary_span_len, 0.8, 0.3], center = true);

                draw_span_side_marker($ps_boundary_span_filled_side);

                translate([0, 0, SPAN_MARKER_Z])
                    // Demo-only display choice:
                    // the core metadata is anchored to the original source-edge direction,
                    // then we flip per span so the rod reads consistently against the
                    // displayed boundary-span parameter direction.
                    draw_dihedral_direction(
                        (($ps_boundary_span_source_t1 > $ps_boundary_span_source_t0) ? -1 : 1)
                        * $ps_boundary_span_adj_face_dir_span_local
                    );

                color("white")
                    translate([0, 0, SPAN_MARKER_Z + DIH_RAY_LEN + 0.4])
                        linear_extrude(height = TXT_H)
                            text(
                                str("se", $ps_boundary_span_source_edge_idx, "/f", $ps_boundary_span_adj_face_idx),
                                size = 1.4,
                                halign = "center",
                                valign = "center"
                            );
            }
        }
    }

    draw_panel_label(label_s);
}

/**
 * Module: Draw one source-edge grouped boundary panel for one selected face.
 * Params: poly (poly descriptor), face_idx (target face index), mode (`"nonzero"` or `"evenodd"`), label_s (panel label)
 * Returns: none
 */
module draw_panel_boundary_source_edges(poly, face_idx, mode, label_s) {
    draw_wireframe(poly, IR, true);

    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            color("gainsboro", 0.20)
                translate([0, 0, RAW_FILL_Z])
                    draw_polygon($ps_face_pts2d);

            draw_source_edge_labels($ps_face_pts2d);

            place_on_face_filled_boundary_source_edges(mode = mode) {
                echo(str(label_s, "draw_panel_boundary_source_edges::", $ps_boundary_source_edge_idx));

                color(source_edge_color($ps_boundary_source_edge_idx), 0.25)
                    translate([0, 0, SOURCE_EDGE_MARKER_Z - 0.35])
                        cube([$ps_boundary_source_edge_len, LINE_R * 1.7, FACE_THK], center = true);

                for (si = [0:1:$ps_boundary_source_edge_span_count-1]) {
                    t_range = $ps_boundary_source_edge_span_t_ranges_local[si];
                    span_filled_side = $ps_boundary_source_edge_span_sides_local[si];
                    side_dir = (span_filled_side >= 0) ? 1 : -1;
                    x0 = (t_range[0] - 0.5) * $ps_boundary_source_edge_len;
                    x1 = (t_range[1] - 0.5) * $ps_boundary_source_edge_len;
                    x_mid = (x0 + x1) / 2;
                    span_len = abs(x1 - x0);

                    color(source_edge_color($ps_boundary_source_edge_idx))
                        translate([x_mid, 0, SOURCE_EDGE_MARKER_Z])
                            cube([span_len, LINE_R * 1.1, FACE_THK], center = true);

                    color(span_filled_side >= 0 ? "limegreen" : "crimson")
                        translate([0, 0, SOURCE_EDGE_MARKER_Z + 0.8])
                            union() {
                                translate([x_mid, side_dir * SOURCE_SIDE_LEN / 2, 0])
                                    cube([LINE_R * 0.45, SOURCE_SIDE_LEN, LINE_R * 0.45], center = true);
                                translate([x_mid, side_dir * SOURCE_SIDE_LEN, 0])
                                    sphere(r = 0.9, $fn = 12);
                            }
                }

                color("white")
                    translate([0, SOURCE_SIDE_LEN + 1.2, SOURCE_EDGE_MARKER_Z + 2.1])
                        linear_extrude(height = TXT_H)
                            text(
                                str(
                                    "se",
                                    $ps_boundary_source_edge_idx,
                                    ":",
                                    $ps_boundary_source_edge_span_count,
                                    " spans"
                                ),
                                size = 1.45,
                                halign = "center",
                                valign = "center"
                            );
            }
        }
    }

    draw_panel_label(label_s);
}

/**
 * Module: Draw the four boundary-model panels for one selected face.
 * Params: poly (poly descriptor), face_idx (target face index), row_s (row label), y (row offset)
 * Returns: none
 */
module draw_boundary_row(poly, face_idx, row_s, y) {
    translate([0, y, 0])
        draw_row_label(row_s);

    translate([-1.5 * PANEL_X, y, 0])
        draw_panel_boundary_model(poly, face_idx, "nonzero", "boundary nonzero");

    translate([-0.5 * PANEL_X, y, 0])
        draw_panel_boundary_model(poly, face_idx, "evenodd", "boundary evenodd");

    translate([0.5 * PANEL_X, y, 0])
        draw_panel_boundary_spans(poly, face_idx, "nonzero", "boundary spans");

    translate([1.5 * PANEL_X, y, 0])
        draw_panel_boundary_source_edges(poly, face_idx, "nonzero", "source edges");
}

draw_boundary_row(STAR_POLY, STAR_FACE_IDX, "star antiprism", PANEL_Y);
draw_boundary_row(ANTI_POLY, ANTI_FACE_IDX, "anti-tet hex", 0);
draw_boundary_row(DODECA_POLY, DODECA_FACE_IDX, "dodeca control", -PANEL_Y);
