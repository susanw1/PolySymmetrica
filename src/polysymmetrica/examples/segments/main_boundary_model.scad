use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/prisms.scad>
use <../../core/segments.scad>
use <../../core/truncation.scad>
use <../../models/tetrahedron.scad>

// Boundary-model demo surface.
// Shows the reconstructed filled boundary loops/spans for:
// - one star face on a star antiprism
// - one self-crossing hex face on an antitruncated tetrahedron

IR = 30;
FACE_THK = 0.35;
LINE_R = 0.9;
TXT_H = 0.4;
TXT_S = 3.6;
PANEL_X = 135;
PANEL_Y = 118;
RAW_FILL_Z = -1.2 * FACE_THK;
BOUNDARY_FILL_Z = 0.15 * FACE_THK;
SPAN_MARKER_Z = 1.2 * FACE_THK;
SPAN_ARROW_LEN = 7;
DIH_RAY_LEN = 10;
FACE_LABEL_Z = 0.8;

STAR_POLY = poly_antiprism(5, 2);
STAR_FACE_IDX = 1;

ANTI_POLY = poly_truncate(tetrahedron(), t = -0.5);
ANTI_FACE_IDX = 0;

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
 * Function: Compute the midpoint of a 2D segment.
 * Params: seg2d (`[[x0,y0],[x1,y1]]`)
 * Returns: midpoint `[x, y]`
 */
function segment_midpoint2d(seg2d) =
    [(seg2d[0][0] + seg2d[1][0]) / 2, (seg2d[0][1] + seg2d[1][1]) / 2];

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
        mid = segment_midpoint2d(seg2d);
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
                    mid = segment_midpoint2d(span[0]);
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

translate([-1 * PANEL_X,  0.5 * PANEL_Y, 0])
    draw_panel_boundary_model(STAR_POLY, STAR_FACE_IDX, "nonzero", "star boundary nonzero");

translate([0 * PANEL_X,  0.5 * PANEL_Y, 0])
    draw_panel_boundary_model(STAR_POLY, STAR_FACE_IDX, "evenodd", "star boundary evenodd");

translate([1 * PANEL_X,  0.5 * PANEL_Y, 0])
    draw_panel_boundary_spans(STAR_POLY, STAR_FACE_IDX, "nonzero", "star boundary spans");

translate([-1 * PANEL_X, -0.5 * PANEL_Y, 0])
    draw_panel_boundary_model(ANTI_POLY, ANTI_FACE_IDX, "nonzero", "anti-tet hex nonzero");

translate([0 * PANEL_X, -0.5 * PANEL_Y, 0])
    draw_panel_boundary_model(ANTI_POLY, ANTI_FACE_IDX, "evenodd", "anti-tet hex evenodd");

translate([1 * PANEL_X, -0.5 * PANEL_Y, 0])
    draw_panel_boundary_spans(ANTI_POLY, ANTI_FACE_IDX, "nonzero", "anti-tet hex spans");
