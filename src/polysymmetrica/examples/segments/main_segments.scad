use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/prisms.scad>
use <../../core/segments.scad>

// Segments demo surface.
// Shows the main face-local segment-analysis modules on one star antiprism:
// parity fill, winding fill, segmented cells, geometry cuts, and visible cells.

IR = 30;
FACE_THK = 0.35;
LINE_R = 1.35;
TXT_H = 0.18;
TXT_S = 3.8;
PANEL_X = 125;

P = poly_antiprism(5, 2);
FACES = poly_faces(P);
STAR_FACE_IDX = 1;
SIDE_FACE_IDX = 8;

/**
 * Function: Pick a strongly differentiated named color for cell visualization.
 * Params: i (cell index)
 * Returns: a named color string
 */
function cell_color(i) =
    (i % 6 == 0) ? "tomato" :
    (i % 6 == 1) ? "gold" :
    (i % 6 == 2) ? "mediumseagreen" :
    (i % 6 == 3) ? "deepskyblue" :
    (i % 6 == 4) ? "orchid" :
    "darkorange";

/**
 * Module: Draw a world-space label below one panel.
 * Params: s (label string)
 * Returns: none
 */
module draw_panel_label(s) {
    translate([0, -55, -34])
        linear_extrude(height = TXT_H)
            text(s, size = TXT_S, halign = "center", valign = "center");
}

/**
 * Module: Draw a thin local face fill from a point loop.
 * Params: pts2d (2D polygon loop), color_name (named color)
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
    linear_extrude(height = FACE_THK * 0.8, center = true)
        hull() {
            translate(seg2d[0]) circle(r = r, $fn = 20);
            translate(seg2d[1]) circle(r = r, $fn = 20);
        }
        translate(seg2d[0]) cylinder(r=0.3, h=5);
        translate(seg2d[1]) cylinder(r=0.2, h=6);
}

/**
 * Module: Draw a light poly wireframe context for the demo poly.
 * Params: poly (poly descriptor), ir (placement scale)
 * Returns: none
 */
module draw_wireframe(poly, ir = IR) {
    color("silver")
        place_on_edges(poly, ir)
            cube([$ps_edge_len, 0.8, 0.8], center = true);

    color("gold")
        place_on_vertices(poly, ir)
            sphere(1.2, $fn = 12);
}

/**
 * Module: Show `ps_polygon(...)` on the star face with specified fill rule.
 * Params: mode (`"evenodd"` or `"nonzero"`), label_s (panel label)
 * Returns: none
 */
module draw_panel_ps_polygon(mode, label_s) {
    draw_wireframe(P, IR);
    
    place_on_faces(P, IR) {
        if ($ps_face_idx == STAR_FACE_IDX) {
            color(mode == "evenodd" ? "lightsalmon" : "deepskyblue")
                linear_extrude(height = FACE_THK, center = true)
                    ps_polygon($ps_face_pts2d, mode = mode);
        }
    }
    draw_panel_label(label_s);
}

/**
 * Module: Show `place_on_face_segments(...)` on the star face.
 * Params: none
 * Returns: none
 */
module draw_panel_face_segments() {
    draw_wireframe(P, IR);
    
    place_on_faces(P, IR) {
        if ($ps_face_idx == STAR_FACE_IDX) {
            place_on_face_segments(mode = "nonzero") {
                // simply draw the polygon in supplied color
                color(cell_color($ps_seg_idx)) 
                    draw_polygon($ps_seg_pts2d);

                // draw a line along each segment
                color("black") 
                    for (seg2d = ps_cyclic_pairs($ps_seg_pts2d))
                        draw_local_segment_stroke(seg2d, r = LINE_R);

                // draw the segment index at the centroid of the segment
                color("white") {
                    centre = ps_centroid2d($ps_seg_pts2d);
                    translate([centre[0], centre[1], FACE_THK / 2 + 0.02])
                        linear_extrude(height = TXT_H)
                            text(str($ps_seg_idx), size = 2.2, halign = "center", valign = "center");
                }
            }
        }
    }
    draw_panel_label("place_on_face_segments");
}

/**
 * Module: Show `place_on_face_geom_cut_segments(...)` on one side face.
 * Params: none
 * Returns: none
 */
module draw_panel_geom_cuts() {
    draw_wireframe(P, IR);
    
    place_on_faces(P, IR) {
        if ($ps_face_idx == SIDE_FACE_IDX) {
            color("gainsboro") draw_polygon($ps_face_pts2d);
            
            place_on_face_geom_cut_segments(mode = "nonzero", filter_parent = true) {
                color("crimson") draw_local_segment_stroke($ps_face_cut_segment2d_local, r = LINE_R * 1.1);

                color("black")
                    for (j = [0:1:1])
                        translate([
                            $ps_face_cut_segment2d_local[j][0],
                            $ps_face_cut_segment2d_local[j][1],
                            FACE_THK / 2 + 0.02
                        ])
                            linear_extrude(height = TXT_H * 0.7)
                                text("x", size = 1.4, halign = "center", valign = "center");

                color("black") {
                    mid = ps_segment_midpoint2d($ps_face_cut_segment2d_local);
                    translate([mid[0], mid[1], FACE_THK / 2 + 0.02])
                        linear_extrude(height = TXT_H)
                            text(str($ps_face_cut_idx), size = 1.8, halign = "center", valign = "center");
                }
            }
        }
    }
    draw_panel_label("place_on_face_geom_cut_segments");
}

/**
 * Module: Show `place_on_face_visible_segments(...)` on one side face.
 * Params: none
 * Returns: none
 */
module draw_panel_visible_segments() {
    draw_wireframe(P, IR);
    
    place_on_faces(P, IR) {
        if ($ps_face_idx == SIDE_FACE_IDX) {
            color("gainsboro", 0.35)
                linear_extrude(height = FACE_THK * 0.4, center = true)
                    polygon(points = $ps_face_pts2d);

            place_on_face_visible_segments(mode = "nonzero", filter_parent = true) {
                color(cell_color($ps_vis_seg_idx))
                    linear_extrude(height = FACE_THK, center = true)
                        polygon(points = $ps_vis_seg_pts2d);

                color("green")
                    for (seg2d = ps_cyclic_pairs($ps_vis_seg_pts2d))
                        draw_local_segment_stroke(seg2d, r = LINE_R * 0.4);

                color("white")
                    translate([ps_centroid2d($ps_vis_seg_pts2d)[0], ps_centroid2d($ps_vis_seg_pts2d)[1], FACE_THK / 2 + 0.02])
                        linear_extrude(height = TXT_H)
                            text(str($ps_vis_seg_idx), size = 2.0, halign = "center", valign = "center");
            }
        }
    }
    draw_panel_label("place_on_face_visible_segments");
}

translate([-2 * PANEL_X, 0, 0]) draw_panel_ps_polygon("evenodd", "ps_polygon evenodd");
translate([-1 * PANEL_X, 0, 0]) draw_panel_ps_polygon("nonzero", "ps_polygon nonzero");
translate([ 0 * PANEL_X, 0, 0]) draw_panel_face_segments();
translate([ 1 * PANEL_X, 0, 0]) draw_panel_geom_cuts();
translate([ 2 * PANEL_X, 0, 0]) draw_panel_visible_segments();
