use <../../core/face_regions.scad>
use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/prisms.scad>
use <../../core/segments.scad>
use <../../core/truncation.scad>
use <../../models/tetrahedron.scad>

// Face anti-interference volume demo.
// Shows the positive face-local shell built from filled boundary spans for:
// - one star face on a pentagrammic antiprism
// - the self-crossing hex face on an antitruncated tetrahedron

IR = 30;
FACE_THK = 0.25;
VOL_Z_MIN = -2;
VOL_Z_MAX = 2;
MAX_PROJECT = 40;
LINE_R = 0.55;
TXT_H = 0.35;
TXT_S = 3.2;
PANEL_X = 130;

STAR_POLY = poly_antiprism(5, 2);
STAR_FACE_IDX = 0;

ANTI_POLY = poly_truncate(tetrahedron(), t = -0.5);
ANTI_FACE_IDX = 0;

/**
 * Module: Draw a world-space panel label.
 * Params: s (label string)
 * Returns: none
 */
module draw_panel_label(s) {
    translate([0, -54, -22])
        linear_extrude(height = TXT_H)
            text(s, size = TXT_S, halign = "center", valign = "center");
}

/**
 * Module: Draw a light poly edge/vertex context.
 * Params: poly (poly descriptor), ir (placement scale)
 * Returns: none
 */
module draw_wireframe(poly, ir = IR) {
    color("silver")
        place_on_edges(poly, ir)
            cube([$ps_edge_len, 0.65, 0.65], center = true);

    color("gold")
        place_on_vertices(poly, ir)
            sphere(1.05, $fn = 12);
}

/**
 * Module: Draw the selected face fill, boundary spans, and generated volume.
 * Params: poly (poly descriptor), face_idx (target face index), label_s (panel label)
 * Returns: none
 */
module draw_volume_panel(poly, face_idx, label_s) {
    draw_wireframe(poly, IR);

    place_on_faces(poly, IR) {
        if ($ps_face_idx == face_idx) {
            color("gainsboro", 0.30)
                translate([0, 0, -0.08])
                    linear_extrude(height = FACE_THK, center = true)
                        ps_polygon($ps_face_pts2d, mode = "nonzero");

            color("deepskyblue", 0.38)
                ps_face_anti_interference_volume_ctx(
                    VOL_Z_MIN,
                    VOL_Z_MAX,
                    mode = "nonzero",
                    max_project = MAX_PROJECT
                );

            color("black")
                place_on_face_boundary_spans(mode = "nonzero")
                    cube([$ps_boundary_span_len, LINE_R, LINE_R], center = true);

            shells = ps_face_anti_interference_shells(
                $ps_face_pts3d_local,
                $ps_face_idx,
                $ps_poly_faces_idx,
                $ps_poly_verts_local,
                $ps_face_neighbors_idx,
                $ps_face_dihedrals,
                VOL_Z_MIN,
                VOL_Z_MAX,
                "nonzero",
                MAX_PROJECT
            );

            color("navy")
                for (shell = shells)
                    for (pt = shell[0])
                        translate(pt)
                            sphere(0.75, $fn = 10);
        }
    }

    draw_panel_label(label_s);
}

translate([-0.5 * PANEL_X, 0, 0])
    draw_volume_panel(STAR_POLY, STAR_FACE_IDX, "star volume");

translate([0.5 * PANEL_X, 0, 0])
    draw_volume_panel(ANTI_POLY, ANTI_FACE_IDX, "anti-tet hex volume");
