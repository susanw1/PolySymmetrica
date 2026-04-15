// Low-level boundary extraction probe for self-crossing/nonconvex faces.
// Use this to verify the true filled perimeter, source-edge spans, and inward-side markers
// before debugging interference or local edge-clearance behavior.

use <../../../core/funcs.scad>
use <../../../core/segments.scad>
use <../../../core/placement.scad>
use <../../../core/prisms.scad>

SC = 1;
IR = 18 * SC;

p = poly_antiprism(5, 2);

SEG_BAR_W = 0.8 * SC;
SEG_BAR_H = 0.8 * SC;
FACE_T = 0.6 * SC;
FACE_SEPARATION = 90 * SC;

module raw_face_plate() {
    linear_extrude(height = FACE_T)
        union() {
            for (cell = ps_face_filled_cells($ps_face_pts3d_local, 1e-8))
                ps_polygon(points = cell[0]);
        }
}

module boundary_seg_marker() {
    color("gold")
        translate([0, 0, FACE_T + SEG_BAR_H / 2])
            cube([$ps_face_boundary_seg_len, SEG_BAR_W, SEG_BAR_H], center = true);

    color($ps_face_boundary_seg_inward_is_positive_ey ? "limegreen" : "tomato")
        translate([0, 1.4 * SC * ($ps_face_boundary_seg_inward_is_positive_ey ? 1 : -1), FACE_T + SEG_BAR_H / 2])
            cube([1.6 * SC, 0.35 * SC, SEG_BAR_H], center = true);
}

module show_face(face_idx, global_shift = [0, 0, 0]) {
    translate(global_shift)
        place_on_faces(p, inter_radius = IR, indices = [face_idx]) {
            color([0.7, 0.7, 0.75, 0.45])
                raw_face_plate();

            place_on_face_filled_boundary_segments()
                boundary_seg_marker();
        }
}

show_face(0, [-FACE_SEPARATION / 2, 0, 0]);
show_face(2, [FACE_SEPARATION / 2, 0, 0]);
