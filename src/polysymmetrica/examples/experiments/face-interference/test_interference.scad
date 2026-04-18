// Face-only anti-interference probe.
// Use this to inspect the face-side interference volume and its cutters in isolation, without
// local edge clearance or the composed proxy carve path.

use <../../../core/funcs.scad>
use <../../../core/segments.scad>
use <../../../core/placement.scad>
use <../../../core/prisms.scad>
use <../../../core/face_regions.scad>

SC = 1;
IR = 20 * SC;

p = poly_antiprism(n = 7, p = 3, angle = 15);

// Keep the default probe on an ordinary convex face. Switch to `0` later when
// testing star-face interference in isolation.
FACE_IDX = 0;

FACE_T = 1.6 * SC;
BASE_Z = -FACE_T / 4 - 1;
PILLOW_THK = 0.4 * SC+0;
FACE_Z0 = BASE_Z;
FACE_Z1 = BASE_Z + FACE_T + PILLOW_THK;

COLUMN_SEPARATION = 60 * SC;
SHOW_CUTTERS = false;

module raw_face_volume() {
    translate([0, 0, FACE_Z0])
        linear_extrude(height = FACE_Z1 - FACE_Z0)
            union() {
                for (cell = ps_face_filled_cells($ps_face_pts3d_local, 1e-8))
                    ps_polygon(points = cell[0]);
            }
}

module show_probe() {
    place_on_faces(p, inter_radius = IR, indices = [FACE_IDX]) {
//        translate([-COLUMN_SEPARATION / 2, 0, 0])
//            color([0.72, 0.72, 0.77, 0.45])
//                raw_face_volume();
//
//        translate([COLUMN_SEPARATION / 2, 0, 0]) {
            color("deepskyblue")
                ps_face_interference_volume_ctx(FACE_Z0, FACE_Z1);

            if (SHOW_CUTTERS)
                color([1, 0.2, 0.2, 0.25])
                    ps_face_interference_cutters_ctx(FACE_Z0, FACE_Z1);
//        }
    }
}

show_probe();
