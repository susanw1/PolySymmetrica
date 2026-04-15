// Primary edge-interference probe for one target source edge.
// Use this to inspect the raw strip, phase 1, phase 2, and crossing-lobe overlap without
// the extra noise of the full proxy/carve stack.

use <../../../core/funcs.scad>
use <../../../core/segments.scad>
use <../../../core/placement.scad>
use <../../../core/prisms.scad>
use <../../../core/face_regions.scad>
use <../../../core/proxy_interaction.scad>

SC = 1;
IR = 20 * SC;

p = poly_antiprism(n = 9, p = 2, angle = 15);
echo (p);
FACE_IDX = 0;
SEG_IDX = 0;
LOBE_IDX = undef;

FACE_T = 1.6 * SC;
BASE_Z = -FACE_T / 4 - 1;
PILLOW_THK = 0.4 * SC;
FACE_Z0 = BASE_Z;
FACE_Z1 = BASE_Z + FACE_T + PILLOW_THK;

// Probe views. Turn on one main colored layer at a time for clarity:
// - SHOW_FACE_CONTEXT: faint finished face result for context
// - SHOW_RAW_LOBES: full crossing-lobe bodies in the target-edge frame
// - SHOW_LOCALIZED_LOBES: lobe overlap with the phase-1 cutter
// - SHOW_RAW_STRIP: the uncropped one-sided edge strip
// - SHOW_EDGE_LINE: the exact source edge line at local y = 0
// - SHOW_PHASE1: target strip after own anti-interference only
// - SHOW_PHASE2: target strip after crossing-lobe subtraction
// - SHOW_PHASE1_REMOVED_BY_PHASE2: what phase 2 removes from phase 1
SHOW_FACE_CONTEXT = true;
SHOW_RAW_LOBES = false;
SHOW_LOCALIZED_LOBES = true;
SHOW_RAW_STRIP = false;
SHOW_EDGE_LINE = true;
SHOW_PHASE1 = false;
SHOW_PHASE2 = false;
SHOW_PHASE1_REMOVED_BY_PHASE2 = false;

PROBE_CUTTER_W = 1.5 * SC;
PROBE_CUTTER_H = 10 * SC;

CROP_TO_WINDOW = false;
CROP_CENTER = [0, 0, FACE_Z0];
CROP_SIZE = [18 * SC, 18 * SC, 8 * SC];
SHOW_EDGE_LABELS = true;
SHOW_VERTEX_MARKERS = true;

module maybe_crop() {
    if (CROP_TO_WINDOW)
        intersection() {
            children();
            translate(CROP_CENTER)
                cube(CROP_SIZE, center = true);
        }
    else
        children();
}

module probe_cutter() {
    y_shift = $ps_face_boundary_seg_inward_is_positive_ey ? -PROBE_CUTTER_W / 2 : PROBE_CUTTER_W / 2;

    translate([0, y_shift, 0])
        cube([$ps_edge_len, PROBE_CUTTER_W, PROBE_CUTTER_H], center = true);
}

module source_edge_line() {
    cube([$ps_edge_len, 0.08 * SC, 0.08 * SC], center = true);
}

module show_probe() {
    place_on_faces(p, inter_radius = IR, indices = [FACE_IDX]) {
        maybe_crop() {
            if (SHOW_FACE_CONTEXT)
                color("lightgrey", 0.9)
                    ps_face_interference_volume_ctx(FACE_Z0, FACE_Z1);

            place_on_face_filled_boundary_source_edges() {
                if ($ps_face_boundary_seg_source_edge_idx == SEG_IDX) {
                    if (SHOW_RAW_LOBES)
                        color("darkblue", 1)
                            _ps_proxy_emit_crossing_local_boundary_lobe_bodies_raw_in_current_edge_frame(
                                [FACE_Z0, FACE_Z1],
                                undef,
                                is_undef(LOBE_IDX) ? undef : [LOBE_IDX],
                                $ps_face_boundary_seg_source_edge_idx
                            );

                    if (SHOW_LOCALIZED_LOBES)
                        color("gold", 1)
                            intersection() {
                                _ps_proxy_emit_crossing_local_boundary_lobe_bodies_raw_in_current_edge_frame(
                                    [FACE_Z0, FACE_Z1],
                                    undef,
                                    is_undef(LOBE_IDX) ? undef : [LOBE_IDX],
                                    $ps_face_boundary_seg_source_edge_idx
                                );

                                ps_clip_local_boundary_edge_by_own_interference_ctx(
                                    undef,
                                    undef,
                                    [FACE_Z0, FACE_Z1]
                                )
                                    probe_cutter();
                            }

                    if (SHOW_RAW_STRIP)
                        color("tomato", 0.32)
                            probe_cutter();

                    if (SHOW_EDGE_LINE)
                        color("black")
                            source_edge_line();

                    if (SHOW_PHASE1)
                        color("green")
                            ps_clip_local_boundary_edge_by_own_interference_ctx(
                                undef,
                                undef,
                                [FACE_Z0, FACE_Z1]
                            )
                                probe_cutter();

                    if (SHOW_PHASE2)
                        color("deepskyblue",1)
                            ps_clip_local_boundary_edge_by_crossing_edges_ctx(
                                undef,
                                undef,
                                [FACE_Z0, FACE_Z1],
                                undef,
                                is_undef(LOBE_IDX) ? undef : [LOBE_IDX]
                            )
                                probe_cutter();

                    if (SHOW_PHASE1_REMOVED_BY_PHASE2)
                        color("orange", 1)
                            difference() {
                                ps_clip_local_boundary_edge_by_own_interference_ctx(
                                    undef,
                                    undef,
                                    [FACE_Z0, FACE_Z1]
                                )
                                    probe_cutter();

                                ps_clip_local_boundary_edge_by_crossing_edges_ctx(
                                    undef,
                                    undef,
                                    [FACE_Z0, FACE_Z1],
                                    undef,
                                    is_undef(LOBE_IDX) ? undef : [LOBE_IDX]
                                )
                                    probe_cutter();
                            }

                }
            }
        }
    }
    if (SHOW_EDGE_LABELS)
        color("darkgray") place_on_edges(p, inter_radius = IR, indices = undef) {
            linear_extrude(1) text(str($ps_edge_idx), size=1, halign="center",valign="center");
        }
    if (SHOW_VERTEX_MARKERS)
        color("darkgray") place_on_vertices(p, inter_radius = IR, indices = undef) {
            cube(0.4, center=true);
        }
}

show_probe();
