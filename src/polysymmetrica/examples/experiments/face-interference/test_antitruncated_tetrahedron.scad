// Antitruncated tetrahedron probe built from the existing truncation operator.
// Use this to inspect face, edge, and vertex placement on a self-intersecting but
// much better-behaved shape than the ad hoc crossing-hexagon shell.

use <../../../core/funcs.scad>
use <../../../core/placement.scad>
use <../../../core/segments.scad>
use <../../../core/validate.scad>
use <../../../core/truncation.scad>
use <../../../core/classify.scad>
use <../../../models/tetrahedron.scad>

SC = 1;
DISPLAY_SCALE = 18 * SC;

TARGET_FACE_IDX = 0;
FACE_SHEET_T = 0.10 * SC;
WALK_STROKE_W = 0.18 * SC;
LABEL_T = 0.08 * SC;
VERT_MARKER = 0.35 * SC;

SHOW_SHELL = true;
SHOW_FACE_LABELS = true;
SHOW_EDGE_LABELS = true;
SHOW_VERTEX_MARKERS = true;
SHOW_FACE_NORMAL = true;
SHOW_EDGE_BISECTORS = true;

/**
 * Function: Build the experiment poly using the existing truncation operator.
 * Params: none.
 * Returns: `poly_truncate(tetrahedron(), t = -0.5)`.
 * Limitations: This shape is expected to be `star_ok` rather than `closed` valid.
 */
function antitruncated_tetrahedron() =
    poly_truncate(tetrahedron(), t = -0.5);

/**
 * Function: Return a readable face fill color for one face family.
 * Params: family_id (face family id from `poly_classify(...)`).
 * Returns: Named color string.
 */
function anti_tet_face_color(family_id) =
    (family_id == 0) ? "khaki" : "crimson";

/**
 * Function: Return a readable edge marker color for one edge family.
 * Params: family_id (edge family id from `poly_classify(...)`).
 * Returns: Named color string.
 */
function anti_tet_edge_color(family_id) =
    (family_id == 0) ? "darkolivegreen" : "teal";

/**
 * Module: Draw a constant-width 2D stroke for one segment.
 * Params: a,b (segment endpoints), w (stroke width).
 * Returns: No value; emits one 2D strip polygon.
 * Limitations: Display helper only; degenerate segments are ignored.
 */
module stroke_segment2d(a, b, w = WALK_STROKE_W) {
    d = b - a;
    L = norm(d);

    if (L > 1e-8) {
        ex = d / L;
        ey = [-ex[1], ex[0]];

        polygon(points = [
            a + ey * (w / 2),
            b + ey * (w / 2),
            b - ey * (w / 2),
            a - ey * (w / 2)
        ]);
    }
}

/**
 * Module: Draw a stroked version of a closed 2D face walk.
 * Params: pts (closed walk vertices), w (stroke width).
 * Returns: No value; emits one strip segment per edge.
 */
module stroke_poly2d(pts, w = WALK_STROKE_W) {
    for (i = [0:1:len(pts)-1])
        stroke_segment2d(pts[i], pts[(i + 1) % len(pts)], w);
}

/**
 * Module: Render filled face sheets and face labels for the experiment poly.
 * Params: poly (poly descriptor), cls (precomputed classify result).
 * Returns: No value; emits face-local sheets and optional labels.
 * Limitations: Uses `ps_polygon(..., mode="nonzero")` for all faces so self-crossing
 *   hexagons render as filled regions rather than raw `polygon(...)` walks.
 */
module render_shell_faces(poly, cls) {
    place_on_faces(poly, edge_len = DISPLAY_SCALE, classify = cls) {
        color(anti_tet_face_color($ps_face_family_id), 0.92)
            linear_extrude(height = FACE_SHEET_T, center = true)
                ps_polygon(points = $ps_face_pts2d, mode = "nonzero");

        color("black")
            linear_extrude(height = FACE_SHEET_T * 1.2, center = true)
                stroke_poly2d($ps_face_pts2d);

        if (SHOW_FACE_LABELS)
            color("black")
                translate([0, 0, 0.7 * SC])
                    linear_extrude(height = LABEL_T, center = true)
                        text(str("f", $ps_face_idx), size = 0.9 * SC, halign = "center", valign = "center");
    }
}

/**
 * Module: Render one simple marker along the target face local +Z direction.
 * Params: poly (poly descriptor), cls (precomputed classify result).
 * Returns: No value; emits one bar centered on the target face frame.
 */
module render_target_face_normal(poly, cls) {
    if (SHOW_FACE_NORMAL)
        color("navy")
            place_on_faces(poly, edge_len = DISPLAY_SCALE, classify = cls, indices = [TARGET_FACE_IDX])
                translate([0, 0, 1.2 * SC])
                    cube([0.18 * SC, 0.18 * SC, 2.4 * SC], center = true);
}

/**
 * Module: Render edge-frame markers and labels for all edges.
 * Params: poly (poly descriptor), cls (precomputed classify result).
 * Returns: No value; emits one marker stack per edge.
 */
module render_edge_markers(poly, cls) {
    if (SHOW_EDGE_BISECTORS)
        place_on_edges(poly, edge_len = DISPLAY_SCALE, classify = cls) {
            color(anti_tet_edge_color($ps_edge_family_id))
                translate([0, 0, 0.85 * SC])
                    cube([0.22 * SC, 0.22 * SC, 1.7 * SC], center = true);
        }

    if (SHOW_EDGE_LABELS)
        color("black")
            place_on_edges(poly, edge_len = DISPLAY_SCALE, classify = cls) {
                translate([0, 0, 1.95 * SC])
                    linear_extrude(height = LABEL_T, center = true)
                        text(str("e", $ps_edge_idx), size = 0.8 * SC, halign = "center", valign = "center");
            }
}

/**
 * Module: Render vertex markers for all vertices.
 * Params: poly (poly descriptor), cls (precomputed classify result).
 * Returns: No value; emits one cube per vertex.
 */
module render_vertex_markers(poly, cls) {
    if (SHOW_VERTEX_MARKERS)
        color("dimgray")
            place_on_vertices(poly, edge_len = DISPLAY_SCALE, classify = cls)
                cube(VERT_MARKER, center = true);
}

/**
 * Module: Assemble the full visual probe scene.
 * Params: poly (poly descriptor), cls (precomputed classify result).
 * Returns: No value; emits the enabled shell and placement overlays.
 */
module show_scene(poly, cls) {
    if (SHOW_SHELL)
        render_shell_faces(poly, cls);

    render_target_face_normal(poly, cls);
    render_edge_markers(poly, cls);
    render_vertex_markers(poly, cls);
}

p = antitruncated_tetrahedron();
cls = poly_classify(p, 1, 1e-6, 1, false);

echo("antitruncated_tetrahedron closed valid:", poly_valid(p, "closed"));
echo("antitruncated_tetrahedron star_ok valid:", poly_valid(p, "star_ok"));
echo("antitruncated_tetrahedron counts:", [len(poly_verts(p)), len(poly_faces(p)), len(poly_edges(p))]);
echo("antitruncated_tetrahedron family_counts:", ps_classify_counts(cls));

show_scene(p, cls);
