use <../../core/funcs.scad>
use <../../core/face_regions.scad>
use <../../core/placement.scad>
use <../../core/prisms.scad>

// Minimum face radius before adding the pillow (mm).
FACE_PLATE_PILLOW_MIN_RAD = 5;
// Pillow inset at the face surface (mm).
FACE_PLATE_PILLOW_INSET = 2;
// Additional pillow inset at the raised height (mm).
FACE_PLATE_PILLOW_RAMP = 1;
// Pillow thickness above the face (mm).
FACE_PLATE_PILLOW_THK = 0.4;

// Clearance height for face sockets (mm) - make larger if face is far inset into a face.
FACE_PLATE_CLEAR_HEIGHT = 10;

/**
 * Module: Emit a face plate clipped by the current face's anti-interference volume.
 * Params: face_thk (plate thickness), idx/face_pts3d_local/poly_faces_idx/poly_verts_local/face_neighbors_idx/face_dihedrals (optional overrides; default from `place_on_faces` context), clear_space (emit clearance cutter), pillow_* (raised pillow sizing), base_z (bottom Z; defaults to `-face_thk` so the top sits on the source face plane), clear_height (clearance height), mode/max_project/eps/convexity (anti-interference controls)
 * Returns: none
 * Limitations/Gotchas: requires `place_on_faces` context or explicit context overrides; the pillow intentionally embosses the source face loop rather than per-shell cut-through loops; true edge insets should be applied by subtractive edge operations around this body
 */
module face_plate(face_thk,
    idx = $ps_face_idx,
    face_pts3d_local = $ps_face_pts3d_local,
    poly_faces_idx = $ps_poly_faces_idx,
    poly_verts_local = $ps_poly_verts_local,
    face_neighbors_idx = $ps_face_neighbors_idx,
    face_dihedrals = $ps_face_dihedrals,
    clear_space=false,
    pillow_min_rad = FACE_PLATE_PILLOW_MIN_RAD,
    pillow_inset = FACE_PLATE_PILLOW_INSET,
    pillow_ramp = FACE_PLATE_PILLOW_RAMP,
    pillow_thk = FACE_PLATE_PILLOW_THK,
    base_z = undef,
    clear_height = FACE_PLATE_CLEAR_HEIGHT,
    mode = "nonzero",
    max_project = undef,
    eps = 1e-4,
    convexity = 6
) {
    assert(!is_undef(idx), "face_plate: idx requires place_on_faces context or an explicit override");
    assert(!is_undef(face_pts3d_local), "face_plate: face_pts3d_local requires place_on_faces context or an explicit override");
    assert(!is_undef(poly_faces_idx), "face_plate: poly_faces_idx requires place_on_faces context or an explicit override");
    assert(!is_undef(poly_verts_local), "face_plate: poly_verts_local requires place_on_faces context or an explicit override");
    assert(!is_undef(face_neighbors_idx), "face_plate: face_neighbors_idx requires place_on_faces context or an explicit override");
    assert(!is_undef(face_dihedrals), "face_plate: face_dihedrals requires place_on_faces context or an explicit override");

    base_z_eff = is_undef(base_z) ? -face_thk : base_z;
    top_z = base_z_eff + face_thk;
    pts = ps_xy(face_pts3d_local);
    shells = ps_face_anti_interference_shells(
        face_pts3d_local,
        idx,
        poly_faces_idx,
        poly_verts_local,
        face_neighbors_idx,
        face_dihedrals,
        base_z_eff,
        top_z,
        mode,
        max_project,
        eps
    );

    color(len(pts) == 3 ? "white" : "red") {
        union() {
            for (shell = shells) {
                if (shell[3] > 0)
                    echo(str("face_plate: capped ", shell[3], " projection(s) on face ", idx, " loop ", shell[2]));

                polyhedron(points = shell[0], faces = shell[1], convexity = convexity);
            }
        }

        // The pillow is source-face decoration: punch-through cuts are incidental
        // and should be handled by later proxy subtraction, not by fragmenting the
        // embossing over anti-interference shell loops.
        loop_centroid = [
            sum([for (p = pts) p[0]]) / len(pts),
            sum([for (p = pts) p[1]]) / len(pts)
        ];
        loop_rad = sum([for (p = pts) norm([p[0] - loop_centroid[0], p[1] - loop_centroid[1]])]) / len(pts);
        if (loop_rad > pillow_min_rad && loop_rad > pillow_inset + pillow_ramp + eps) {
            s0 = max(0, 1 - pillow_inset / loop_rad);
            s1 = max(0, 1 - (pillow_inset + pillow_ramp) / loop_rad);
            p0 = [for (p = pts) [
                loop_centroid[0] + (p[0] - loop_centroid[0]) * s0,
                loop_centroid[1] + (p[1] - loop_centroid[1]) * s0
            ]];
            scale_xy = s0 <= eps ? 0 : (s1 / s0);
            translate([0, 0, top_z])
                linear_extrude(height = pillow_thk, scale = scale_xy)
                    ps_polygon(points = p0, mode = "nonzero");
        }
    }

    if (clear_space)
        color("magenta")
            translate([0, 0, top_z - eps])
                linear_extrude(height = clear_height)
                    ps_polygon(points = pts, mode = "nonzero");
}

// Direct smoke demo: subtract one placed star face cutter from a cube.
_demo_poly = poly_antiprism(n=5, p=2, angle=15);
difference() {
    translate([0, -15, -15]) cube(30);

    place_on_faces(_demo_poly, 17)
        if ($ps_face_idx == 1)
            #face_plate(face_thk=1.2, clear_space=false, max_project=10);
}
