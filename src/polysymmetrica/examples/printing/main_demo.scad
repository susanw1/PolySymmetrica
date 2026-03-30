use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/prisms.scad>
use <../../core/construction.scad>

use <../../models/platonics_all.scad>
use <../../models/archimedians_all.scad>
use <../../models/johnsons_all.scad>

use <edge_seg.scad>
use <face_plate.scad>

SC = 1;
IR = 20 * SC;

//p = tetrahedron();
//p = poly_truncate(dodecahedron());
//p = dodecahedron();
//base = icosahedron();
//sol = solve_cantitruncate_trig(base);
//s = poly_cantellate_norm(base, 0.5);
//p = poly_dual(great_rhombicuboctahedron());
//p = poly_attach(octahedron(), icosahedron(), f1=[0,7]);
//p = poly_attach(octahedron(), icosahedron(), f1=[0,1,2,3,4,5,6,7]);
//p = poly_attach(p1, icosahedron(), f1=0);

//AP_CANT_ROWS = [
//    ["face", "family", 0, ["df", 0.30]],
//    ["face", "family", 1, ["df", 0.30]]
//];

//p = poly_cantellate(poly_antiprism(6), params_overrides=AP_CANT_ROWS);
//p = poly_prism(5);
//p = poly_antiprism(5);
//p = poly_prism(n=5, p=2);
//p = poly_antiprism(n=5, p=2, angle = 0);
//p = poly_antiprism(n=7, p=2, angle = 0);
p = poly_antiprism(n=7, p=3, angle = 15);

//p = j1_square_pyramid();
//p = poly_dual(j2_pentagonal_pyramid());

EDGE_T = 3.5 * SC;
FACE_T = 1.6 * SC;
FIN_T = 1 * SC;
INSET = 1.1 * SC;
FACET_BASE_T = 1;
FACET_BASE_W = 2.2;
BASE_Z = -FACE_T / 4;

SHOW_FACES = undef;
CLEAR_AIRSPACE = true;

// Experimental scale/profile:
//IR = 12 * SC;
//EDGE_T = 3.0 * SC;
//FACE_T = 1.2 * SC;

/**
* Generate skeletal frame:
    show_faces = undef;

* Generate single face, eg face#8:
    show_faces = [8];

* Generate multiple faces, eg face#0,#9, and #23:
    show_faces = [0, 9, 23];
*/
module model(show_faces = undef, clear_airspace = true) {
    let(inter_radius = IR)
    union() {
        *union() {
            color("gray")
            place_on_edges(p, IR) {
                edge_seg($ps_edge_pts_local, $ps_poly_center_local, edge_t = EDGE_T);
            }

            color("blue")
            place_on_faces(p, IR) {
                translate([0, 0, BASE_Z - FACET_BASE_T])
                    linear_extrude(FACET_BASE_T)
                        difference() {
                            ps_polygon(points = $ps_face_pts2d);
                            offset(-FACET_BASE_W)
                                ps_polygon(points = $ps_face_pts2d);
                        }
            }
        }

        place_on_faces(p, IR) {
            color($ps_face_idx < 2 ? "red" : ["yellow", "green", "blue", "white", "orange"][$ps_face_idx % 5], 1)
            if (is_undef(show_faces) || len(search($ps_face_idx, [for (i = show_faces) i])) > 0) {
                face_plate_visible(
                    $ps_face_idx,
                    $ps_face_pts2d,
                    FACE_T,
                    $ps_face_dihedrals,
                    undef,
                    clear_airspace,
                    edge_inset = INSET,
                    base_z = BASE_Z,
                    clear_height = 0.6,
                    seg_apply_cut_bands = true
                );
            }
        }
    }
}

model(SHOW_FACES, CLEAR_AIRSPACE);
