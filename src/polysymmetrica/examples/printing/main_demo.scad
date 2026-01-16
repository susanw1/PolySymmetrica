use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <edge_seg.scad>
use <face_plate.scad>

IR = 20;

//p = (tetrahedron());
//p = poly_truncate(octahedron());
//p = (dodecahedron());
p = poly_truncate(dodecahedron());

EDGE_T = 3;
FACE_T = 0.8;
FIN_T = 0.8;

SHOW_FACES = undef; //[0:5]; // [0,1]; // [0:50];
SHOW_EDGES = true;

difference() {
    if (is_undef(SHOW_EDGES) || SHOW_EDGES) {
        color("gray")
        place_on_edges(p, IR) {
            edge_seg($ps_edge_len, $ps_poly_center_local, $ps_edge_midradius);
        }
    } else {
        group();
    }
    
    place_on_faces(p, IR) {
        if (is_undef(SHOW_FACES) || len(search($ps_facet_idx, [for (i=SHOW_FACES) i])) > 0) {
            face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_facet_dihedrals, undef, is_undef(SHOW_FACES), edge_inset = 1.5);
        }
    }

    sphere(r = IR * 5/6);
}

