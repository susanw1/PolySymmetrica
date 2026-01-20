use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>

use <edge_seg.scad>
use <face_plate.scad>

SC = 1.5;
IR = 20 * SC;

//p = (tetrahedron());
//p = poly_truncate(octahedron());
//p = (dodecahedron());
p = (poly_truncate(dodecahedron()));

EDGE_T = 3.5 * SC;
FACE_T = 1.5 * SC;
FIN_T = 1 * SC;
INSET = 1.2 * SC;

/**
* Generate skeletal frame:
    SHOW_FACES = undef;

* Generate single facet, eg facet#8:
    SHOW_FACES = [8];
    and place '!'

* Generate multiple facets, eg facet#0,#9, and #23:
    SHOW_FACES = [0, 9, 23];
    and place '!'
*/

SHOW_FACES = undef; //[0:5]; // [0,1]; // [0:50];

// This is what we printed, but x1.5, not x2

module model() {
    difference() {
        // Constructs the edge-based frame
        color("gray")
        place_on_edges(p, IR) {
            edge_seg($ps_edge_len, $ps_poly_center_local, edge_t = EDGE_T);
        }

        // Constructs facets, removes them from frame to create facet-fitting sockets.
        place_on_faces(p, IR) {
            // add '!' here to force facets-only:
            if (is_undef(SHOW_FACES) || len(search($ps_facet_idx, [for (i=SHOW_FACES) i])) > 0) {
                face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_facet_dihedrals, undef, is_undef(SHOW_FACES), 
                    edge_inset = INSET, base_z = -FACE_T/2);
            }
        }
    }
}

model();
