use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>
use <../../core/classify.scad>
use <../../core/prisms.scad>
use <../../core/construction.scad>
use <../../core/face_regions.scad>

use <../../models/platonics_all.scad>
use <../../models/archimedians_all.scad>
use <../../models/johnsons_all.scad>

use <edge_seg.scad>
use <face_plate.scad>

SC = 1;
IR = 20 * SC;

//p = (tetrahedron());
//p = poly_truncate(dodecahedron());
//p = (dodecahedron());
//base = icosahedron();
//sol = solve_cantitruncate_trig(base);
//s = poly_cantellate_norm(base, 0.5);
//p = poly_dual(great_rhombicuboctahedron());
//p = poly_truncate(poly_dual(poly_truncate(hexahedron())), t=0, params_overrides=[["vert", "id", [0,1,2,3,4,55], ["t",0.5001]]]);
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
p = poly_antiprism(n=5, p=2, angle = 15);

//p = j1_square_pyramid();
//p = poly_dual(j2_pentagonal_pyramid());

EDGE_T = 3.5 * SC; // 3.5
FACE_T = 1.6 * SC; // 1.6 * SC;
INSET = 1.1 * SC;
FACET_BASE_T = 1;
FACET_BASE_W = 2.2;
BASE_Z = -FACE_T / 2;

//// Or, experimental:
//IR = 12 * SC;
//EDGE_T = 3.0 * SC; // 3.5
//FACE_T = 1.2 * SC; // 1.6 * SC;

/**
* Generate skeletal frame:
    show_faces = undef;

* Generate single face, eg face#8:
    show_faces = [8];
    and place '!'

* Generate multiple faces, eg face#0,#9, and #23:
    show_faces = [0, 9, 23];
    and place '!'
*/
module model(show_faces = undef, clear_airspace = true) {
    let (inter_radius = IR)
    difference() {
        union() {
            // Constructs the edge-based frame
            color("gray")
            place_on_edges(p, IR) {
                edge_seg($ps_edge_pts_local, $ps_poly_center_local, edge_t = EDGE_T);
            }

            // mounting plate - the edge frame isn't quite substantial enough
            color("blue")
            place_on_faces(p, IR) {
                translate([0,0, BASE_Z - FACET_BASE_T]) linear_extrude(FACET_BASE_T)
                    difference() {
                        ps_polygon(points = $ps_face_pts2d);
                        offset(-FACET_BASE_W) ps_polygon(points = $ps_face_pts2d);
                    }
            }
        }
        // Constructs faces, removes them from frame to create face-fitting sockets.
        place_on_faces(p, IR) {
            // add '!' here to force faces-only:
            if (is_undef(show_faces) || len(search($ps_face_idx, [for (i=show_faces) i])) > 0) {
                face_plate(face_thk = FACE_T, base_z = BASE_Z, max_project = 10, 
                        clear_space = clear_airspace, clear_height = 0.6);
            }
        }
    }
}

//model();
//poly_render(p, 20);

module demo_face() {
    face_plate(base_z = BASE_Z, face_thk = FACE_T, clear_space = false, clear_height = 0.1, max_project = 10);
}
module demo_edge(endPoints, edge_t) {
    hull() {
        translate(endPoints[0]) sphere(d = edge_t, $fn = 40);
        translate(endPoints[1]) sphere(d = edge_t, $fn = 40);
    }
}

module demo_vert() {
    cylinder(r=3, $fn = $ps_vertex_valence);
}


place_on_faces(p, IR, indices = [2]) {
    union() {
        demo_face();
        place_on_face_foreign_proxy_sites() {
            demo_face();
            demo_edge($ps_edge_pts_local, edge_t = EDGE_T);
            demo_vert();
        }
    }
}


