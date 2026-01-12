use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>
use <../../core/render.scad>

use <../../models/regular_all.scad>


t = undef;
IR = 20;

LAYER1 = -100;

//p = (dodecahedron());
p = poly_truncate(dodecahedron());

EDGE_T = 3.5;
FACE_T = 0.8;
FIN_T = 0.8;

SINGLE_FACE = 20;

module edge_seg(len, pc, r_poly) {
    translate([0, 0, 0]) {
        // builds the edge bar itself
        hull() {
            translate([len/2, 0, 0]) sphere(d = EDGE_T, $fn = 40);
            translate([-len/2, 0, 0]) sphere(d = EDGE_T, $fn = 40);
        }
        
        // creates the support fin under the edge
        hull() {
            translate([len/2,0,0]) cube(FIN_T, center = true);
            translate([-len/2,0,0]) cube(FIN_T, center = true);
            translate(pc) cube(FIN_T, center = true);
        }
    }    
} 

//    color(len(pts) == 3? "white" : "red") { 

module face_plate(idx, pts, ht, pc, clear_space) {
    
    module p(h) {
        linear_extrude(height = h) polygon(points = pts);
    }
    
    // facet edges are bevelled underneath by intersecting with a pyramid towards the poly-centre 
    // to prevent over-erosion of the frame. Might need to be cleverer for facet planes that aren't normal to the poly-radius. The bevelling should really worry about dihedral angles.
    render() {
        intersection() {
            translate([0, 0, -ht/2]) p(ht);
            hull() {
                translate([0, 0, ht/2]) p(0.01);
                translate(pc) cube(FIN_T, center = true);
            }
        }
        // Clear the airspace above the polygon.
        if (clear_space) {
            translate([0, 0, ht/2]) p(EDGE_T);
        }
    }
}

difference() {
    if (is_undef(SINGLE_FACE)) {
        color("gray")
        place_on_edges(p, IR) {
            edge_seg($ps_edge_len, $ps_poly_center_local, $ps_edge_midradius);
        }
    }
    
    place_on_faces(p, IR) {
        echo($ps_vertex_count, $ps_facet_radius, $ps_face_pts2d);

       !if (is_undef(SINGLE_FACE) || $ps_facet_idx == SINGLE_FACE)
       face_plate($ps_facet_idx, $ps_face_pts2d, FACE_T, $ps_poly_center_local, is_undef(SINGLE_FACE));
    }

//    sphere(r = IR * 3/4);
}

*edge_seg(10, [0,0,-20], 20);
