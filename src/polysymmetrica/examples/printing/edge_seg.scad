// Edge strut diameter (mm).
EDGE_SEG_EDGE_T = 3;
// Support fin diameter (mm).
EDGE_SEG_FIN_T = 0.8;
EDGE_SEG_SCALE_Y = 2;

module edge_seg(len, pc, edge_t = EDGE_SEG_EDGE_T, fin_t = EDGE_SEG_FIN_T) {
    translate([0, 0, 0]) {
        // builds the edge bar itself
        hull() {
            translate([len/2, 0, 0]) sphere(d = edge_t, $fn = 40);
            translate([-len/2, 0, 0]) sphere(d = edge_t, $fn = 40);
        }
        
        // creates the support fin under the edge - goes to centre of poly
        difference() {
            hull() {
                translate([len/2,0,0]) sphere(d = fin_t, $fn = 10);
                translate([-len/2,0,0]) sphere(d = fin_t, $fn = 10);
                translate(pc) sphere(d = fin_t, $fn = 10);
            }
            translate(pc) sphere(r = norm(pc) * 5/6, $fn = 50);
        }
    }
} 

edge_seg(10, [0,0,-20], 20);
