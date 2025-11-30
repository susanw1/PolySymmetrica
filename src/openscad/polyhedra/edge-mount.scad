use <funcs.scad>
use <polygon-sym.scad>

EDGE_DIAM = 5;
EDGE_KEEL = 5;
EDGE_KEEL_T = 2;

$fn = 20;


module edge_mount(len, edge_diam = EDGE_DIAM) {
    hull() {
        sphere(d = edge_diam);
        translate([len, 0, 0]) sphere(d = edge_diam);
    }
    translate([EDGE_KEEL-1, -EDGE_KEEL_T/2, -EDGE_KEEL]) 
        cube([len-EDGE_KEEL*2+2, EDGE_KEEL_T, EDGE_KEEL]);
}

module regular_polygon_frame(n_vertex, edge_len) {
    polygon_rad = calc_radius(n_vertex, edge_len);
    apply_polygon_sym(n_vertex, polygon_rad) {
        edge_mount(edge_len);
    } 
}


regular_polygon_frame(5, 40);