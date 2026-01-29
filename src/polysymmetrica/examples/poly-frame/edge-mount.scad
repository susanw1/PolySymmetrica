use <../../core/funcs.scad>
use <polygon-sym.scad>

// Tools and examples for building various face styles


EDGE_DIAM = 5;
EDGE_KEEL = 5;
EDGE_KEEL_IN = 2;
EDGE_KEEL_T = 2;

$fn = 20;

// edge-component of a polygon to be placed on polygonal faces
module edge_mount(len, edge_diam = EDGE_DIAM) {
    hull() {
        sphere(d = edge_diam);
        translate([len, 0, 0]) sphere(d = edge_diam);
    }
    translate([EDGE_KEEL_IN-1, -EDGE_KEEL_T/2, -EDGE_KEEL])
        cube([len-EDGE_KEEL_IN*2+2, EDGE_KEEL_T, EDGE_KEEL]);
}

// replicates the edge-mount around a polygon
module regular_polygon_frame(n_vertex, edge_len) {
    polygon_rad = ps_calc_radius(n_vertex, edge_len);
    apply_polygon_sym(n_vertex, polygon_rad) {
        edge_mount(edge_len);
    }
}


// simple 2d circle
//module regular_polygon_2d(n_vertex, edge_len) {
//    polygon_rad = ps_calc_radius(n_vertex, edge_len);
//    circle(r = polygon_rad, $fn = n_vertex);
//}


// face with an appealing cushion appearance
module face_cushion(n_vertex, edge_len, top_h = 2.5, ledge_h = 1.2, ledge_w = 4) {
    polygon_rad = ps_calc_radius(n_vertex, edge_len);
    hull() {
        apply_polygon_sym(n_vertex, polygon_rad) {
            cube([0.01, 0.01, ledge_h]);
        }
    }
    hull() {
        apply_polygon_sym(n_vertex, polygon_rad - ledge_w * 2) {
            H = top_h;
            intersection() {
                sphere(r = H, $fn = 25);
                translate([-H, -H, 0]) cube([H * 2, H * 2, H]);
            }
        }
    }
}

regular_polygon_frame(5, 40);

translate([100,0,0]) edge_mount(40);

translate([200,0,0]) face_cushion(6, 40);
