use <../../core/funcs.scad>

/**
Applies polygonal symmetry to the supplied children in the X-Y plane. The child objects are repeated to create an instance per vertex and rotating each so their X-axis points along the edge.
*/
module apply_polygon_sym(n_vertex, polygon_rad) {
    rot_angle = 360 / n_vertex;
    turn_angle = 90 + rot_angle / 2;
    for (i = [0: n_vertex-1]) {
        $pgn_edge_idx = i;
        rotate([0, 0, rot_angle * i]) {
            translate([polygon_rad, 0, 0]) {
                rotate([0, 0, turn_angle]) {
                    children();
                }
            }
        }
    }
}


module regular_polygon_2d(n_vertex, edge_len) {
    polygon_rad = calc_radius(n_vertex, edge_len);
    hull() {
        apply_polygon_sym(n_vertex, polygon_rad) {
            square(0.001);
        }
    }
}

regular_polygon_2d(5, 20);