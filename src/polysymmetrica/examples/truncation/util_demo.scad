use <../../core/funcs.scad>
use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/render.scad>

IR = 30;
T = 0.01;

COLORS = [ "", "", "", "yellow", "red", "green", "blue", "gray", "red", "white", "red" ];

// Robust face cap:
// - planar faces: single polygon extrusion
// - non-planar faces: triangle-fan fallback to avoid full-loop polygon dropouts
module _face_fill_plate(thk = T) {
    n = len($ps_face_pts2d);
    if (n >= 3) {
        pts3d = is_undef($ps_face_pts3d_local)
            ? [for (p = $ps_face_pts2d) [p[0], p[1], 0]]
            : $ps_face_pts3d_local;
        pts2d = [for (p = pts3d) [p[0], p[1]]];
        is_planar = is_undef($ps_face_is_planar) ? true : $ps_face_is_planar;
        if (is_planar) {
            linear_extrude(height = thk) polygon(points = pts2d);
        } else {
            cx = sum([for (p = pts3d) p[0]]) / n;
            cy = sum([for (p = pts3d) p[1]]) / n;
            c = [cx, cy];
            union() {
                for (i = [0:1:n-1]) {
                    linear_extrude(height = thk)
                        polygon(points = [c, pts2d[i], pts2d[(i+1)%n]]);
                }
            }
        }
    }
}

module demo(p, ir = IR, detail = 0, name = undef) {
    poly_describe(p, detail = detail, name = name);
    place_on_faces(p, ir) {
        let (col = COLORS[$ps_vertex_count]) {
            color(col) {
                _face_fill_plate(T);
//                translate([0,0,3]) text(str($ps_face_idx), halign="center",valign="center", size=4);
            }
        }
    }

    color("silver")
    place_on_edges(p, ir) {
        cube([$ps_edge_len, 1, 1], center = true);
//        translate($ps_poly_center_local) cube([1,1,$ps_edge_midradius], center = false);
//        translate([0,0,2]) text(str($ps_edge_idx), halign="center",valign="center", size=3);
    }

    color("gold")
    place_on_vertices(p, ir) {
        sphere(1.5, $fn=15);
//        translate([0,0,2]) text(str($ps_vertex_idx), halign="center",valign="center", size=2);
    }
}

module combo(p, scale_f = function(p,d) scale_dual_edge_cross(p,d, 0)) {
    d = poly_dual(p);
    m = scale_f(p, d);       
    
    color("blue", 1)
    place_on_faces(p, IR) {
        translate([0,0,-T]) _face_fill_plate(T);
    }
    color("yellow", 1)
    place_on_faces(d, IR * m) {
        translate([0,0,-T]) _face_fill_plate(T);
    }
}

// Turns a list into a list of [i+start, a[i]] tuples.
function with_index(a, start=0) = [ for (i = [0:1:len(a)-1]) [ i + start, a[i] ] ];
