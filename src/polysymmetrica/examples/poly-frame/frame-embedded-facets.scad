use <../../core/placement.scad>
use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>
use <edge-mount.scad>

use <../../core/duals.scad>


// Uncomment a value of 'd' as required, then update N to number of edges on facet

//d = tetrahedron();
//d = octahedron();
d = icosahedron();
//
//d = poly_dual(octahedron());
//d = poly_dual(icosahedron());

difference() {
    N = 3;

    place_on_faces(d, 30) {
        regular_polygon_frame(N, $ps_edge_len);
    }
    place_on_faces(d, 30) {
        // main facet base
        translate([0, 0, 1]) {
            linear_extrude(height = 5)
                circle($fn = N, r = $ps_facet_radius);
            %facet_cushion(N, $ps_edge_len);
        }
    }
}
