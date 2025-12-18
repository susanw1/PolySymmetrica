use <../../core/placement.scad>
use <../../core/duals.scad>
use <../../core/truncation.scad>

use <../../models/tetrahedron.scad>
use <../../models/octahedron.scad>
use <../../models/icosahedron.scad>

use <../../core/render.scad>


//p = poly_truncate(octahedron(), 0.3);

//p = poly_truncate(icosahedron(), 0.1);
//p1 = poly_truncate(p, 0.3);

p = poly_truncate(tetrahedron(), 0.25);

place_on_faces_ir(p, 30) {
    let (col = $ps_vertex_count == 3? "white" : "blue") {
        color(col)
            polygon(points = $ps_face_pts2d);
        echo($ps_face_pts2d);
    }
}

//render_poly(p, 30);