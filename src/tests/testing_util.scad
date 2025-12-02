use <../polysymmetrica/core/placement.scad>

RAD = 40;

module placement_tester(poly, rad=RAD, face_sides=3) {
    color("yellow") 
    place_on_faces_ir(poly, rad)
        cylinder(r = $ps_facet_radius, h = 2, $fn = face_sides);
        
    color("red") 
    place_on_vertices_ir(poly, rad)
        sphere(5);
        
    color("blue") 
    place_on_edges_ir(poly, rad)
        rotate([0,90,0]) cylinder(r1 = 5, r2 = 1, h = $ps_edge_len/2);
}

