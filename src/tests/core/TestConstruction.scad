use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/construction.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/johnsons_all.scad>

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_true(v, msg="") {
    assert(v == true, str(msg, " expected=true got=", v));
}

module assert_false(v, msg="") {
    assert(v == false, str(msg, " expected=false got=", v));
}

module test_poly_boundary_loops__cube_missing_one_face() {
    p = poly_delete_faces(hexahedron(), 0, cap=false, cleanup=false);
    loops = poly_boundary_loops(p);
    assert_int_eq(len(loops), 1, "cube delete one face boundary loop count");
    assert_int_eq(len(loops[0]), 4, "cube delete one face loop arity");
}

module test_poly_delete_faces__cube_missing_one_face_open() {
    p = poly_delete_faces(hexahedron(), 0, cap=false, cleanup=false);
    assert_false(poly_valid(p, "closed"), "cube with deleted face should be open");
    assert_true(poly_valid(p, "struct"), "cube with deleted face should remain structurally valid");
    assert_int_eq(len(poly_faces(p)), 5, "cube delete one face face count");
    assert_int_eq(len(poly_edges(p)), 12, "cube delete one face edge count");
}

module test_poly_cap_loops__cube_restores_closed_mesh() {
    p0 = poly_delete_faces(hexahedron(), 0, cap=false, cleanup=false);
    q = poly_cap_loops(p0);
    assert_true(poly_valid(q, "closed"), "cap loops should restore cube closure");
    assert_int_eq(len(poly_faces(q)), 6, "cap loops restored face count");
    assert_int_eq(len(poly_edges(q)), 12, "cap loops restored edge count");
}

module test_poly_delete_faces__cube_cap_true_restores_closed_mesh() {
    q = poly_delete_faces(hexahedron(), 0, cap=true);
    assert_true(poly_valid(q, "closed"), "delete face with cap should be closed");
    assert_int_eq(len(poly_faces(q)), 6, "delete face with cap face count");
    assert_int_eq(len(poly_edges(q)), 12, "delete face with cap edge count");
}

module test_poly_boundary_loops__cube_missing_two_faces() {
    p = poly_delete_faces(hexahedron(), [0, 1], cap=false, cleanup=false);
    loops = poly_boundary_loops(p);
    assert_int_eq(len(loops), 2, "cube delete two faces boundary loop count");
    assert_int_eq(len(loops[0]), 4, "cube delete two faces loop0 arity");
    assert_int_eq(len(loops[1]), 4, "cube delete two faces loop1 arity");
}

module test_poly_slice__cube_midplane_open() {
    p = poly_slice(hexahedron(), [0,0,0], [0,0,1], keep="above", cap=false);
    loops = poly_boundary_loops(p);
    assert_false(poly_valid(p, "closed"), "cube mid-slice open should not be closed");
    assert_true(poly_valid(p, "struct"), "cube mid-slice open should remain structurally valid");
    assert_int_eq(len(loops), 1, "cube mid-slice open boundary loop count");
    assert_int_eq(len(loops[0]), 4, "cube mid-slice open loop arity");
    assert_int_eq(len(poly_faces(p)), 5, "cube mid-slice open face count");
}

module test_poly_slice__cube_midplane_capped() {
    p = poly_slice(hexahedron(), [0,0,0], [0,0,1], keep="above", cap=true);
    assert_true(poly_valid(p, "closed"), "cube mid-slice capped should be closed");
    assert_int_eq(len(poly_faces(p)), 6, "cube mid-slice capped face count");
    assert_int_eq(len(poly_edges(p)), 12, "cube mid-slice capped edge count");
}

module test_poly_pyramid__square_counts_match_j1() {
    p = poly_pyramid(4);
    q = j1_square_pyramid();
    assert_true(poly_valid(p, "closed"), "square pyramid should be closed");
    assert_int_eq(len(poly_verts(p)), 5, "square pyramid vertex count");
    assert_int_eq(len(poly_faces(p)), 5, "square pyramid face count");
    assert_int_eq(len(poly_edges(p)), 8, "square pyramid edge count");
    assert_int_eq(len(poly_verts(q)), 5, "J1 wrapper vertex count");
    assert_int_eq(len(poly_faces(q)), 5, "J1 wrapper face count");
    assert_int_eq(len(poly_edges(q)), 8, "J1 wrapper edge count");
}

module test_poly_pyramid__pentagonal_counts_match_j2() {
    p = poly_pyramid(5);
    q = j2_pentagonal_pyramid();
    assert_true(poly_valid(p, "closed"), "pentagonal pyramid should be closed");
    assert_int_eq(len(poly_verts(p)), 6, "pentagonal pyramid vertex count");
    assert_int_eq(len(poly_faces(p)), 6, "pentagonal pyramid face count");
    assert_int_eq(len(poly_edges(p)), 10, "pentagonal pyramid edge count");
    assert_int_eq(len(poly_verts(q)), 6, "J2 wrapper vertex count");
    assert_int_eq(len(poly_faces(q)), 6, "J2 wrapper face count");
    assert_int_eq(len(poly_edges(q)), 10, "J2 wrapper edge count");
}

module test_poly_cupola__triangular_counts_match_j3() {
    p = poly_cupola(3);
    q = j3_triangular_cupola();
    assert_true(poly_valid(p, "closed"), "triangular cupola should be closed");
    assert_int_eq(len(poly_verts(p)), 9, "triangular cupola vertex count");
    assert_int_eq(len(poly_faces(p)), 8, "triangular cupola face count");
    assert_int_eq(len(poly_edges(p)), 15, "triangular cupola edge count");
    assert_int_eq(len(poly_verts(q)), 9, "J3 wrapper vertex count");
    assert_int_eq(len(poly_faces(q)), 8, "J3 wrapper face count");
    assert_int_eq(len(poly_edges(q)), 15, "J3 wrapper edge count");
}

module test_poly_cupola__square_counts_match_j4() {
    p = poly_cupola(4);
    q = j4_square_cupola();
    assert_true(poly_valid(p, "closed"), "square cupola should be closed");
    assert_int_eq(len(poly_verts(p)), 12, "square cupola vertex count");
    assert_int_eq(len(poly_faces(p)), 10, "square cupola face count");
    assert_int_eq(len(poly_edges(p)), 20, "square cupola edge count");
    assert_int_eq(len(poly_verts(q)), 12, "J4 wrapper vertex count");
    assert_int_eq(len(poly_faces(q)), 10, "J4 wrapper face count");
    assert_int_eq(len(poly_edges(q)), 20, "J4 wrapper edge count");
}

module test_poly_cupola__pentagonal_counts_match_j5() {
    p = poly_cupola(5);
    q = j5_pentagonal_cupola();
    assert_true(poly_valid(p, "closed"), "pentagonal cupola should be closed");
    assert_int_eq(len(poly_verts(p)), 15, "pentagonal cupola vertex count");
    assert_int_eq(len(poly_faces(p)), 12, "pentagonal cupola face count");
    assert_int_eq(len(poly_edges(p)), 25, "pentagonal cupola edge count");
    assert_int_eq(len(poly_verts(q)), 15, "J5 wrapper vertex count");
    assert_int_eq(len(poly_faces(q)), 12, "J5 wrapper face count");
    assert_int_eq(len(poly_edges(q)), 25, "J5 wrapper edge count");
}

module run_TestConstruction() {
    test_poly_boundary_loops__cube_missing_one_face();
    test_poly_delete_faces__cube_missing_one_face_open();
    test_poly_cap_loops__cube_restores_closed_mesh();
    test_poly_delete_faces__cube_cap_true_restores_closed_mesh();
    test_poly_boundary_loops__cube_missing_two_faces();
    test_poly_slice__cube_midplane_open();
    test_poly_slice__cube_midplane_capped();
    test_poly_pyramid__square_counts_match_j1();
    test_poly_pyramid__pentagonal_counts_match_j2();
    test_poly_cupola__triangular_counts_match_j3();
    test_poly_cupola__square_counts_match_j4();
    test_poly_cupola__pentagonal_counts_match_j5();
}

run_TestConstruction();
