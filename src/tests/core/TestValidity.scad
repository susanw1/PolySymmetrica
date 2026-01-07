use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/regular_all.scad>
use <../testing_util.scad>

EPS = 1e-7;

module assert_true(v, msg="") {
    assert(v == true, str(msg, " expected=true got=", v));
}

module assert_false(v, msg="") {
    assert(v == false, str(msg, " expected=false got=", v));
}

function _star_pyramid_poly() =
    let(
        base = [ for (i = [0:1:4]) [cos(72*i), sin(72*i), 0] ],
        apex = [0, 0, 1],
        verts = concat(base, [apex]),
        face_star = [0, 2, 4, 1, 3],
        faces = [
            face_star,
            [0, 2, 5],
            [2, 4, 5],
            [4, 1, 5],
            [1, 3, 5],
            [3, 0, 5]
        ]
    )
    [verts, faces, 1];

module test_validity__regulars_closed() {
    assert_poly_valid_mode(tetrahedron(), "closed");
    assert_poly_valid_mode(octahedron(), "closed");
    assert_poly_valid_mode(hexahedron(), "closed");
}

module test_validity__star_face_modes() {
    p = _star_pyramid_poly();
    assert_false(poly_valid(p, "closed"), "star face disallowed in closed mode");
    assert_true(poly_valid(p, "star_ok"), "star face allowed in star_ok mode");
}

module test_validity__duplicate_indices() {
    verts = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]];
    faces = [[0,1,2],[0,3,1],[0,2,3],[1,1,2]]; // duplicate index
    p = [verts, faces, 1];
    assert_false(poly_valid(p, "closed"), "duplicate indices should fail");
}

module run_TestValidity() {
    test_validity__regulars_closed();
    test_validity__star_face_modes();
    test_validity__duplicate_indices();
}

run_TestValidity();
