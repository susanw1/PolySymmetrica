use <../../polysymmetrica/core/render.scad>
use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/validate.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/core/prisms.scad>

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module test_ps_render_mesh__auto_keeps_raw_on_convex_cube() {
    p = hexahedron();
    mesh = ps_render_mesh(p, 1, "auto");
    assert_int_eq(len(mesh[1]), len(poly_faces(p)), "cube auto should keep raw faces");
    assert(max([for (f = mesh[1]) len(f)]) > 3, "cube raw faces should remain non-triangulated");
}

module test_ps_render_mesh__always_triangulates_cube() {
    p = hexahedron();
    mesh = ps_render_mesh(p, 1, "always");
    assert_int_eq(len(mesh[1]), 12, "cube triangulated face count");
    assert(min([for (f = mesh[1]) len(f)]) == 3, "cube triangulated min arity");
    assert(max([for (f = mesh[1]) len(f)]) == 3, "cube triangulated max arity");
}

module test_ps_render_mesh__auto_triangulates_star_faces() {
    p = poly_antiprism(5, 2);
    mesh_raw = ps_render_mesh(p, 1, "raw");
    mesh_auto = ps_render_mesh(p, 1, "auto");

    assert_int_eq(len(mesh_raw[1]), len(poly_faces(p)), "star raw face count");
    assert(max([for (f = mesh_raw[1]) len(f)]) > 3, "star raw includes non-triangle faces");

    assert(min([for (f = mesh_auto[1]) len(f)]) == 3, "star auto triangulated min arity");
    assert(max([for (f = mesh_auto[1]) len(f)]) == 3, "star auto triangulated max arity");
    assert(len(mesh_auto[1]) > len(mesh_raw[1]), "star auto should produce more triangles than raw faces");

    q_auto = make_poly(mesh_auto[0], mesh_auto[1]);
    assert(_ps_edges_manifold(mesh_auto[0], mesh_auto[1]), "star auto triangulated mesh should be edge-manifold");
    assert(poly_valid(q_auto, "star_ok"), "star auto triangulated mesh should be star_ok valid");
}

module run_TestRender() {
    test_ps_render_mesh__auto_keeps_raw_on_convex_cube();
    test_ps_render_mesh__always_triangulates_cube();
    test_ps_render_mesh__auto_triangulates_star_faces();
}

run_TestRender();
