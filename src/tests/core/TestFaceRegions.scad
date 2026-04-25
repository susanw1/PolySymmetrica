use <../../polysymmetrica/core/face_regions.scad>
use <../../polysymmetrica/core/funcs.scad>
use <../../polysymmetrica/core/placement.scad>
use <../../polysymmetrica/core/prisms.scad>
use <../../polysymmetrica/core/segments.scad>
use <../../polysymmetrica/core/truncation.scad>
use <../../polysymmetrica/models/platonics_all.scad>
use <../../polysymmetrica/models/tetrahedron.scad>

EPS = 1e-8;

module assert_int_eq(a, b, msg="") {
    assert(a == b, str(msg, " expected=", b, " got=", a));
}

module assert_near(a, b, eps=EPS, msg="") {
    assert(abs(a - b) <= eps, str(msg, " expected=", b, " got=", a));
}

function _test_face_site(poly, face_idx) =
    ps_face_sites(poly)[face_idx];

function _test_shell_points_are_finite(points) =
    len([
        for (p = points)
            if (!is_undef(p) && len(p) == 3 && norm(p) < 1e9)
                1
    ]) == len(points);

module test_ps_face_anti_interference_shells__cube_face_single_quad_shell() {
    p = hexahedron();
    site = _test_face_site(p, 0);
    shells = ps_face_anti_interference_shells(
        site[11],
        site[0],
        site[13],
        site[12],
        site[20],
        site[21],
        -0.4,
        0.6
    );

    assert_int_eq(len(shells), 1, "cube face should produce one filled boundary shell");
    assert_int_eq(len(shells[0][0]), 8, "cube quad shell should have bottom+top vertices");
    assert_int_eq(len(shells[0][1]), 8, "quad shell should have two triangulated caps plus four sides");
    assert(_test_shell_points_are_finite(shells[0][0]), "cube shell points should be finite");

    zs = [for (p = shells[0][0]) p[2]];
    assert_near(min(zs), -0.4, EPS, "cube shell min z");
    assert_near(max(zs), 0.6, EPS, "cube shell max z");
}

module test_ps_face_anti_interference_shells__matches_boundary_loop_count() {
    p = poly_antiprism(5, 2);
    site = _test_face_site(p, 1);
    bm = ps_face_boundary_model(site[11], "nonzero");
    shells = ps_face_anti_interference_shells(
        site[11],
        site[0],
        site[13],
        site[12],
        site[20],
        site[21],
        -0.3,
        0.3,
        "nonzero"
    );

    assert_int_eq(len(shells), len(bm[2]), "one shell per filled boundary loop");
    for (i = [0:1:len(shells)-1])
        assert_int_eq(len(shells[i][0]), 2 * len(bm[2][i][0]), "shell vertices should match loop arity");
}

module test_ps_face_anti_interference_shells__pentagram_zmax_expands_outward() {
    p = poly_antiprism(5, 2);
    site = _test_face_site(p, 1);
    shells = ps_face_anti_interference_shells(
        site[11],
        site[0],
        site[13],
        site[12],
        site[20],
        site[21],
        -0.5,
        0.5,
        "nonzero"
    );

    assert_int_eq(len(shells), 1, "pentagram cap should produce one shell");
    area_zmin = abs(_ps_seg_poly_area2(shells[0][4]));
    area_zmax = abs(_ps_seg_poly_area2(shells[0][5]));
    assert(
        area_zmax > area_zmin,
        str("pentagram +Z cap should expand outward area_zmin=", area_zmin, " area_zmax=", area_zmax)
    );
}

module test_ps_face_anti_interference_shells__anti_tet_hex_is_finite() {
    p = poly_truncate(tetrahedron(), t = -0.5);
    site = _test_face_site(p, 0);
    shells = ps_face_anti_interference_shells(
        site[11],
        site[0],
        site[13],
        site[12],
        site[20],
        site[21],
        -0.25,
        0.25,
        "nonzero",
        20
    );

    assert(len(shells) >= 1, "anti-truncated tetrahedron hex should produce at least one shell");
    for (shell = shells) {
        assert(len(shell[0]) >= 6, "anti-tet shell should have points");
        assert(len(shell[1]) >= 4, "anti-tet shell should have faces");
        assert(_test_shell_points_are_finite(shell[0]), "anti-tet shell points should be finite");
    }
}

module test_ps_face_anti_interference_shells__anti_tet_winding_splits_z_direction() {
    p = poly_truncate(tetrahedron(), t = -0.5);
    site = _test_face_site(p, 0);
    shells = ps_face_anti_interference_shells(
        site[11],
        site[0],
        site[13],
        site[12],
        site[20],
        site[21],
        -0.5,
        0.5,
        "nonzero",
        20
    );

    assert_int_eq(len(shells), 4, "anti-tet hex should split into one centre and three corner shells");

    area_deltas = [
        for (shell = shells)
            abs(_ps_seg_poly_area2(shell[5])) - abs(_ps_seg_poly_area2(shell[4]))
    ];
    expanded_count = sum([for (d = area_deltas) d > EPS ? 1 : 0]);
    shrunk_count = sum([for (d = area_deltas) d < -EPS ? 1 : 0]);

    assert_int_eq(expanded_count, 1, "anti-tet same-winding centre shell should expand toward +Z");
    assert_int_eq(shrunk_count, 3, "anti-tet opposite-winding corner shells should shrink toward +Z");
}

module test_ps_face_anti_interference_projection_cap__limits_offset() {
    assert_near(_ps_fr_project_offset(10, 0.5, 3), 3, EPS, "positive projection cap");
    assert_near(_ps_fr_project_offset(-10, 0.5, 3), -3, EPS, "negative projection cap");
    assert_near(_ps_fr_project_offset(10, 0.5, undef), 20, EPS, "uncapped projection");
    assert_near(_ps_fr_project_offset(10, 0, 3), 3, EPS, "near-flat projection uses cap");
}

module run_TestFaceRegions() {
    test_ps_face_anti_interference_shells__cube_face_single_quad_shell();
    test_ps_face_anti_interference_shells__matches_boundary_loop_count();
    test_ps_face_anti_interference_shells__pentagram_zmax_expands_outward();
    test_ps_face_anti_interference_shells__anti_tet_hex_is_finite();
    test_ps_face_anti_interference_shells__anti_tet_winding_splits_z_direction();
    test_ps_face_anti_interference_projection_cap__limits_offset();
}

run_TestFaceRegions();
