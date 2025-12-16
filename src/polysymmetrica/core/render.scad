// ---------------------------------------------------------------------------
// PolySymmetrica - Polyhedral Geometry Engine
// Version: 0.1.0
// Copyright 2025 Susan Witts
// SPDX-License-Identifier: MIT

use <funcs.scad>

module render_poly(poly, inter_radius = 1) {
    let(
        k = (inter_radius * poly_e_over_ir(poly)) / poly_unit_edge(poly),
        pts = poly_verts(poly) * k
    )
    polyhedron(
        points = pts,
        faces  = orient_all_faces_outward(pts, poly_faces(poly))
    );
}

