# Proxy Interaction

Proxy interaction is the punch-through layer that replays deliberate proxy
geometry for foreign placement sites. It exists because OpenSCAD cannot inspect
arbitrary geometry after a child module has emitted it.

## Contract

Use `place_on_face_foreign_proxy_sites(...)` inside `place_on_faces(...)`:

```scad
place_on_faces(poly) {
    if ($ps_face_idx == target_face_idx) {
        difference() {
            my_target_face_geometry();

            place_on_face_foreign_proxy_sites() {
                my_foreign_face_proxy();   // child slot 0
                my_foreign_edge_proxy();   // child slot 1, reserved
                my_foreign_vertex_proxy(); // child slot 2, reserved
            }
        }
    }
}
```

The module dispatches to child slots by source kind:

- `face_child=0` receives foreign face proxy callbacks.
- `edge_child=1` receives foreign edge proxy callbacks.
- `vertex_child=2` receives foreign vertex proxy callbacks.

Foreign face callbacks are exact face-plane intrusions. Foreign edge and vertex
callbacks are conservative provenance-driven candidates derived from those face
intrusions: they identify source edges and vertices implicated by a known
punch-through. For each exact intruding face, all of its boundary edges and
vertices are emitted as candidates, then deduplicated by source kind/index.
They are not distance/proximity envelope tests.

## Coordinate Modes

`coords="element"` is the default. Children run in the replayed foreign site
frame, transformed into the current target face-local coordinate system. This is
the normal mode for proxy bodies authored like placement children.

`coords="parent"` leaves children in the target face-local coordinate system.
Use this for debugging or for manually applying `$ps_proxy_center_local` and the
`$ps_proxy_*_local` axes.

## Metadata

Each callback exposes:

- `$ps_proxy_idx`, `$ps_proxy_count`
- `$ps_proxy_kind`
- `$ps_proxy_source_kind`, `$ps_proxy_source_idx`
- `$ps_proxy_target_face_idx`, `$ps_proxy_child_idx`
- `$ps_proxy_center_local`, `$ps_proxy_ex_local`, `$ps_proxy_ey_local`, `$ps_proxy_ez_local`
- `$ps_proxy_intrusion_record`, `$ps_proxy_intrusion_segment2d_local`, `$ps_proxy_intrusion_dihedral`, `$ps_proxy_intrusion_confidence`
- `$ps_proxy_face_pts2d`, `$ps_proxy_face_pts3d_local`, `$ps_proxy_face_verts_idx`
- `$ps_proxy_edge_pts_local`, `$ps_proxy_edge_verts_idx`, `$ps_proxy_edge_adj_faces_idx`
- `$ps_proxy_vertex_valence`, `$ps_proxy_vertex_neighbors_idx`, `$ps_proxy_vertex_neighbor_pts_local`
- `$ps_proxy_poly_verts_local`, `$ps_proxy_poly_center_local`

The corresponding `$ps_replay_*` variables are also exposed for compatibility
with the lower-level replay iterator.

With `coords="element"`, the child also receives the normal placement context
for its source kind:

- face callbacks expose the usual `$ps_face_*` variables;
- edge callbacks expose the usual `$ps_edge_*` variables;
- vertex callbacks expose the usual `$ps_vertex_*` variables.

## Responsibilities

The proxy child must emit the boolean body that represents possible collision
from that foreign source. For reliable later subtraction, that body should be
closed and deliberately oversized/toleranced where appropriate.

The library only identifies and positions replay contexts. It cannot discover
already-rendered arbitrary geometry, and it does not infer a closed punch-through
body from filtered face-plane cut segments.

## Printable Face Plates

`examples/printing/face_plate.scad` provides an opt-in printable wrapper:

```scad
place_on_faces(poly) {
    if ($ps_face_idx == target_face_idx) {
        face_plate_minus_foreign_proxies(
            $ps_face_idx,
            $ps_face_pts2d,
            face_thk,
            $ps_face_dihedrals,
            undef,
            false
        ) {
            my_foreign_face_proxy();   // child slot 0
            my_foreign_edge_proxy();   // child slot 1
            my_foreign_vertex_proxy(); // child slot 2
        }
    }
}
```

This is just:

```scad
difference() {
    face_plate(...);
    place_on_face_foreign_proxy_sites(...) children();
}
```

The wrapper keeps the proxy layer separate from the positive face-body builder:
it subtracts supplied foreign proxy bodies from `face_plate(...)`; it does not
imply automatic clearance of arbitrary geometry.
