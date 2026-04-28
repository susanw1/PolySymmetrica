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
- `edge_child=1` is reserved for foreign edge proxy callbacks.
- `vertex_child=2` is reserved for foreign vertex proxy callbacks.

Currently only exact foreign face intrusions are discovered. Edge and vertex
slots are part of the public contract so proxy modules can be written against a
stable shape while later candidate APIs are added.

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
- `$ps_proxy_poly_verts_local`, `$ps_proxy_poly_center_local`

The corresponding `$ps_replay_*` variables are also exposed for compatibility
with the lower-level replay iterator.

## Responsibilities

The proxy child must emit the boolean body that represents possible collision
from that foreign source. For reliable later subtraction, that body should be
closed and deliberately oversized/toleranced where appropriate.

The library only identifies and positions replay contexts. It cannot discover
already-rendered arbitrary geometry, and it does not infer a closed punch-through
body from filtered face-plane cut segments.

