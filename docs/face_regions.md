# Face Region Volumes

`src/polysymmetrica/core/face_regions.scad` contains local face-frame volume helpers for:

- the admissible region around a real poly face
- visible-cell keep volumes for segmented faces
- bounded cut-band volumes derived from segmentation edges

These helpers are deliberately local. They depend on:

- the face polygon in face-local 2D
- the real face-edge dihedrals
- cut-edge metadata from `segments.scad`

They do **not** depend on any global polyhedral centre convention.

## Purpose

This file exists to keep face-intersection / segmentation geometry in `core`, rather than growing example-specific logic inside `examples/printing/face_plate.scad`.

The intended layering is:

1. `segments.scad` determines visible cells and cut provenance.
2. `face_regions.scad` turns that local information into keep/cut volumes.
3. Example code such as `face_plate.scad` intersects arbitrary child geometry with those volumes.

## Main Helpers

### `ps_face_region_volume(face_pts2d, face_diheds, z0, z1, mode="nonzero", eps=1e-8)`

Builds the local face-region volume between `z0` and `z1`.

- Real face edges are split by the usual dihedral/2 rule.
- Positive `z` is outward from the face plane.
- `mode` is passed through to `ps_face_segments(...)`, so star/self-intersecting faces can use `"evenodd"` or `"nonzero"` semantics consistently.

### `ps_face_region_volume_ctx(z0, z1, mode="nonzero", eps=1e-8)`

Context wrapper for use inside `place_on_faces(...)`.

Requires:

- `$ps_face_pts2d`
- `$ps_face_dihedrals`

### `ps_clip_to_face_region(...)` / `ps_clip_to_face_region_ctx(...)`

Child-consuming wrappers that intersect arbitrary geometry with the ordinary
face-region volume.

These are the generic counterparts to example-specific code such as
`face_plate(...)`.

### `ps_face_visible_cell_volume(cell, z0, z1, cut_clearance=0, eps=1e-8)`

Builds a simple keep volume for one visible face cell returned by `ps_face_visible_segments(...)`.

- Parent edges stay on the cell boundary.
- Cut edges are pulled inward by `cut_clearance`.

This is the current stable baseline used by `face_plate_visible(...)`.

### `ps_face_visible_cell_volume_ctx(z0, z1, cut_clearance=0, eps=1e-8)`

Context wrapper for use inside `place_on_face_visible_segments(...)`.

Requires:

- `$ps_vis_seg_pts2d`
- `$ps_vis_seg_edge_kinds`
- optionally `$ps_vis_seg_cut_entry_ids`

### `ps_clip_to_visible_face_cell_ctx(...)`

Child-consuming wrapper for one visible segmented face cell.

- Baseline behavior: clip to the visible cell volume only.
- Optional `apply_cut_bands=true` enables subtraction of bounded cut-band
  volumes derived from segmentation edges.

### `ps_clip_to_visible_face_segments_ctx(...)`

Child-consuming wrapper that unions the clipped result over all visible cells
 of the current face.

This is the generic segmented-face consumer now used by
`examples/printing/face_plate.scad`.

Consumers can keep the stable baseline by leaving `apply_cut_bands=false`, or
opt into experimental segmentation-edge join relief by passing
`apply_cut_bands=true`.

When `apply_cut_bands=true`, the visible-cell keep volume no longer applies its
own cut-edge clearance. The cut-band subtraction owns cut-edge relief
completely in that mode.

### `ps_face_cut_band_volume_profiled(seg2d, inward_n, profile2d, along_pad=0, eps=1e-8)`

Builds a bounded cut-band volume from an explicit cross-section profile.

- `profile2d` is a list of `[u, z]` points
- `u` is inward/outward offset from the segmentation edge in the face plane
- `z` is face-local height

This is the route used by the current cut-band path.

### `ps_face_visible_segment_cut_bands_ctx(z0, z1, cut_clearance=0, along_pad=0, mode="nonzero", eps=1e-8, band_overcut=1e-3)`

Context wrapper that emits cut-band volumes for the current visible cell.

Requires:

- `place_on_faces(...)` context
- `place_on_face_visible_segments(...)` context

## Notes

- The current printing path uses the cut-band path in core.
- Treat `cut_clearance` as the segmentation join-gap control.
- Treat `along_pad` as along-edge overreach only, not gap width.
- Treat `band_overcut` as a robustness epsilon for exact-coincidence cleanup, not as a primary design parameter.
