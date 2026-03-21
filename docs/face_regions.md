# Face Region Volumes

This document is the companion to [segments.md](segments.md).

- [segments.md](segments.md) explains how to split faces, derive cut segments,
  and enumerate visible cells.
- `face_regions.scad` takes that data and turns it into 3D keep/cut volumes for
  arbitrary child geometry.

If you are deciding where to start:

- start with `segments.scad` if you need **data**
- start with `face_regions.scad` if you need **clipping volumes**

`src/polysymmetrica/core/face_regions.scad` contains local face-frame volume helpers for:

- the admissible region around a real poly face
- visible-cell keep volumes for segmented faces
- segmentation-edge join regions derived from intersecting faces

These helpers are deliberately local. They depend on:

- the face polygon in face-local 2D
- the real face-edge dihedrals
- cut-edge metadata from `segments.scad`

They do **not** depend on any global polyhedral centre convention.

## Philosophy

The important design decision is that face-region logic must be driven by the
**underlying polyhedral geometry**, not by the detailed form of whatever child
geometry the user places on the face.

That gives a stable contract:

1. `segments.scad` determines the visible face cells and the cut provenance.
2. `face_regions.scad` turns that into admissible 3D regions.
3. User/example geometry is clipped to those regions.

This avoids a circular problem.

It is **not** generally possible to say:

- "take arbitrary decorated geometry on every face, look at all final mutual
  intersections, and automatically infer the perfect result in one pass"

because each final clipped face would depend recursively on the others.

It **is** possible to say:

- "derive the legal keep region from the polyhedron alone, then clip arbitrary
  child geometry to that region"

That is the supported model here.

### What Is Supported

- ordinary face regions bounded by true face-edge dihedral/2 planes
- segmented visible cells derived from intersecting faces
- segmentation joins derived from the underlying cutter faces
- arbitrary child geometry clipped to those admissible regions

### What Is Not Generically Solvable

- arbitrary child geometry that intentionally protrudes outside the admissible
  region and then expects the system to recursively reconcile those extra
  intersections

For example, if a user places a tall looping shape that leaves one face,
crosses a neighboring face, and comes back elsewhere, that becomes a bespoke
constructive-modeling problem rather than a generic face-placement problem.

So the guarantee here is:

- the admissible region is correct and stable

not:

- every arbitrary clipped child shape will look aesthetically simple after
  clipping

That distinction matters for examples such as `face_plate(...)`, which layer
their own body/roof/pillow structure on top of the generic region model.

## Purpose

This file exists to keep face-intersection / segmentation geometry in `core`, rather than growing example-specific logic inside `examples/printing/face_plate.scad`.

The intended layering is:

1. `segments.scad` determines visible cells and cut provenance.
2. `face_regions.scad` turns that local information into keep/cut volumes.
3. Example code such as `face_plate.scad` intersects arbitrary child geometry with those volumes.

That split is deliberate:

- `segments.scad` answers "what are the cells and cuts?"
- `face_regions.scad` answers "what volume is legal to keep?"

## Current Plan

The current direction is deliberately narrow and geometric:

1. Keep `segments.scad` focused on **analysis only**:
   - split loops
   - derive cut entries
   - determine visible cells
   - propagate cut provenance

2. Keep `face_regions.scad` focused on **region construction only**:
   - ordinary face regions from parent-edge dihedral/2 planes
   - segmented visible-cell regions from all cell side lines together

3. For segmented cells, prefer a **direct region build** over "keep prism minus
   a union of independent cut bands".
   This matters when multiple cut edges meet in one place, for example when a
   true poly edge punches through another face.

4. Keep examples such as `face_plate.scad` thin:
   - build example geometry
   - clip it to the core face/visible-cell regions
   - do not solve generic segmentation inside the example layer

This is the route intended to support more difficult future cases such as
stellations.

## Main Helpers

### `ps_face_region_volume(face_pts2d, face_diheds, z0, z1, mode="nonzero", eps=1e-8)`

Builds the local face-region volume between `z0` and `z1`.

- Real face edges are split by the usual dihedral/2 rule.
- Positive `z` is outward from the face plane.
- `mode` is passed through to `ps_face_segments(...)`, so star/self-intersecting faces can use `"nonzero"` (default solid fill) or `"evenodd"` semantics consistently.

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

Typical use:

```scad
place_on_faces(p, edge_len = 20) {
    ps_clip_to_face_region_ctx(-2, 4) {
        children();
    }
}
```

That means "allow arbitrary geometry between `z=-2` and `z=4`, clipped by the
real face-edge dihedral planes."

### `ps_face_visible_cell_volume(cell, z0, z1, cut_clearance=0, eps=1e-8)`

Builds a simple keep volume for one visible face cell returned by `ps_face_visible_segments(...)`.

- Parent edges stay on the cell boundary.
- Cut edges are pulled inward by `cut_clearance`.

This is the simple segmented baseline used when you only want vertical
cut-clearance handling.

The `cut_clearance` here is the nominal segmentation gap. It is independent of
ordinary outer face-edge beveling.

### `ps_face_visible_cell_volume_ctx(z0, z1, cut_clearance=0, eps=1e-8)`

Context wrapper for use inside `place_on_face_visible_segments(...)`.

Requires:

- `$ps_vis_seg_pts2d`
- `$ps_vis_seg_edge_kinds`
- optionally `$ps_vis_seg_cut_entry_ids`

### `ps_face_visible_cell_region_volume(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8)`

Builds a segmented visible-cell region directly from **all** of the cell's side
constraints together.

- Parent edges use the true face dihedral/2 rule.
- Cut edges use cutter-derived join geometry.
- The loop at each `z` is formed by intersecting adjacent side lines in order.

This direct region construction is the important model for complex local cases
where multiple cut edges meet in one place.

### `ps_face_visible_cell_region_volume_ctx(...)`

Context wrapper for the direct visible-cell region builder.

### `ps_clip_to_visible_face_cell_ctx(...)`

Child-consuming wrapper for one visible segmented face cell.

- Baseline behavior: clip to the simple visible-cell volume only.
- With `apply_cut_bands=true`, use the direct segmented visible-cell region
  instead.

Typical use:

```scad
place_on_faces(p, edge_len = 20)
    place_on_face_visible_segments(mode = "nonzero") {
        ps_clip_to_visible_face_cell_ctx(
            z0 = -2,
            z1 = 4,
            cut_clearance = 0.5
        ) {
            children();
        }
    }
```

### `ps_clip_to_visible_face_segments_ctx(...)`

Child-consuming wrapper that unions the clipped result over all visible cells
of the current face.

This is the generic segmented-face consumer now used by
`examples/printing/face_plate.scad`.

Consumers can keep the simple baseline by leaving `apply_cut_bands=false`, or
use the direct segmented visible-cell region by passing `apply_cut_bands=true`.

Despite the name, the preferred main path is now the direct visible-cell region
construction, not a literal "subtract a union of independent cut bands" model.

## Recommended Usage Patterns

### 1. Clip arbitrary face-mounted geometry to a real face region

```scad
place_on_faces(poly_cupola(5), edge_len = 20) {
    ps_clip_to_face_region_ctx(z0 = -1.5, z1 = 3) {
        union() {
            sphere(3, $fn = 48);
            translate([0, 0, 1]) cylinder(h = 4, r = 1.2);
        }
    }
}
```

### 2. Clip arbitrary geometry to the visible segmented pieces of a face

```scad
place_on_faces(poly_prism(n=7, p=3), edge_len = 20) {
    if ($ps_face_idx == 2)
        ps_clip_to_visible_face_segments_ctx(
            z0 = -2,
            z1 = 4,
            cut_clearance = 0.6,
            mode = "nonzero",
            apply_cut_bands = true
        ) {
            linear_extrude(height = 1, center = true)
                offset(delta = 0.4)
                    polygon(points = $ps_face_pts2d);
        }
}
```

This is the generic version of what `face_plate_visible(...)` now does.

### 3. Keep the separation straight

- Want split cells or cut lines as arrays? Use [segments.md](segments.md).
- Want to clip arbitrary 3D geometry to those cells? Use `face_regions.scad`.
- Want a quick strip-style subtraction helper for a demo? `face_cut_stencil(...)`
  in `segments.scad` still exists, but it is not the preferred generic path.

## Notes

- The current printing path is built on the core wrappers in this file. Example
  code should stay a thin consumer of these region builders.
- Treat `cut_clearance` as the segmentation join-gap control.
- The old edge-band `along_pad` tuning belongs to the retired subtraction-band
  approach and is no longer part of the active segmented-region path.
