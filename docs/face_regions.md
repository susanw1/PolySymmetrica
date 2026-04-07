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
- split cut spans bounded by finite run endpoints, not treated as infinite
  cutter joins
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

There is also a plausible future **proxy interaction** path for fabrication:

- keep the analytic layer for visibility / provenance / feature classification
- but let user-supplied face / edge / vertex proxy solids drive the actual
  interference subtraction through ordinary CSG
- the current signature sketch for that work lives in
  `core/proxy_interaction.scad`

That would be a separate higher-level construction path, not a replacement for
the analytic topology layer described here.

That proxy path should explicitly distinguish:

- **occupancy**
  - the material that really exists in the finished face / edge / vertex
    structure
- **clearance**
  - the empty space intentionally reserved for fitting, seating, or tolerance

The intended rule is:

- foreign intersecting polys should cut each other using **occupancy**
- local seats / inset gaps / skeleton strips should be created using
  **clearance**

So a future proxy workflow should avoid reusing local seat-cutting strips as
the inter-poly cutter. Those are different semantics.

The next proxy step should be built around three layers:

- **shielded face**
  - supplied face geometry, clipped by the face slab and by one dihedral/2
    half-space per adjacent edge
- **foreign occupied volume**
  - non-adjacent intersecting faces / edges / vertices only
- **local clearance**
  - the target face's own seating / inset subtraction

In compact notation:

- `F_shield(i) = F_raw(i) ∩ Z(i) ∩ ⋂[e in boundary(i)] B(i,e)`
- `V(i) = ⋃[x in I(i)] Occ(x)`
- `F_final(i) = F_shield(i) - V(i) - C_local(i)`

The important rule is:

- adjacent faces should not be treated as foreign cutters once the bisector
  shields are active

That proxy model is now the design basis. The older corridor/ownership
experiments should be treated as failed diagnostics, not active design.

One important complication remains for star and other self-intersecting faces:

- a global intersection of edge-derived bisector half-spaces produces the
  convex kernel, not the intended filled face shape

So proxy shielding for those faces must not operate on the original
self-crossing walk directly. Instead, it should:

1. split face edges at self-intersections
2. create pseudo-vertices at crossing points
3. classify the filled planar arrangement using the chosen fill rule
4. keep only the subsegments that form the true filled boundary
5. apply shields/clearance from those boundary subsegments, not from internal
   crossing lines

Because the filled result may still be nonconvex, the shielded face should then
be built as a union of convex atom/cell results rather than one monolithic
half-space intersection.

## What Ended Up Working

The important practical lesson from the recent segmentation work is:

- visible cells from `segments.scad` can still be **nonconvex**
- a half-plane-clipped loop at one `z` is only trustworthy for a convex cell
  or convex atom
- for a nonconvex visible cell, that clipped loop is at best the **convex
  kernel** of the cell, and may even be empty once half-plane orientation is
  handled strictly

That means a direct "replace the live loop with the clipped loop" fix is wrong
for nonconvex cells. It shrinks the cell to its kernel and destroys real
visible area.

The working remedy lives in `face_regions.scad`:

1. keep `segments.scad` as the analysis layer
2. detect nonconvex visible cells in `face_regions.scad`
3. decompose those cells into convex atoms
4. build one admissible region per convex atom
5. union the atom regions back together

This keeps the example layer thin, preserves nonconvex visible pieces, and
still lets the region builder operate on convex pieces locally.

It also explains why the bad `poly_antiprism(7,3, angle=15)` `f2/c0` case was
hard:

- the problematic blue visible cell was nonconvex
- the orange live loop was too large because adjacent-line intersection was too
  loose
- the cyan clipped loop looked "right" locally because it was the convex kernel
  near the bad corner
- but using that cyan loop directly for the whole cell broke the rest of the
  piece

So the current live strategy is:

- build segmented visible-cell regions from the full visible cell shape
- but decompose nonconvex cells into convex atoms before lofting

One important correction from the later debugging passes:

- `z0` and `z1` are only bounds on the admissible region
- they do **not** imply different geometric rules at the "top" and "bottom"

In core region logic, a parent edge or cut edge should contribute one side
plane across the full `z0..z1` span. If a cut span is finite, its start and end
should likewise be bounded by full-depth run-end planes. Any tapering,
styling, pillow/roof behaviour, or cap-specific shaping belongs in the
consumer, not in `face_regions.scad`.

That remains the target model. The current live path still carries `band_z0` /
`band_z1` compatibility plumbing, and recent attempts to remove that asymmetry
directly caused broader region-builder regressions. So treat "full-depth
planes across `z0..z1`" as the design direction, not yet a completed live-path
fact.

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
- optionally `$ps_vis_seg_cut_pair_ids`
- optionally `$ps_vis_seg_cut_run_ids`

### `ps_face_visible_cell_region_volume(cell, face_diheds, cut_entries, poly_faces_idx, poly_verts_local, z0, z1, band_z0=undef, band_z1=undef, cut_clearance=0, eps=1e-8)`

Builds the segmented visible-cell region used by the live clipping path.

- Parent edges use the true face dihedral/2 rule.
- Cut edges use cutter-derived join geometry.
- Convex cells are lofted directly.
- Nonconvex cells are first decomposed into convex atoms, and those atom
  regions are then unioned.

That convex-atom decomposition is the key to handling punch-through cases and
future stellation-style intersections without collapsing the visible region to a
convex kernel.

### `ps_face_visible_cell_region_volume_ctx(...)`

Context wrapper for the direct visible-cell region builder.

### `ps_face_visible_cell_region_planes(...)`

Returns the exact side-plane set for one visible cell:

- bottom slab plane
- top slab plane
- one kept plane per parent edge
- one kept plane per cut edge

This is primarily a debug/probe helper for cross-section reasoning and future
exact-region work.

With `include_run_ends=true`, it can also append full-depth run-end planes for
split cut spans. That remains an analysis/debug path for now; the live builder
does not yet consume those extra planes.

### `ps_face_visible_cell_cut_run_end_entries(...)`

Returns the finite run-end boundaries implied by `cut_run_id` on one visible
cell boundary.

Each entry describes one full-depth endpoint plane for a split cut span:

- local 2D line/plane data
- source edge index on the visible cell
- whether it is the run start or run end
- the local `cut_run_id`

This is currently an analysis/debug helper. It exists so future endpoint
handling can be keyed to finite cut spans rather than treating one cutter face
as an infinite join across the whole cell.

### `ps_clip_to_visible_face_cell_ctx(...)`

Child-consuming wrapper for one visible segmented face cell.

- Baseline behavior: clip to the simple visible-cell volume only.
- With `apply_cut_bands=true`, use the direct segmented visible-cell region
  instead.

This direct region path is the current preferred segmented path.

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
